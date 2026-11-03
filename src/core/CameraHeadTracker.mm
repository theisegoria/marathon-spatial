// SPDX-License-Identifier: Apache-2.0
// Where the listener sits, from the front camera: Vision finds the box round each person's upper body, the largest is followed, a
// pinhole model turns its place and size into a position (ms_head_position_from_body), and a little smoothing steadies it. The request
// returns boxes only, so nothing estimates which way anyone faces. Frames are looked at and dropped; nothing is kept.
#import <AVFoundation/AVFoundation.h>
#import <Vision/Vision.h>

#include <atomic>
#include <mutex>

#include "marathon_spatial/marathon_spatial_output.h"

namespace {
constexpr double kStaleSeconds = 0.5;

struct State {
    std::mutex mutex;
    ms_head_position head{};
    double seenAt = -1.0;  // CACurrentMediaTime of the last person seen.
    ms_head_tracker_config config{};
    double lastFrame = 0.0;
};
}  // namespace

@interface MSHeadTrackerDelegate : NSObject <AVCaptureVideoDataOutputSampleBufferDelegate>
@property(nonatomic) State* state;
@end

@implementation MSHeadTrackerDelegate
- (void)captureOutput:(AVCaptureOutput*)output didOutputSampleBuffer:(CMSampleBufferRef)sampleBuffer fromConnection:(AVCaptureConnection*)connection {
    (void)output;
    (void)connection;
    State* s = self.state;
    const double now = CACurrentMediaTime();
    if (now - s->lastFrame < 1.0 / s->config.frames_per_second) return;
    s->lastFrame = now;
    CVPixelBufferRef pixels = CMSampleBufferGetImageBuffer(sampleBuffer);
    if (!pixels) return;
    VNDetectHumanRectanglesRequest* request = [[VNDetectHumanRectanglesRequest alloc] init];
    request.upperBodyOnly = YES;
    VNImageRequestHandler* handler = [[VNImageRequestHandler alloc] initWithCVPixelBuffer:pixels options:@{}];
    if (![handler performRequests:@[ request ] error:nil]) return;
    // The nearest person, by the size of their box. Others in view change nothing.
    VNHumanObservation* largest = nil;
    for (VNHumanObservation* person in request.results)
        if (!largest || person.boundingBox.size.width > largest.boundingBox.size.width) largest = person;
    if (!largest) return;
    const CGRect box = largest.boundingBox;  // Normalised to the image, as the camera sees it.
    ms_head_position reading{};
    if (ms_head_position_from_body(float(CGRectGetMidX(box)), float(box.size.width), s->config.horizontal_fov_degrees, s->config.body_width_m,
                                   &reading) != MS_OK)
        return;
    std::lock_guard lock(s->mutex);
    const float a = s->head.valid && now - s->seenAt < kStaleSeconds ? s->config.smoothing : 1.0f;
    s->head.x_m += a * (reading.x_m - s->head.x_m);
    s->head.z_m += a * (reading.z_m - s->head.z_m);
    s->head.valid = 1;
    s->seenAt = now;
}
@end

struct ms_head_tracker {
    State state;
    AVCaptureSession* session = nil;
    AVCaptureDevice* camera = nil;
    MSHeadTrackerDelegate* delegate = nil;
    dispatch_queue_t queue = nullptr;
    bool configured = false;
};

namespace {
bool configure(ms_head_tracker* t) {
    if (t->configured) return true;
    NSError* error = nil;
    AVCaptureDeviceInput* input = [AVCaptureDeviceInput deviceInputWithDevice:t->camera error:&error];
    if (!input) return false;
    AVCaptureSession* session = [[AVCaptureSession alloc] init];
    session.sessionPreset = AVCaptureSessionPreset640x480;
    if (![session canAddInput:input]) return false;
    [session addInput:input];
    AVCaptureVideoDataOutput* output = [[AVCaptureVideoDataOutput alloc] init];
    output.alwaysDiscardsLateVideoFrames = YES;
    [output setSampleBufferDelegate:t->delegate queue:t->queue];
    if (![session canAddOutput:output]) return false;
    [session addOutput:output];
    t->session = session;
    t->configured = true;
    return true;
}
}  // namespace

extern "C" {

ms_status ms_head_tracker_create(const ms_head_tracker_config* config, ms_head_tracker** tracker) {
    if (!tracker) return MS_ERR_INVALID_ARGUMENT;
    *tracker = nullptr;
    AVCaptureDevice* camera = [AVCaptureDevice defaultDeviceWithMediaType:AVMediaTypeVideo];
    if (!camera) return MS_ERR_DEVICE;
    auto* t = new ms_head_tracker;
    ms_head_tracker_config c = config ? *config : ms_head_tracker_config{};
    if (!(c.horizontal_fov_degrees > 0.0f)) c.horizontal_fov_degrees = 70.0f;
    if (!(c.body_width_m > 0.0f)) c.body_width_m = 0.42f;
    if (!(c.frames_per_second > 0.0f)) c.frames_per_second = 15.0f;
    if (!(c.smoothing > 0.0f && c.smoothing <= 1.0f)) c.smoothing = 0.35f;
    t->state.config = c;
    t->camera = camera;
    t->delegate = [[MSHeadTrackerDelegate alloc] init];
    t->delegate.state = &t->state;
    t->queue = dispatch_queue_create("marathon.spatial.head-tracker", DISPATCH_QUEUE_SERIAL);
    *tracker = t;
    return MS_OK;
}

ms_status ms_head_tracker_start(ms_head_tracker* t) {
    if (!t) return MS_ERR_INVALID_ARGUMENT;
    switch ([AVCaptureDevice authorizationStatusForMediaType:AVMediaTypeVideo]) {
    case AVAuthorizationStatusDenied:
    case AVAuthorizationStatusRestricted: return MS_ERR_UNSUPPORTED;
    case AVAuthorizationStatusNotDetermined:
        // Starts once the user answers, if they allow it.
        [AVCaptureDevice requestAccessForMediaType:AVMediaTypeVideo
                                 completionHandler:^(BOOL granted) {
                                     if (!granted) return;
                                     dispatch_async(t->queue, ^{
                                         if (configure(t)) [t->session startRunning];
                                     });
                                 }];
        return MS_OK;
    default: break;
    }
    dispatch_async(t->queue, ^{
        if (configure(t) && !t->session.running) [t->session startRunning];
    });
    return MS_OK;
}

void ms_head_tracker_stop(ms_head_tracker* t) {
    if (!t) return;
    dispatch_sync(t->queue, ^{
        if (t->session.running) [t->session stopRunning];
    });
    std::lock_guard lock(t->state.mutex);
    t->state.head.valid = 0;
}

void ms_head_tracker_destroy(ms_head_tracker* t) {
    if (!t) return;
    ms_head_tracker_stop(t);
    delete t;
}

void ms_head_tracker_get(const ms_head_tracker* t, ms_head_position* head) {
    if (!head) return;
    *head = ms_head_position{};
    if (!t) return;
    auto& s = const_cast<ms_head_tracker*>(t)->state;
    std::lock_guard lock(s.mutex);
    *head = s.head;
    if (CACurrentMediaTime() - s.seenAt > kStaleSeconds) head->valid = 0;
}

}  // extern "C"
