import CStudioBridge
import Foundation
import StudioKit

/// Where the listener sits, from the front camera through the library's tracker, steering 3D speakers. Position only: the library finds
/// the box round the nearest person's upper body and never estimates which way anyone faces, and others in view change nothing. The
/// camera starts only when asked, and no frame is kept. A position goes to the engine when it has moved a step worth a new design, and
/// once a second regardless, so a rebuilt engine picks it up.
@MainActor
final class CameraHeadTracker {
    nonisolated(unsafe) private var tracker: OpaquePointer?
    private var timer: Timer?
    private var sent: TrackedHeadPosition?
    private var sentValid = false
    private var sinceSend = 0.0
    private(set) var isAvailable: Bool
    private(set) var isRefused = false

    init() {
        var made: OpaquePointer?
        isAvailable = ms_head_tracker_create(nil, &made) == MS_OK
        tracker = made
    }

    deinit { if let tracker { ms_head_tracker_destroy(tracker) } }

    func start(engine: SpatialEngine) {
        guard let tracker, timer == nil else { return }
        if ms_head_tracker_start(tracker) == MS_ERR_UNSUPPORTED { isRefused = true; return }
        isRefused = false
        timer = Timer.scheduledTimer(withTimeInterval: 0.05, repeats: true) { [weak self] _ in
            MainActor.assumeIsolated { self?.poll(engine: engine) }
        }
    }

    func stop(engine: SpatialEngine) {
        timer?.invalidate()
        timer = nil
        if let tracker { ms_head_tracker_stop(tracker) }
        sent = nil
        sentValid = false
        engine.setTrackedHead(nil)
    }

    private func poll(engine: SpatialEngine) {
        guard let tracker else { return }
        var h = ms_head_position()
        ms_head_tracker_get(tracker, &h)
        sinceSend += 0.05
        let head = h.valid != 0 ? TrackedHeadPosition(xMetres: Double(h.x_m), zMetres: Double(h.z_m)) : nil
        let moved: Bool
        switch (head, sent) {
        case let (now?, last?):
            moved = abs(now.xMetres - last.xMetres) >= 0.015 || abs(now.zMetres - last.zMetres) >= 0.015
        case (nil, nil): moved = sentValid
        default: moved = true
        }
        guard moved || sinceSend >= 1 else { return }
        sinceSend = 0
        sent = head
        sentValid = head != nil
        engine.setTrackedHead(head)
    }
}
