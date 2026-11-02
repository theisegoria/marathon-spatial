import CoreMotion
import Foundation
import simd

/// Optional head tracking from headphones with motion sensors, through Core Motion. It feeds the engine directly, so an
/// orientation reaches the next audio block without a commit. The host app needs NSMotionUsageDescription.
public final class HeadTracker: @unchecked Sendable {
    private let manager = CMHeadphoneMotionManager()
    private let queue = OperationQueue()
    private let engine: SpatialEngine
    private var reference: simd_quatf?

    public init(engine: SpatialEngine) {
        self.engine = engine
        queue.maxConcurrentOperationCount = 1
        queue.qualityOfService = .userInteractive
    }

    public var isAvailable: Bool { manager.isDeviceMotionAvailable }

    /// Starts tracking. The pose at the first update, or at the next `recentre()`, becomes straight ahead.
    public func start() {
        manager.startDeviceMotionUpdates(to: queue) { [weak self] motion, _ in
            guard let self, let q = motion?.attitude.quaternion else { return }
            // Core Motion's headphone frame: x right, y forward, z up. The engine's head frame: x right, y up, z back.
            let attitude = simd_quatf(ix: Float(q.x), iy: Float(q.z), iz: Float(-q.y), r: Float(q.w))
            let origin = self.reference ?? attitude
            self.reference = origin
            self.engine.setHeadOrientation(origin.inverse * attitude)
        }
    }

    public func recentre() { queue.addOperation { [weak self] in self?.reference = nil } }

    public func stop() {
        manager.stopDeviceMotionUpdates()
        engine.setHeadOrientation(simd_quatf(ix: 0, iy: 0, iz: 0, r: 1))
    }
}
