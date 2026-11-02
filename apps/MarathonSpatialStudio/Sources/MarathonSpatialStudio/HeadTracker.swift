import CoreMotion
import StudioKit

/// Head tracking from headphones that report motion, such as AirPods Pro. The rotation goes straight to the library.
@MainActor
final class HeadTracker {
    private let manager = CMHeadphoneMotionManager()
    private var reference: CMQuaternion?
    private var latest: CMQuaternion?

    var isAvailable: Bool { manager.isDeviceMotionAvailable }

    func start(engine: SpatialEngine) {
        guard manager.isDeviceMotionAvailable, !manager.isDeviceMotionActive else { return }
        manager.startDeviceMotionUpdates(to: .main) { [weak self] motion, _ in
            guard let q = motion?.attitude.quaternion else { return }
            MainActor.assumeIsolated {
                guard let self else { return }
                self.latest = q
                let r = self.reference ?? q
                self.reference = r
                engine.setHeadOrientation(Self.engineOrientation(Self.relative(q, to: r)))
            }
        }
    }

    func stop(engine: SpatialEngine) {
        manager.stopDeviceMotionUpdates()
        reference = nil
        engine.setHeadOrientation(SIMD4<Float>(0, 0, 0, 1))
    }

    /// Make the current pose straight ahead.
    func recenter() { reference = latest }

    private static func relative(_ q: CMQuaternion, to r: CMQuaternion) -> CMQuaternion {
        // conjugate(r) * q
        let (ax, ay, az, aw) = (-r.x, -r.y, -r.z, r.w)
        return CMQuaternion(x: aw * q.x + ax * q.w + ay * q.z - az * q.y,
                            y: aw * q.y - ax * q.z + ay * q.w + az * q.x,
                            z: aw * q.z + ax * q.y - ay * q.x + az * q.w,
                            w: aw * q.w - ax * q.x - ay * q.y - az * q.z)
    }

    /// Core Motion's headphone frame is X right, Y forward, Z up. The library's is X right, Y up, -Z forward.
    private static func engineOrientation(_ q: CMQuaternion) -> SIMD4<Float> {
        SIMD4<Float>(Float(q.x), Float(q.z), Float(-q.y), Float(q.w))
    }
}
