#if MARATHON_SPATIAL_LINKED
import CStudioBridge
import Foundation
import StudioKit
import StudioRender

/// The real engine, on the library's C API: live output around the shared `EngineInstance`.
/// The real engine, on the library's C API.
final class LibraryEngine: SpatialEngine, @unchecked Sendable {
    let status: EngineStatus
    let preferredSampleRate: Double

    private let lock = NSLock()
    private var live: EngineInstance?
    private var output: OpaquePointer?
    private var configuration: EngineConfiguration?
    private var lastTick = ContinuousClock.now

    init() {
        var rate: Float = 48000
        let found = ms_output_default_sample_rate(&rate) == MS_OK && rate > 0
        preferredSampleRate = found ? Double(rate) : 48000
        status = EngineStatus(name: "Marathon Spatial", version: String(cString: ms_version_string()), canPlay: true, canRender: true,
                              canStrike: true, hasLateReverb: true, limitation: nil)
    }

    deinit {
        if let output {
            ms_output_stop(output)
            ms_output_destroy(output)
        }
    }

    func hrtfSets() async -> [HRTFSetInfo] {
        var out = [HRTFSetInfo(id: HRTFCatalog.sphereID, name: "Rigid sphere", detail: "Analytic head model, needs no data", licence: "Not applicable")]
        for url in HRTFCatalog.packs() {
            var h: OpaquePointer?
            guard ms_hrtf_load(url.path, 0, &h) == MS_OK, let h else { continue }
            out.append(info(id: url.lastPathComponent, h))
            ms_hrtf_release(h)
        }
        return out
    }

    private func info(id: String, _ h: OpaquePointer) -> HRTFSetInfo {
        var i = ms_hrtf_info()
        ms_hrtf_get_info(h, &i)
        func text<T>(_ tuple: T) -> String {
            withUnsafeBytes(of: tuple) { String(cString: $0.bindMemory(to: CChar.self).baseAddress!) }
        }
        return HRTFSetInfo(id: id, name: text(i.name), detail: "\(i.directions) directions, measured at \(String(format: "%.1f", i.measured_distance_m)) m", licence: text(i.licence))
    }

    func importHRTF(_ url: URL) async throws -> HRTFSetInfo {
        let scoped = url.startAccessingSecurityScopedResource()
        defer { if scoped { url.stopAccessingSecurityScopedResource() } }
        var o = ms_hrtf_import_options()
        ms_hrtf_import_options_default(&o)
        o.diffuse_field_eq = 1
        var h: OpaquePointer?
        try check(ms_hrtf_import_sofa(url.path, &o, &h), "Importing \(url.lastPathComponent)")
        defer { ms_hrtf_release(h) }
        try FileManager.default.createDirectory(at: HRTFCatalog.userDirectory, withIntermediateDirectories: true)
        let dest = HRTFCatalog.userDirectory.appendingPathComponent(url.deletingPathExtension().lastPathComponent + ".mshrtf")
        try check(ms_hrtf_save(h, dest.path), "Saving the HRTF set")
        return info(id: dest.lastPathComponent, h!)
    }

    func configure(_ c: EngineConfiguration) async throws {
        try lock.withLock {
            configuration = c
            if let instance = live, instance.limiter == c.render.limiter, instance.rate == c.sampleRate {
                try instance.apply(c)
            } else {
                closeOutput()
                let metallib = Bundle.main.url(forResource: "acoustic_trace", withExtension: "metallib")
                let instance = try EngineInstance(rate: c.sampleRate, limiter: c.render.limiter, metallib: metallib)
                try instance.apply(c)
                live = instance
            }
        }
    }

    private func closeOutput() {
        if let output {
            ms_output_stop(output)
            ms_output_destroy(output)
        }
        output = nil
    }

    func play(from seconds: Double) async throws {
        try lock.withLock {
            guard let instance = live else { throw UnavailableEngine.Missing(reason: "The scene is still loading.") }
            var o: OpaquePointer?
            try check(ms_output_create(instance.engine, nil, &o), "Opening the output device")
            output = o
            try instance.start(from: seconds, base: instance.framesRendered())
            try check(ms_output_start(o), "Starting output")
        }
    }

    func stop() async {
        lock.withLock {
            closeOutput()
            try? live?.halt()
        }
    }

    /// Skips rather than waits while a configure holds the lock; the next tick carries the newer pose anyway.
    func update(_ frame: EngineFrame) {
        guard lock.try() else { return }
        defer { lock.unlock() }
        guard let live else { return }
        try? live.publish(frame)
        let now = ContinuousClock.now
        let elapsed = lastTick.duration(to: now)
        if elapsed >= .milliseconds(16) {
            lastTick = now
            let seconds = Double(elapsed.components.seconds) + Double(elapsed.components.attoseconds) * 1e-18
            try? live.tickAcoustics(frame, seconds: min(seconds, 0.1))
        }
        try? live.commit()
    }

    func measuredRoom() async -> MeasuredRoom? {
        lock.withLock { live?.measuredRoom() }
    }

    func setHeadOrientation(_ q: SIMD4<Float>) {
        guard lock.try() else { return }
        defer { lock.unlock() }
        if let live { ms_head_set_orientation(live.engine, ms_quat(x: q.x, y: q.y, z: q.z, w: q.w)) }
    }

    func strike(_ strike: Strike) {
        lock.withLock { live?.strike(strike, delaySeconds: 0) }
    }

    func scrape(body: UUID, u: Double, v: Double, speedMetresPerSecond: Double, normalForceN: Double) {
        lock.withLock { live?.scrape(body: body, u: u, v: v, speed: speedMetresPerSecond, force: normalForceN) }
    }

    func endScrape(body: UUID) {
        lock.withLock { live?.endScrape(body: body) }
    }

    func meters() -> MeterFrame {
        guard lock.try() else { return MeterFrame() }
        defer { lock.unlock() }
        return live?.meters() ?? MeterFrame()
    }

    func modes(of body: UUID) async -> [BodyMode] {
        lock.withLock { live?.modes(of: body) ?? [] }
    }

    /// Its own engine instance, so an export neither opens a device nor disturbs what is playing.
    func render(_ configuration: EngineConfiguration, _ request: RenderRequest) -> AsyncThrowingStream<RenderBlock, Error> {
        AsyncThrowingStream { continuation in
            let task = Task.detached(priority: .userInitiated) {
                do {
                    try OfflineRender.run(configuration, request) { continuation.yield($0) }
                    continuation.finish()
                } catch {
                    continuation.finish(throwing: error)
                }
            }
            continuation.onTermination = { _ in task.cancel() }
        }
    }
}
#endif
