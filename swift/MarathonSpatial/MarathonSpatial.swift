import CMarathonSpatial
import CMarathonSpatialOutput
import Foundation
import simd

/// A failed call into the library.
public struct SpatialError: Error, LocalizedError, Equatable {
    public let status: Int32
    public var errorDescription: String? { String(cString: ms_status_string(ms_status(rawValue: status))) }
}

@inline(__always)
private func check(_ status: ms_status) throws {
    if status != MS_OK { throw SpatialError(status: status.rawValue) }
}

@inline(__always)
private func vec(_ v: SIMD3<Float>) -> ms_vec3 { ms_vec3(x: v.x, y: v.y, z: v.z) }

/// A finished HRTF set. It is immutable and may serve several engines.
public final class HRTF: @unchecked Sendable {
    let handle: OpaquePointer

    private init(handle: OpaquePointer) { self.handle = handle }
    deinit { ms_hrtf_release(handle) }

    /// A cooked .mshrtf pack. Pass the engine's rate to resample while loading.
    public convenience init(contentsOf url: URL, sampleRate: Float = 0) throws {
        var result: OpaquePointer?
        try check(url.withUnsafeFileSystemRepresentation { ms_hrtf_load($0, sampleRate, &result) })
        self.init(handle: result!)
    }

    /// A user's AES69 SOFA file of type SimpleFreeFieldHRIR. Takes a moment; call it off the main thread.
    public convenience init(sofa url: URL, sampleRate: Float = 0, diffuseFieldEqualised: Bool = true) throws {
        var options = ms_hrtf_import_options()
        ms_hrtf_import_options_default(&options)
        options.sample_rate = sampleRate
        options.diffuse_field_eq = diffuseFieldEqualised ? 1 : 0
        var result: OpaquePointer?
        try check(url.withUnsafeFileSystemRepresentation { ms_hrtf_import_sofa($0, &options, &result) })
        self.init(handle: result!)
    }

    /// The analytic rigid sphere. Needs no data.
    public static func sphere(headRadius: Float = 0.0875, sampleRate: Float = 48000) throws -> HRTF {
        var result: OpaquePointer?
        try check(ms_hrtf_create_sphere(headRadius, sampleRate, &result))
        return HRTF(handle: result!)
    }

    public struct Info: Sendable, Equatable {
        public let name: String, licence: String
        public let sampleRate: Float, headRadius: Float, measuredDistance: Float
        public let directions: Int, taps: Int
    }

    public var info: Info {
        var i = ms_hrtf_info()
        ms_hrtf_get_info(handle, &i)
        let string: (UnsafeRawPointer) -> String = { String(cString: $0.assumingMemoryBound(to: CChar.self)) }
        return Info(name: withUnsafeBytes(of: &i.name) { string($0.baseAddress!) },
                    licence: withUnsafeBytes(of: &i.licence) { string($0.baseAddress!) },
                    sampleRate: i.sample_rate, headRadius: i.head_radius_m, measuredDistance: i.measured_distance_m,
                    directions: Int(i.directions), taps: Int(i.taps))
    }

    public func save(to url: URL) throws {
        try check(url.withUnsafeFileSystemRepresentation { ms_hrtf_save(handle, $0) })
    }
}

public typealias SourceParameters = ms_source_params
public typealias ReflectionTap = ms_reflection_tap

/// The spatial engine. Setters write a working copy and `commit()` publishes it to the render side as one snapshot.
/// All setters may be called from any thread. One thread at a time renders.
public final class SpatialEngine: @unchecked Sendable {
    public struct Configuration: Sendable {
        public var description: ms_engine_desc
        public init(sampleRate: Float = 48000, blockFrames: UInt32 = 128, maxSources: UInt32 = 256, maxDirectVoices: UInt32 = 64,
                    ambisonicOrder: UInt32 = 3, limiter: Bool = true) {
            var d = ms_engine_desc()
            ms_engine_desc_default(&d)
            d.sample_rate = sampleRate
            d.block_frames = blockFrames
            d.max_sources = maxSources
            d.max_hrtf_voices = maxDirectVoices
            d.ambisonic_order = ambisonicOrder
            d.limiter = limiter ? 1 : 0
            description = d
        }
    }

    public struct Source: Hashable, Sendable { public let id: ms_source }
    public struct Buffer: Hashable, Sendable { public let id: ms_buffer }
    public struct Bus: Hashable, Sendable {
        public let id: ms_bus
        public static let master = Bus(id: 1)
    }
    public struct Bed: Hashable, Sendable { public let id: ms_bed }

    let handle: OpaquePointer
    private var hrtf: HRTF?

    public init(_ configuration: Configuration = Configuration()) throws {
        var d = configuration.description
        var result: OpaquePointer?
        try check(ms_engine_create(&d, &result))
        handle = result!
    }

    deinit { ms_engine_destroy(handle) }

    public var sampleRate: Float { ms_engine_sample_rate(handle) }
    public var blockFrames: UInt32 { ms_engine_block_frames(handle) }

    public func commit() throws { try check(ms_engine_commit(handle)) }

    public func setHRTF(_ set: HRTF) throws {
        try check(ms_engine_set_hrtf(handle, set.handle))
        hrtf = set
    }

    public func setListener(position: SIMD3<Float>, velocity: SIMD3<Float> = .zero, forward: SIMD3<Float>, up: SIMD3<Float>) throws {
        try check(ms_listener_set_pose(handle, vec(position), vec(velocity), vec(forward), vec(up)))
    }

    /// Head orientation relative to the listener pose. Takes effect at the next block, with no commit.
    public func setHeadOrientation(_ q: simd_quatf) {
        ms_head_set_orientation(handle, ms_quat(x: q.imag.x, y: q.imag.y, z: q.imag.z, w: q.real))
    }

    public func makeBuffer(interleaved samples: [Float], channels: UInt32, sampleRate: Float) throws -> Buffer {
        var id: ms_buffer = 0
        try samples.withUnsafeBufferPointer {
            try check(ms_buffer_create(handle, $0.baseAddress, UInt64(samples.count) / UInt64(channels), channels, sampleRate, &id))
        }
        return Buffer(id: id)
    }

    /// Decodes any file Core Audio reads, on the calling thread.
    public func makeBuffer(contentsOf url: URL) throws -> Buffer {
        var id: ms_buffer = 0
        try check(url.withUnsafeFileSystemRepresentation { ms_buffer_load_file(handle, $0, &id) })
        return Buffer(id: id)
    }

    public func destroy(_ buffer: Buffer) throws { try check(ms_buffer_destroy(handle, buffer.id)) }

    public func makeSource() throws -> Source {
        var id: ms_source = 0
        try check(ms_source_create(handle, &id))
        return Source(id: id)
    }

    public func destroy(_ source: Source) throws { try check(ms_source_destroy(handle, source.id)) }

    public static var defaultSourceParameters: SourceParameters {
        var p = ms_source_params()
        ms_source_params_default(&p)
        return p
    }

    public func set(_ source: Source, parameters: SourceParameters) throws {
        var p = parameters
        try check(ms_source_set_params(handle, source.id, &p))
    }

    public func set(_ source: Source, position: SIMD3<Float>, velocity: SIMD3<Float> = .zero) throws {
        try check(ms_source_set_pose(handle, source.id, vec(position), vec(velocity)))
    }

    public func attach(_ buffer: Buffer, to source: Source, channel: UInt32 = 0, loop: Bool = false) throws {
        try check(ms_source_set_buffer(handle, source.id, buffer.id, channel, loop ? 1 : 0))
    }

    /// Pulls mono frames on the render side. The callback must not allocate, lock or block.
    public func attach(callback: @escaping ms_source_callback, user: UnsafeMutableRawPointer?, to source: Source) throws {
        try check(ms_source_set_callback(handle, source.id, callback, user))
    }

    public func play(_ source: Source, atFrame frame: UInt64 = 0) throws { try check(ms_source_play(handle, source.id, frame)) }
    public func stop(_ source: Source) throws { try check(ms_source_stop(handle, source.id)) }
    public func seek(_ source: Source, toFrame frame: Double) throws { try check(ms_source_seek(handle, source.id, frame)) }

    public func setReflections(_ taps: [ReflectionTap], for source: Source) throws {
        try taps.withUnsafeBufferPointer { try check(ms_source_set_reflections(handle, source.id, $0.baseAddress, UInt32(taps.count))) }
    }

    public func state(of source: Source) throws -> ms_source_state {
        var s = ms_source_state()
        try check(ms_source_get_state(handle, source.id, &s))
        return s
    }

    public func makeBed() throws -> Bed {
        var id: ms_bed = 0
        try check(ms_bed_create(handle, &id))
        return Bed(id: id)
    }

    public func destroy(_ bed: Bed) throws { try check(ms_bed_destroy(handle, bed.id)) }
    public func attach(_ buffer: Buffer, to bed: Bed, loop: Bool = true) throws { try check(ms_bed_set_buffer(handle, bed.id, buffer.id, loop ? 1 : 0)) }
    public func setGain(_ gain: Float, of bed: Bed) throws { try check(ms_bed_set_gain(handle, bed.id, gain)) }
    public func setOrientation(_ q: simd_quatf, of bed: Bed, headLocked: Bool = false) throws {
        try check(ms_bed_set_orientation(handle, bed.id, ms_quat(x: q.imag.x, y: q.imag.y, z: q.imag.z, w: q.real), headLocked ? 1 : 0))
    }
    public func play(_ bed: Bed, atFrame frame: UInt64 = 0) throws { try check(ms_bed_play(handle, bed.id, frame)) }
    public func stop(_ bed: Bed) throws { try check(ms_bed_stop(handle, bed.id)) }

    public func makeBus(parent: Bus = .master) throws -> Bus {
        var id: ms_bus = 0
        try check(ms_bus_create(handle, parent.id, &id))
        return Bus(id: id)
    }

    public func setGain(_ gain: Float, of bus: Bus) throws { try check(ms_bus_set_gain(handle, bus.id, gain)) }
    public func setFilter(highpass: Float, lowpass: Float, of bus: Bus) throws { try check(ms_bus_set_filter(handle, bus.id, highpass, lowpass)) }
    public func setProcessor(_ processor: ms_bus_processor?, user: UnsafeMutableRawPointer?, of bus: Bus) throws {
        try check(ms_bus_set_processor(handle, bus.id, processor, user))
    }

    /// Renders into memory. This is the same entry the live output calls, so an offline render is the live result.
    @discardableResult
    public func render(left: UnsafeMutablePointer<Float>, right: UnsafeMutablePointer<Float>, frames: UInt32) -> UInt64 {
        ms_engine_render(handle, left, right, frames)
    }

    /// Convenience offline render into interleaved stereo.
    public func renderInterleaved(frames: Int) -> [Float] {
        var out = [Float](repeating: 0, count: frames * 2)
        out.withUnsafeMutableBufferPointer { _ = ms_engine_render_interleaved(handle, $0.baseAddress, UInt32(frames)) }
        return out
    }

    public var stats: ms_engine_stats {
        var s = ms_engine_stats()
        ms_engine_get_stats(handle, &s)
        return s
    }

    public var operationCounts: ms_op_counts {
        var c = ms_op_counts()
        ms_engine_get_op_counts(handle, &c)
        return c
    }
}

/// Live output on a Core Audio device. Creating one binds to the device; only `start()` makes it run.
public final class SpatialOutput: @unchecked Sendable {
    private let handle: OpaquePointer
    private let engine: SpatialEngine

    public static func defaultDeviceSampleRate() throws -> Float {
        var rate: Float = 0
        try check(ms_output_default_sample_rate(&rate))
        return rate
    }

    /// `deviceUID` nil means the default output device. The engine must already run at the device's rate.
    public init(engine: SpatialEngine, deviceUID: String? = nil) throws {
        var result: OpaquePointer?
        if let deviceUID {
            try check(deviceUID.withCString { ms_output_create(engine.handle, $0, &result) })
        } else {
            try check(ms_output_create(engine.handle, nil, &result))
        }
        handle = result!
        self.engine = engine
    }

    deinit { ms_output_destroy(handle) }

    public func start() throws { try check(ms_output_start(handle)) }
    public func stop() throws { try check(ms_output_stop(handle)) }

    public var stats: ms_output_stats {
        var s = ms_output_stats()
        ms_output_get_stats(handle, &s)
        return s
    }
}
