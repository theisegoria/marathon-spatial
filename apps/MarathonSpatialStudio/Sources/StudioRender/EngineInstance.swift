import CStudioBridge
import Foundation
import StudioKit

public struct LibraryError: Error, LocalizedError {
    public var status: ms_status
    public var context: String
    public var errorDescription: String? { "\(context): \(String(cString: ms_status_string(status)))" }
}

public func check(_ status: ms_status, _ context: String) throws {
    if status != MS_OK { throw LibraryError(status: status, context: context) }
}

private func vec(_ v: Vec3) -> ms_vec3 { ms_vec3(x: Float(v.x), y: Float(v.y), z: Float(v.z)) }
private func linear(_ dB: Double) -> Float { Float(pow(10, dB / 20)) }
private let masterBus: ms_bus = 1

/// One ms_engine with the handles of everything loaded into it. Not thread-safe: LibraryEngine guards the live one, and an
/// export owns its own.
/// With no output library behind it this type cannot open a device, which is why the render tests link it alone.
public final class EngineInstance {
    enum Content: Hashable {
        case clip(UUID)
        case signal(TestSignal)
    }
    struct SourceSlot {
        var handle: ms_source
        var setup: SoundSource
        var content: SourceContent
        var params: ms_source_params
        var buffer: ms_buffer
        var frames: UInt64
        /// Has been sent a play since it last changed content.
        var started = false
    }
    struct BodySlot {
        var handle: ms_source
        /// One physics instance per body: its single voice is this body's source, so the sound stays where the body is.
        var physics: OpaquePointer
        var model: ms_modal_model
        var surface: ms_surface
        var contact: ms_physics_material
        var body: SoundingBody
        var params: ms_source_params
        var material: PhysicalMaterial
        var modes: [BodyMode]
    }
    struct BedSlot {
        var handle: ms_bed
        var bed: AmbisonicBed
    }

    public let engine: OpaquePointer
    public let rate: Double
    public let limiter: Bool
    private var sources: [UUID: SourceSlot] = [:]
    private var bodies: [UUID: BodySlot] = [:]
    private var beds: [UUID: BedSlot] = [:]
    private var buffers: [Content: (buffer: ms_buffer, frames: UInt64)] = [:]
    private var bedBuffers: [UUID: ms_buffer] = [:]
    /// Physics instances stay alive until the engine is destroyed: the render side may still hold their callback.
    private var physics: [OpaquePointer] = []
    private var hrtfID: String??
    private var stereo = false
    /// The canceller as the document sets it, before any tracked head, or nil when the output is not 3D speakers.
    private var crosstalk: ms_crosstalk_params?
    private var trackedHead = ms_head_position()
    /// The room's simulator, its late reverb and the bus that runs it. One room, one network.
    private var acoustics: OpaquePointer?
    private var reverb: OpaquePointer?
    private var reverbBus: ms_bus = 0
    private var acousticKey: AcousticKey?
    private var acousticIDs: [UUID: UInt32] = [:]
    private var earlyReflections = true
    private var lateReverb = true
    private let metallib: URL?
    /// Replaced simulators the render side may still be reading, with the frame after which it cannot be.
    private var retired: [(acoustics: OpaquePointer, reverb: OpaquePointer?, safeAt: UInt64)] = []
    /// Timeline origin while playing: a source added or unmuted mid-play joins at the right place.
    private var running: (seconds: Double, base: UInt64)?

    /// `metallib` runs reflection rays on the GPU. An export passes nil: GPU rays finish on their own schedule, and it must repeat.
    public init(rate: Double, limiter: Bool, metallib: URL? = nil) throws {
        self.metallib = metallib
        var d = ms_engine_desc()
        ms_engine_desc_default(&d)
        d.sample_rate = Float(rate)
        d.block_frames = 256
        d.max_sources = 128
        d.max_hrtf_voices = 32
        d.max_beds = 8
        d.max_buffers = 256
        d.ambisonic_order = 3
        d.limiter = limiter ? 1 : 0
        var e: OpaquePointer?
        try check(ms_engine_create(&d, &e), "Creating the engine")
        engine = e!
        self.rate = rate
        self.limiter = limiter
        try check(ms_bus_create(engine, masterBus, &reverbBus), "Creating the reverb bus")
    }

    deinit {
        ms_engine_destroy(engine)
        physics.forEach { ms_physics_destroy($0) }
        retired.forEach { destroyAcoustics($0.acoustics, $0.reverb) }
        if let acoustics { destroyAcoustics(acoustics, reverb) }
    }

    public func apply(_ c: EngineConfiguration) throws {
        stereo = c.render.output == .stereo
        try applyCrosstalk(c.render)
        earlyReflections = c.render.earlyReflections
        lateReverb = c.render.lateReverb
        try selectHRTF(c.render.hrtfSet)
        try applyAcoustics(c)
        try applySources(c)
        try applyBodies(c)
        try applyBeds(c)
        syncAcousticSources(c)
        try check(ms_bus_set_gain(engine, reverbBus, linear(c.render.reverbTrimDB)), "Setting the reverb trim")
        try check(ms_bus_set_gain(engine, masterBus, linear(c.render.masterGainDB)), "Setting the master gain")
        if let running {
            let now = running.seconds + Double(framesRendered() - running.base) / rate
            for id in sources.keys where sources[id]?.started == false && sources[id]?.buffer != 0 {
                try startSource(id, from: now, base: 0)
            }
        }
        try commit()
    }

    public func commit() throws { try check(ms_engine_commit(engine), "Publishing the scene") }

    private func applyCrosstalk(_ render: RenderSettings) throws {
        guard render.output == .speakers3D else {
            if crosstalk != nil { try check(ms_engine_set_crosstalk(engine, nil), "Turning crosstalk cancellation off") }
            crosstalk = nil
            return
        }
        var p = ms_crosstalk_params()
        ms_crosstalk_params_default(&p)
        p.speaker_span_degrees = Float(render.crosstalkSettings.spanDegrees)
        p.listener_distance_m = Float(render.crosstalkSettings.distanceMetres)
        if let crosstalk, crosstalk.speaker_span_degrees == p.speaker_span_degrees, crosstalk.listener_distance_m == p.listener_distance_m { return }
        crosstalk = p
        try steer()
    }

    private func steer() throws {
        guard var p = crosstalk else { return }
        var head = trackedHead
        ms_crosstalk_params_follow_head(&p, &head)
        try check(ms_engine_set_crosstalk(engine, &p), "Setting crosstalk cancellation")
    }

    /// Re-steers the canceller to where the listener sits, or back to the centre with nil. Designs on this thread, a few milliseconds.
    public func setTrackedHead(_ head: TrackedHeadPosition?) throws {
        trackedHead = head.map { ms_head_position(x_m: Float($0.xMetres), z_m: Float($0.zMetres), valid: 1) } ?? ms_head_position()
        guard crosstalk != nil else { return }
        try steer()
        try commit()
    }

    /// The canceller's added delay while it runs, in frames.
    public var crosstalkLatencyFrames: UInt32 { ms_engine_crosstalk_latency_frames(engine) }

    private func destroyAcoustics(_ a: OpaquePointer, _ r: OpaquePointer?) {
        if let r { studio_reverb_destroy(r) }
        ms_acoustics_destroy(a)
    }

    /// What changes the simulator's scene. Moving sources and the listener do not.
    private struct AcousticKey: Equatable {
        var mesh: EngineMesh
        var materials: [AcousticMaterial]
        var room: Room
    }

    /// A new scene is built whole and swapped in, so the render thread is never touched: it reads the old one until the bus changes over.
    private func applyAcoustics(_ c: EngineConfiguration) throws {
        let frames = framesRendered()
        retired.removeAll { old in
            guard frames >= old.safeAt else { return false }
            destroyAcoustics(old.acoustics, old.reverb)
            return true
        }
        let key = AcousticKey(mesh: c.mesh, materials: c.materials, room: c.room)
        if key == acousticKey { return }
        guard !c.mesh.indices.isEmpty else { return }

        let bands = AcousticMaterial.bandsHz.count
        let absorption = c.materials.flatMap { $0.absorption.map(Float.init) }
        let scattering = c.materials.map { Float($0.scattering) }
        // The document holds transmission in dB of energy.
        let transmission = c.materials.map { Float(min(1, pow(10, $0.transmissionDB / 10))) }
        let surfaces = c.room.surfaces()
        var absorptionArea = [Float](repeating: 0, count: bands)
        for s in surfaces {
            let alpha = c.materials.first { $0.id == s.material }?.absorption ?? [Double](repeating: 0.1, count: bands)
            for b in 0..<bands { absorptionArea[b] += Float(s.area * alpha[b]) }
        }
        let xs = c.room.outline.map(\.x), zs = c.room.outline.map(\.z)
        let centre = c.room.centroid
        var desc = studio_room_desc()
        desc.bounds_min = ms_vec3(x: Float(xs.min() ?? 0), y: 0, z: Float(zs.min() ?? 0))
        desc.bounds_max = ms_vec3(x: Float(xs.max() ?? 0), y: Float(c.room.height), z: Float(zs.max() ?? 0))
        desc.probe = ms_vec3(x: Float(centre.x), y: Float(c.room.height / 2), z: Float(centre.z))
        desc.volume = Float(c.room.volume)
        desc.surface_area = Float(surfaces.reduce(0) { $0 + $1.area })
        desc.material_count = UInt32(c.materials.count)
        desc.vertex_count = UInt32(c.mesh.positions.count / 3)
        desc.index_count = UInt32(c.mesh.indices.count)

        let path = metallib?.path
        let created: OpaquePointer? = c.mesh.positions.withUnsafeBufferPointer { positions in
            c.mesh.indices.withUnsafeBufferPointer { indices in
                c.mesh.triangleMaterials.withUnsafeBufferPointer { triangles in
                    absorption.withUnsafeBufferPointer { a in
                        scattering.withUnsafeBufferPointer { s in
                            transmission.withUnsafeBufferPointer { t in
                                absorptionArea.withUnsafeBufferPointer { area in
                                    desc.positions = positions.baseAddress
                                    desc.indices = indices.baseAddress
                                    desc.triangle_materials = triangles.baseAddress
                                    desc.absorption = a.baseAddress
                                    desc.scattering = s.baseAddress
                                    desc.transmission = t.baseAddress
                                    desc.absorption_area = area.baseAddress
                                    return studio_acoustics_create(&desc, path, Float(rate))
                                }
                            }
                        }
                    }
                }
            }
        }
        guard let created else { throw LibraryError(status: MS_ERR_INVALID_ARGUMENT, context: "Building the room's acoustics") }
        guard let network = studio_reverb_create(created, 4096) else {
            ms_acoustics_destroy(created)
            throw LibraryError(status: MS_ERR_INVALID_ARGUMENT, context: "Building the late reverb")
        }
        try check(ms_bus_set_processor(engine, reverbBus, studio_reverb_bus, UnsafeMutableRawPointer(network)), "Connecting the late reverb")
        try commit()
        if let acoustics { retired.append((acoustics, reverb, framesRendered() + 1024)) }
        acoustics = created
        reverb = network
        acousticKey = key
        acousticIDs = [:]
    }

    /// Sources and bodies alike reflect, up to the simulator's limit. Past it a voice simply has no room response.
    private func syncAcousticSources(_ c: EngineConfiguration) {
        guard let acoustics else { return }
        let places = c.sources.map { ($0.source.id, $0.source.position) } + c.bodies.map { ($0.body.id, $0.body.position) }
        let wanted = Set(places.map(\.0))
        for (id, slot) in acousticIDs where !wanted.contains(id) {
            ms_acoustics_source_remove(acoustics, slot)
            acousticIDs[id] = nil
        }
        for (id, position) in places where acousticIDs[id] == nil && acousticIDs.count < Int(MS_ACOUSTIC_MAX_SOURCES) {
            acousticIDs[id] = ms_acoustics_source_add(acoustics, vec(position), 0.1)
        }
    }

    /// One step of the simulator, then its reflections, occlusion and reverb send for every voice. Control thread.
    public func tickAcoustics(_ f: EngineFrame, seconds: Double) throws {
        guard let acoustics else { return }
        ms_acoustics_set_listener(acoustics, vec(f.listenerPosition))
        // The convolution reverb's stereo decoding turns with the listener.
        ms_acoustics_set_listener_orientation(acoustics, vec(f.listenerForward), vec(f.listenerUp))
        for pose in f.sources + f.bodies {
            if let id = acousticIDs[pose.id] { ms_acoustics_source_set(acoustics, id, vec(pose.position), 0.1) }
        }
        ms_acoustics_tick(acoustics, Float(seconds))
        for pose in f.sources {
            guard var slot = sources[pose.id] else { continue }
            try respond(pose, f, handle: slot.handle, params: &slot.params, reflects: slot.setup.reflections)
            sources[pose.id] = slot
        }
        for pose in f.bodies {
            guard var slot = bodies[pose.id] else { continue }
            try respond(pose, f, handle: slot.handle, params: &slot.params, reflects: true)
            bodies[pose.id] = slot
        }
    }

    private func respond(_ pose: EngineSourcePose, _ f: EngineFrame, handle: ms_source, params p: inout ms_source_params, reflects: Bool) throws {
        var view = ms_acoustic_source_view()
        var taps = [ms_reflection_tap](repeating: ms_reflection_tap(), count: Int(MS_MAX_REFLECTION_TAPS))
        var count: UInt32 = 0
        if let acoustics, let id = acousticIDs[pose.id] {
            ms_acoustics_source_view_get(acoustics, id, &view)
            // The pan path renders no reflections.
            if reflects, earlyReflections, !stereo { count = ms_acoustics_source_taps(acoustics, id, &taps) }
        }
        p.position = vec(pose.position)
        p.velocity = vec(pose.velocity)
        p.eq_low_db = view.active != 0 ? view.eq_low_db : 0
        p.eq_mid_db = view.active != 0 ? view.eq_mid_db : 0
        p.eq_high_db = view.active != 0 ? view.eq_high_db : 0
        // A send carries the voice after its distance gain, and the room wants the source at 1 m.
        let d = pose.position - f.listenerPosition
        p.send_bus.0 = lateReverb ? reverbBus : 0
        p.send_gain.0 = lateReverb ? Float(max(((d * d).sum() as Double).squareRoot(), 1)) : 0
        try check(ms_source_set_params(engine, handle, &p), "Setting a source's acoustics")
        try check(ms_source_set_reflections(engine, handle, taps, count), "Publishing reflections")
    }

    /// The simulator's own view of the room, for the inspector.
    public func measuredRoom() -> MeasuredRoom? {
        guard let acoustics else { return nil }
        var view = ms_acoustic_room_view()
        ms_acoustics_room_view_get(acoustics, 0, &view)
        guard view.active != 0 else { return nil }
        let t60 = withUnsafeBytes(of: view.t60) { Array($0.bindMemory(to: Float.self)) }
        var energy = ms_acoustic_energy_view()
        ms_acoustics_energy_view_get(acoustics, &energy)
        let traced = energy.seconds > 0 ? withUnsafeBytes(of: energy.t60) { Array($0.bindMemory(to: Float.self)) }.map(Double.init) : []
        return MeasuredRoom(volume: Double(view.volume), surfaceArea: Double(view.surface_area), t60Seconds: t60.map(Double.init),
                            tracedT60Seconds: traced)
    }

    private func selectHRTF(_ id: String?) throws {
        if case let .some(current) = hrtfID, current == id { return }
        var h: OpaquePointer?
        if let id, let url = HRTFCatalog.url(for: id) {
            try check(ms_hrtf_load(url.path, Float(rate), &h), "Loading the HRTF set")
        } else if let url = id == nil ? HRTFCatalog.defaultPack() : nil {
            try check(ms_hrtf_load(url.path, Float(rate), &h), "Loading the HRTF set")
        } else {
            try check(ms_hrtf_create_sphere(0.0875, Float(rate), &h), "Building the sphere HRTF")
        }
        // The engine keeps its own reference.
        defer { ms_hrtf_release(h) }
        try check(ms_engine_set_hrtf(engine, h), "Selecting the HRTF set")
        hrtfID = .some(id)
    }

    private func params(for s: SoundSource) -> ms_source_params {
        var p = ms_source_params()
        ms_source_params_default(&p)
        let a = s.directivity.yawDegrees * .pi / 180, b = s.directivity.pitchDegrees * .pi / 180
        p.forward = ms_vec3(x: Float(sin(a) * cos(b)), y: Float(sin(b)), z: Float(-cos(a) * cos(b)))
        p.gain = linear(s.gainDB)
        p.mode = (stereo ? MS_SPATIAL_PAN : MS_SPATIAL_HRTF).rawValue
        p.distance_model = MS_DISTANCE_INVERSE.rawValue
        p.ref_distance = 1
        p.max_distance = 1000
        p.rolloff = 1
        p.doppler = s.doppler ? 1 : 0
        p.air_absorption = s.airAbsorption ? 1 : 0
        p.near_field = s.nearField ? 1 : 0
        p.propagation_delay = 1
        p.directivity_alpha = Float(s.directivity.dipoleWeight)
        p.directivity_sharpness = Float(s.directivity.power)
        return p
    }

    private func buffer(for setup: EngineSourceSetup) throws -> (buffer: ms_buffer, frames: UInt64)? {
        let key: Content
        switch setup.source.content {
        case .none: return nil
        case .clip(let id): key = .clip(id)
        case .signal(let s): key = .signal(s)
        }
        if let cached = buffers[key] { return cached }
        guard let clip = setup.clip, clip.frames > 0 else { return nil }
        let mono = clip.monoMix
        var b: ms_buffer = 0
        try check(ms_buffer_create(engine, mono, UInt64(mono.count), 1, Float(clip.sampleRate), &b), "Loading audio")
        var frames: UInt64 = 0, channels: UInt32 = 0
        try check(ms_buffer_get_info(engine, b, &frames, &channels), "Reading audio")
        buffers[key] = (b, frames)
        return (b, frames)
    }

    private func applySources(_ c: EngineConfiguration) throws {
        let wanted = Set(c.sources.map(\.source.id))
        for (id, slot) in sources where !wanted.contains(id) {
            try check(ms_source_destroy(engine, slot.handle), "Removing a source")
            sources[id] = nil
        }
        for setup in c.sources {
            let s = setup.source
            var slot: SourceSlot
            if let existing = sources[s.id] {
                slot = existing
            } else {
                var h: ms_source = 0
                try check(ms_source_create(engine, &h), "Creating a source")
                slot = SourceSlot(handle: h, setup: s, content: .none, params: ms_source_params(), buffer: 0, frames: 0)
            }
            if slot.content != s.content || slot.setup.loops != s.loops || slot.buffer == 0 {
                if let loaded = try buffer(for: setup) {
                    try check(ms_source_set_buffer(engine, slot.handle, loaded.buffer, 0, s.loops ? 1 : 0), "Assigning audio")
                    slot.buffer = loaded.buffer
                    slot.frames = loaded.frames
                } else {
                    slot.buffer = 0
                }
                slot.content = s.content
                slot.started = false
            }
            slot.setup = s
            var p = params(for: s)
            if sources[s.id] != nil {
                // An edit during playback keeps where the voice is and what the room is doing to it until the next tick.
                p.position = slot.params.position
                p.velocity = slot.params.velocity
                p.eq_low_db = slot.params.eq_low_db; p.eq_mid_db = slot.params.eq_mid_db; p.eq_high_db = slot.params.eq_high_db
                p.send_bus = slot.params.send_bus; p.send_gain = slot.params.send_gain
            } else {
                p.position = vec(s.position)
            }
            slot.params = p
            try check(ms_source_set_params(engine, slot.handle, &slot.params), "Setting a source")
            sources[s.id] = slot
        }
    }

    private func applyBodies(_ c: EngineConfiguration) throws {
        let wanted = Set(c.bodies.map(\.body.id))
        for (id, slot) in bodies where !wanted.contains(id) {
            try check(ms_source_destroy(engine, slot.handle), "Removing a body")
            bodies[id] = nil
        }
        for setup in c.bodies {
            let b = setup.body
            if let existing = bodies[b.id], existing.body.size == b.size, existing.body.wallThicknessM == b.wallThicknessM, existing.material == setup.material {
                var p = params(forBody: b)
                p.eq_low_db = existing.params.eq_low_db; p.eq_mid_db = existing.params.eq_mid_db; p.eq_high_db = existing.params.eq_high_db
                p.send_bus = existing.params.send_bus; p.send_gain = existing.params.send_gain
                try check(ms_source_set_params(engine, existing.handle, &p), "Setting a body")
                bodies[b.id]?.body = b
                bodies[b.id]?.params = p
                continue
            }
            if let old = bodies[b.id] { try check(ms_source_destroy(engine, old.handle), "Replacing a body") }
            bodies[b.id] = nil

            var material = ms_physics_material()
            if ms_physics_material_preset(setup.material.id, &material) != MS_OK {
                try check(ms_physics_material_preset("steel", &material), "Finding a material preset")
            }
            material.young_modulus_pa = Float(setup.material.youngsModulusPa)
            material.density_kg_m3 = Float(setup.material.densityKgM3)
            material.poisson_ratio = Float(setup.material.poissonRatio)
            material.loss_factor = Float(setup.material.lossFactor)

            var d = ms_physics_desc()
            ms_physics_desc_default(&d)
            d.sample_rate = Float(rate)
            d.voices = 1
            d.contacts = 1
            var created: OpaquePointer?
            try check(ms_physics_create(&d, &created), "Creating physics for \(b.name)")
            physics.append(created!)
            let size = ms_vec3(x: Float(b.size.x), y: Float(b.size.y), z: Float(b.size.z))
            var model: ms_modal_model = 0, surface: ms_surface = 0
            try check(ms_physics_model_box(created, &material, size, Float(b.wallThicknessM), &model), "Analysing \(b.name)")
            try check(ms_physics_surface_create(created, material.corner_wavelength_m, 1, &surface), "Texturing \(b.name)")

            var h: ms_source = 0
            try check(ms_source_create(engine, &h), "Creating a body")
            try check(ms_source_set_callback(engine, h, ms_physics_source_callback, ms_physics_voice(created, 0)), "Connecting a body")
            var p = params(forBody: b)
            try check(ms_source_set_params(engine, h, &p), "Setting a body")
            try check(ms_source_play(engine, h, 0), "Starting a body")

            var hz = [Double](repeating: 0, count: 64), damping = hz
            let count = Int(studio_box_modes(&material, size, Float(b.wallThicknessM), &hz, &damping, 64))
            let modes = (0..<count).map { BodyMode(frequencyHz: hz[$0], dampingRatio: damping[$0]) }
            bodies[b.id] = BodySlot(handle: h, physics: created!, model: model, surface: surface, contact: material, body: b, params: p,
                                    material: setup.material, modes: modes)
        }
    }

    private func params(forBody b: SoundingBody) -> ms_source_params {
        var s = SoundSource(name: b.name, position: b.position)
        s.gainDB = b.gainDB
        s.doppler = false
        var p = params(for: s)
        p.position = vec(b.position)
        return p
    }

    private func applyBeds(_ c: EngineConfiguration) throws {
        let wanted = Set(c.beds.map(\.bed.id))
        for (id, slot) in beds where !wanted.contains(id) {
            try check(ms_bed_destroy(engine, slot.handle), "Removing a bed")
            beds[id] = nil
        }
        for setup in c.beds {
            let bed = setup.bed
            var slot: BedSlot
            if let existing = beds[bed.id] {
                slot = existing
            } else {
                var h: ms_bed = 0
                try check(ms_bed_create(engine, &h), "Creating a bed")
                slot = BedSlot(handle: h, bed: bed)
            }
            var b = bedBuffers[bed.clip] ?? 0
            if b == 0 {
                let data = setup.clip.interleaved
                try check(ms_buffer_create(engine, data, UInt64(setup.clip.frames), UInt32(setup.clip.channels.count), Float(setup.clip.sampleRate), &b), "Loading a bed")
                bedBuffers[bed.clip] = b
            }
            try check(ms_bed_set_buffer(engine, slot.handle, b, bed.loops ? 1 : 0), "Assigning a bed")
            try check(ms_bed_set_gain(engine, slot.handle, linear(bed.gainDB)), "Setting a bed")
            let half = bed.yawDegrees * .pi / 360
            let q = ms_quat(x: 0, y: Float(-sin(half)), z: 0, w: Float(cos(half)))
            try check(ms_bed_set_orientation(engine, slot.handle, q, bed.headLocked ? 1 : 0), "Turning a bed")
            slot.bed = bed
            beds[bed.id] = slot
        }
    }

    /// Everything starts relative to absolute engine frame `base`, from timeline position `seconds`.
    public func start(from seconds: Double, base: UInt64) throws {
        running = (seconds, base)
        for id in sources.keys where sources[id]?.buffer != 0 { try startSource(id, from: seconds, base: base) }
        for slot in beds.values { try check(ms_bed_play(engine, slot.handle, base), "Starting a bed") }
        try commit()
    }

    /// A `base` of 0 means the next block.
    private func startSource(_ id: UUID, from seconds: Double, base: UInt64) throws {
        guard let slot = sources[id] else { return }
        sources[id]?.started = true
        let offset = seconds - slot.setup.startSeconds
        if offset >= 0 {
            var at = offset * rate
            if slot.setup.loops, slot.frames > 0 { at = fmod(at, Double(slot.frames)) }
            if !slot.setup.loops, at >= Double(slot.frames) { return }
            try check(ms_source_seek(engine, slot.handle, at), "Seeking")
            try check(ms_source_play(engine, slot.handle, base), "Starting a source")
        } else {
            try check(ms_source_seek(engine, slot.handle, 0), "Seeking")
            try check(ms_source_play(engine, slot.handle, max(base, framesRendered()) + UInt64(-offset * rate)), "Starting a source")
        }
    }

    public func halt() throws {
        running = nil
        for id in sources.keys { sources[id]?.started = false }
        for slot in sources.values { try check(ms_source_stop(engine, slot.handle), "Stopping a source") }
        for slot in beds.values { try check(ms_bed_stop(engine, slot.handle), "Stopping a bed") }
        try commit()
    }

    public func framesRendered() -> UInt64 {
        var s = ms_engine_stats()
        ms_engine_get_stats(engine, &s)
        return s.frames_rendered
    }

    public func publish(_ f: EngineFrame) throws {
        try check(ms_listener_set_pose(engine, vec(f.listenerPosition), vec(f.listenerVelocity), vec(f.listenerForward), vec(f.listenerUp)), "Moving the listener")
        // Right is forward x up, for a right-handed frame with -Z forward.
        let r = Vec3(f.listenerForward.y * f.listenerUp.z - f.listenerForward.z * f.listenerUp.y,
                     f.listenerForward.z * f.listenerUp.x - f.listenerForward.x * f.listenerUp.z,
                     f.listenerForward.x * f.listenerUp.y - f.listenerForward.y * f.listenerUp.x)
        for pose in f.sources {
            guard var slot = sources[pose.id] else { continue }
            if stereo {
                let d = pose.position - f.listenerPosition
                let l = max(((d * d).sum() as Double).squareRoot(), 1e-3)
                slot.params.position = vec(pose.position)
                slot.params.velocity = vec(pose.velocity)
                slot.params.pan = Float(((d * r).sum() as Double) / l)
                // The library's pan mode is head locked and has no distance law, so the stereo path applies it here.
                slot.params.gain = linear(slot.setup.gainDB) / Float(max(l, 1))
                try check(ms_source_set_params(engine, slot.handle, &slot.params), "Panning a source")
                sources[pose.id] = slot
            } else {
                try check(ms_source_set_pose(engine, slot.handle, vec(pose.position), vec(pose.velocity)), "Moving a source")
            }
        }
        for pose in f.bodies {
            if let slot = bodies[pose.id] {
                try check(ms_source_set_pose(engine, slot.handle, vec(pose.position), vec(pose.velocity)), "Moving a body")
            }
        }
    }

    /// `delaySeconds` counts from the next rendered frame, so an export lands a strike on its exact frame.
    public func strike(_ s: Strike, delaySeconds: Float) {
        guard let slot = bodies[s.body] else { return }
        var d = ms_physics_impact_desc()
        d.instance = 1
        d.model = slot.model
        d.local_position = vec(slot.body.localPosition(u: s.u, v: s.v))
        d.impulse_ns = Float(s.impulseNs)
        d.contact_seconds = Float(s.contactSeconds)
        d.weight = 1
        d.delay_seconds = delaySeconds
        _ = ms_physics_impact(slot.physics, &d)
    }

    public func scrape(body: UUID, u: Double, v: Double, speed: Double, force: Double) {
        guard let slot = bodies[body] else { return }
        var d = ms_physics_contact_desc()
        ms_physics_contact_desc_default(&d)
        d.contact = 1
        d.model = slot.model
        d.surface = slot.surface
        d.local_position = vec(slot.body.localPosition(u: u, v: v))
        d.kind = MS_CONTACT_SLIDE
        d.speed_m_s = Float(speed)
        d.normal_force_n = Float(force)
        d.dynamic_friction = slot.contact.dynamic_friction
        d.static_friction = slot.contact.static_friction
        d.fluctuation = slot.contact.slide_fluctuation
        d.stick_stiffness_n_m = slot.contact.stick_stiffness_n_m
        _ = ms_physics_contact_update(slot.physics, &d)
    }

    public func endScrape(body: UUID) {
        if let slot = bodies[body] { _ = ms_physics_contact_end(slot.physics, 1, 0.05) }
    }

    public func modes(of id: UUID) -> [BodyMode] { bodies[id]?.modes ?? [] }

    public func meters() -> MeterFrame {
        var s = ms_engine_stats()
        ms_engine_get_stats(engine, &s)
        var m = MeterFrame()
        m.peak = [s.peak.0, s.peak.1]
        m.rms = [s.rms.0, s.rms.1]
        m.overs = Int(clamping: s.overs)
        for (id, slot) in sources {
            var state = ms_source_state()
            if ms_source_get_state(engine, slot.handle, &state) == MS_OK { m.sources[id] = state.level }
        }
        return m
    }
}

/// Renders a whole request faster than real time into memory, on its own engine.
public enum OfflineRender {
    static let tickFrames = 2048.0
    static let settleTicks = 30

    public static func run(_ c: EngineConfiguration, _ r: RenderRequest, yield: (RenderBlock) -> Void) throws {
        let instance = try EngineInstance(rate: c.sampleRate, limiter: c.render.limiter)
        try instance.apply(c)
        try instance.publish(r.frame(r.startSeconds))
        // Early reflections build up over ticks. Settle them first, so an export has them from its first frame.
        for _ in 0..<settleTicks { try instance.tickAcoustics(r.frame(r.startSeconds), seconds: tickFrames / c.sampleRate) }
        try instance.start(from: r.startSeconds, base: 0)

        let total = Int64(r.durationSeconds * c.sampleRate)
        let block = max(1, r.blockFrames)
        var done: Int64 = 0
        var nextStrike = 0
        var sinceTaps = 0
        var samples = [Float](repeating: 0, count: block * 2)
        while done < total {
            try Task.checkCancellation()
            let frames = Int(min(Int64(block), total - done))
            let time = r.startSeconds + Double(done) / c.sampleRate
            let frame = r.frame(time)
            try instance.publish(frame)
            sinceTaps += frames
            if Double(sinceTaps) >= tickFrames {
                try instance.tickAcoustics(frame, seconds: Double(sinceTaps) / c.sampleRate)
                sinceTaps = 0
            }
            let end = done + Int64(frames)
            while nextStrike < r.strikes.count {
                let s = r.strikes[nextStrike]
                let at = Int64((s.timeSeconds - r.startSeconds) * c.sampleRate)
                guard at < end else { break }
                instance.strike(s, delaySeconds: Float(Double(max(0, at - done)) / c.sampleRate))
                nextStrike += 1
            }
            try instance.commit()
            ms_engine_render_interleaved(instance.engine, &samples, UInt32(frames))
            done = end
            yield(RenderBlock(samples: Array(samples[0..<(frames * 2)]), framesDone: done, framesTotal: total))
        }
    }
}

/// Where HRTF packs live: shipped in the app, or imported by the user.
public enum HRTFCatalog {
    public static let sphereID = "sphere"

    public static var userDirectory: URL {
        FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask)[0]
            .appendingPathComponent("Marathon Spatial Studio/HRTF")
    }

    public static func packs() -> [URL] {
        let shipped = Bundle.main.resourceURL?.appendingPathComponent("HRTF")
        return ([shipped, userDirectory].compactMap { $0 }).flatMap { dir in
            ((try? FileManager.default.contentsOfDirectory(at: dir, includingPropertiesForKeys: nil)) ?? [])
                .filter { $0.pathExtension == "mshrtf" }
                .sorted { $0.lastPathComponent < $1.lastPathComponent }
        }
    }

    public static func url(for id: String) -> URL? { packs().first { $0.lastPathComponent == id } }
    public static func defaultPack() -> URL? { packs().first }
}
