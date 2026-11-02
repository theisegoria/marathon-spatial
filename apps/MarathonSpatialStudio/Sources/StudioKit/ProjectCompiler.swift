import Foundation
import simd

/// Turns a document into what the engine consumes.
public enum ProjectCompiler {
    /// Test signals loop, so a few seconds are enough; the sweep is long enough to hear each decade.
    static let signalSeconds = 8.0

    /// Triangulate the room and blocks. Winding follows each surface's normal, so back-face tests hold.
    public static func mesh(room: Room, materials: [AcousticMaterial]) -> EngineMesh {
        var mesh = EngineMesh()
        func materialIndex(_ id: MaterialID) -> UInt32 { UInt32(materials.firstIndex { $0.id == id } ?? 0) }
        func vertex(_ p: Vec3) -> UInt32 {
            let i = UInt32(mesh.positions.count / 3)
            mesh.positions.append(contentsOf: [Float(p.x), Float(p.y), Float(p.z)])
            return i
        }
        func triangle(_ a: UInt32, _ b: UInt32, _ c: UInt32, _ normal: Vec3, _ p: [Vec3], _ m: UInt32) {
            let n = simd_cross(p[1] - p[0], p[2] - p[0])
            mesh.indices.append(contentsOf: dot(n, normal) >= 0 ? [a, b, c] : [a, c, b])
            mesh.triangleMaterials.append(m)
        }
        for s in room.surfaces() {
            let m = materialIndex(s.material)
            let ids = s.vertices.map(vertex)
            if s.vertices.count == 4, s.ref != .floor, s.ref != .ceiling {
                triangle(ids[0], ids[1], ids[2], s.normal, [s.vertices[0], s.vertices[1], s.vertices[2]], m)
                triangle(ids[0], ids[2], ids[3], s.normal, [s.vertices[0], s.vertices[2], s.vertices[3]], m)
            } else {
                for t in Polygon2.triangulate(room.outline) {
                    triangle(ids[t.0], ids[t.1], ids[t.2], s.normal, [s.vertices[t.0], s.vertices[t.1], s.vertices[t.2]], m)
                }
            }
        }
        return mesh
    }

    /// Everything the engine needs for a project, decoding clips and generating signals. Off the main actor.
    public static func configuration(project: SpatialProject, store: ClipStore, sampleRate: Double) async throws -> EngineConfiguration {
        var decoded: [UUID: DecodedClip] = [:]
        var setups: [EngineSourceSetup] = []
        for source in project.audibleSources {
            var clip: DecodedClip?
            switch source.content {
            case .none: break
            case .signal(let signal):
                clip = DecodedClip(channels: [signal.generate(sampleRate: sampleRate, seconds: signalSeconds)], sampleRate: sampleRate)
            case .clip(let id):
                if let cached = decoded[id] {
                    clip = cached
                } else if let ref = project.clip(id) {
                    let d = try await store.decode(ref)
                    decoded[id] = d
                    clip = d
                }
            }
            setups.append(EngineSourceSetup(source: source, clip: clip))
        }
        let bodies = project.bodies.compactMap { body in
            project.physicalMaterials.first { $0.id == body.material }.map { EngineBodySetup(body: body, material: $0) }
        }
        var beds: [EngineBedSetup] = []
        for bed in project.beds {
            guard let ref = project.clip(bed.clip), AmbisonicBed.order(channels: ref.channels) != nil else { continue }
            let d: DecodedClip
            if let cached = decoded[bed.clip] {
                d = cached
            } else {
                d = try await store.decode(ref)
                decoded[bed.clip] = d
            }
            beds.append(EngineBedSetup(bed: bed, clip: d))
        }
        return EngineConfiguration(room: project.room, mesh: mesh(room: project.room, materials: project.materials), materials: project.materials,
                                   sources: setups, bodies: bodies, beds: beds, render: project.render, sampleRate: sampleRate)
    }

    /// Poses at `time`, with velocities by central difference so scripted motion drives Doppler. Pure: export calls it per block.
    public static func frame(project: SpatialProject, time: Double) -> EngineFrame {
        let h = 0.005
        func pose(_ id: UUID, _ home: Vec3, _ motion: Motion) -> EngineSourcePose {
            let p = motion.position(home: home, time: time)
            let v = (motion.position(home: home, time: time + h) - motion.position(home: home, time: time - h)) / (2 * h)
            return EngineSourcePose(id: id, position: p, velocity: v)
        }
        let l = project.listener
        let lp = pose(UUID(uuid: UUID_NULL), l.position, l.motion)
        return EngineFrame(timeSeconds: time, listenerPosition: lp.position, listenerVelocity: lp.velocity,
                           listenerForward: l.forward, listenerUp: l.up,
                           sources: project.audibleSources.map { pose($0.id, $0.position, $0.motion) },
                           bodies: project.bodies.map { pose($0.id, $0.position, .still) })
    }
}

private let UUID_NULL: uuid_t = (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)

public struct ExportOptions: Sendable, Equatable {
    public var format: WAVFormat = .pcm24
    public var sampleRate = 48000
    public var startSeconds = 0.0
    public var durationSeconds = 30.0
    /// Rendered after the end so reverb and modes ring out.
    public var tailSeconds = 2.0
    public init() {}
}

public struct ExportResult: Sendable, Equatable {
    public var summary: WAVSummary
    public var seconds: Double
    public var url: URL
}

public enum Exporter {
    /// Render through the engine into a WAV. A cancelled or failed export leaves no partial file behind.
    public static func export(project: SpatialProject, engine: SpatialEngine, store: ClipStore, options: ExportOptions, to url: URL,
                              progress: @Sendable (Double) -> Void) async throws -> ExportResult {
        let configuration = try await ProjectCompiler.configuration(project: project, store: store, sampleRate: Double(options.sampleRate))
        let range = options.startSeconds..<(options.startSeconds + options.durationSeconds)
        let snapshot = project
        let request = RenderRequest(startSeconds: options.startSeconds, durationSeconds: options.durationSeconds + options.tailSeconds,
                                    blockFrames: 1024,
                                    frame: { ProjectCompiler.frame(project: snapshot, time: $0) },
                                    strikes: project.strikes.filter { range.contains($0.timeSeconds) }.sorted { $0.timeSeconds < $1.timeSeconds })
        let writer = try WAVWriter(url: url, format: options.format, sampleRate: options.sampleRate, channels: 2)
        do {
            for try await block in engine.render(configuration, request) {
                try Task.checkCancellation()
                try block.samples.withUnsafeBufferPointer { try writer.append($0) }
                progress(block.framesTotal > 0 ? Double(block.framesDone) / Double(block.framesTotal) : 0)
            }
            let summary = try writer.finish()
            return ExportResult(summary: summary, seconds: Double(summary.frames) / Double(options.sampleRate), url: url)
        } catch {
            _ = try? writer.finish()
            try? FileManager.default.removeItem(at: url)
            throw error
        }
    }
}
