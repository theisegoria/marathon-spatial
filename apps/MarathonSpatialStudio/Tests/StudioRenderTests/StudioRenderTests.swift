import Foundation
import StudioKit
import StudioRender
import Testing

// Everything here renders into memory through the library. This target links no output library, so it cannot open a device.
// Each test compares two renders of one scene, so it pins what the wiring does and not how anything sounds.

private let rate = 48000.0

private func store() -> ClipStore {
    ClipStore(directory: FileManager.default.temporaryDirectory.appendingPathComponent("studio-render-tests-\(UUID().uuidString)"))
}

private func render(_ project: SpatialProject, seconds: Double) async throws -> [Float] {
    let configuration = try await ProjectCompiler.configuration(project: project, store: store(), sampleRate: rate)
    let snapshot = project
    let request = RenderRequest(startSeconds: 0, durationSeconds: seconds, blockFrames: 1024,
                                frame: { ProjectCompiler.frame(project: snapshot, time: $0) },
                                strikes: project.strikes.sorted { $0.timeSeconds < $1.timeSeconds })
    var out: [Float] = []
    try OfflineRender.run(configuration, request) { out += $0.samples }
    return out
}

private func energy(_ stereo: [Float], from: Double, to: Double) -> Double {
    let a = Int(from * rate) * 2, b = min(Int(to * rate) * 2, stereo.count)
    return stride(from: a, to: b, by: 1).reduce(0.0) { $0 + Double(stereo[$1]) * Double(stereo[$1]) }
}

private func firstFrame(_ stereo: [Float], above level: Float) -> Int? {
    stereo.firstIndex { abs($0) > level }.map { $0 / 2 }
}

/// One click from a corner of a hard room.
private func room(width: Double = 12, depth: Double = 16, height: Double = 6) -> SpatialProject {
    var p = SpatialProject(room: .shoebox(width: width, depth: depth, height: height, walls: "concrete", floor: "concrete", ceiling: "concrete"),
                           listener: Listener(position: Vec3(0, 1.6, 4)))
    var s = SoundSource(name: "click", position: Vec3(2, 1.6, -3), content: .signal(.clickTrain))
    s.loops = false
    s.doppler = false
    p.sources.append(s)
    p.render.limiter = false
    return p
}

@Test func lateReverbAddsATailAndTheTrimScalesIt() async throws {
    var dry = room()
    dry.render.lateReverb = false
    dry.render.earlyReflections = false
    var wet = dry
    wet.render.lateReverb = true
    var quiet = wet
    quiet.render.reverbTrimDB = -12

    let a = try await render(dry, seconds: 0.45), b = try await render(wet, seconds: 0.45), c = try await render(quiet, seconds: 0.45)
    // The first click is over well before this window, and the second has not reached it.
    let tailDry = energy(a, from: 0.2, to: 0.24), tailWet = energy(b, from: 0.2, to: 0.24), tailQuiet = energy(c, from: 0.2, to: 0.24)
    #expect(tailDry < 1e-9)
    #expect(tailWet > 1e-6)
    #expect(abs(10 * log10(tailQuiet / tailWet) + 12) < 1)
}

@Test func reflectionsCountFromTheDirectSoundAtAnyDistance() async throws {
    // 34 m apart. A delay counted from emission would put the first reflection about 100 ms behind the direct sound.
    var off = SpatialProject(room: .shoebox(width: 5, depth: 60, height: 4, walls: "concrete", floor: "concrete", ceiling: "concrete"),
                             listener: Listener(position: Vec3(0, 1.6, 17)))
    var s = SoundSource(name: "click", position: Vec3(0, 1.6, -17), content: .signal(.clickTrain))
    s.loops = false
    s.doppler = false
    s.airAbsorption = false
    off.sources.append(s)
    off.render.limiter = false
    off.render.lateReverb = false
    off.render.earlyReflections = false
    var on = off
    on.render.earlyReflections = true

    let dry = try await render(off, seconds: 0.2), wet = try await render(on, seconds: 0.2)
    let direct = try #require(firstFrame(dry, above: 1e-4))
    let added = zip(wet, dry).map { $0 - $1 }
    let first = try #require(firstFrame(added, above: 1e-5))
    #expect(direct > Int(0.09 * rate))
    #expect(first >= direct - 256)
    #expect(Double(first - direct) / rate < 0.02)
}

@Test func aBlockInTheWayQuietensTheDirectSound() async throws {
    var clear = room()
    clear.render.lateReverb = false
    clear.render.earlyReflections = false
    var blocked = clear
    blocked.room.blocks.append(Block(name: "wall", center: Vec3(1, 1.6, 0.5), size: Vec3(6, 3, 0.4), material: "concrete"))

    let a = try await render(clear, seconds: 0.1), b = try await render(blocked, seconds: 0.1)
    #expect(energy(b, from: 0, to: 0.1) < 0.5 * energy(a, from: 0, to: 0.1))
}

@Test func aStruckBodySoundsAtItsTimeAndExcitesTheRoom() async throws {
    var dry = room()
    dry.sources = []
    let body = SoundingBody(name: "plate", material: "steel", position: Vec3(1, 1.2, 1), size: Vec3(0.6, 0.4, 0.004))
    dry.bodies = [body]
    dry.strikes = [Strike(body: body.id, timeSeconds: 0.1, u: 0.3, v: 0.4, impulseNs: 1, contactSeconds: 1e-4)]
    dry.render.lateReverb = false
    dry.render.earlyReflections = false
    var wet = dry
    wet.render.lateReverb = true
    wet.render.earlyReflections = true

    let a = try await render(dry, seconds: 0.5), b = try await render(wet, seconds: 0.5)
    #expect(energy(a, from: 0, to: 0.09) == 0)
    #expect(energy(a, from: 0.1, to: 0.2) > 1e-6)
    #expect(energy(b, from: 0.2, to: 0.5) > energy(a, from: 0.2, to: 0.5))
}

/// To -120 dB, not to the bit: two runs in one process have differed by about 1e-7 at most, with the room's simulator off as well.
@Test func anExportRepeats() async throws {
    let p = room()
    let a = try await render(p, seconds: 0.5), b = try await render(p, seconds: 0.5)
    #expect(a.count == b.count)
    #expect((zip(a, b).map { abs($0 - $1) }.max() ?? 1) < 1e-6)
}

@Test func theSimulatorSeesTheDocumentsRoom() async throws {
    let p = room(width: 10, depth: 8, height: 4)
    let configuration = try await ProjectCompiler.configuration(project: p, store: store(), sampleRate: rate)
    let instance = try EngineInstance(rate: rate, limiter: false)
    try instance.apply(configuration)
    try instance.tickAcoustics(ProjectCompiler.frame(project: p, time: 0), seconds: 1.0 / 60)
    let measured = try #require(instance.measuredRoom())
    #expect(abs(measured.volume - 320) < 1)
    #expect(measured.t60Seconds.count == AcousticMaterial.bandsHz.count)
    let estimate = RoomAcoustics.reverb(room: p.room, materials: p.materials)
    #expect(measured.midSeconds > 0.5 * estimate.midSeconds && measured.midSeconds < 2 * estimate.midSeconds)
}
