import Foundation
import Testing
@testable import StudioKit

private func tempURL(_ name: String) -> URL {
    FileManager.default.temporaryDirectory.appendingPathComponent("\(UUID().uuidString)-\(name)")
}

@Suite struct GeometryTests {
    @Test func triangulationCoversAnLShapedRoom() {
        let l = [Vec2(0, 0), Vec2(4, 0), Vec2(4, 2), Vec2(2, 2), Vec2(2, 4), Vec2(0, 4)]
        let tris = Polygon2.triangulate(l)
        #expect(tris.count == 4)
        let area = tris.reduce(0.0) { $0 + abs(Polygon2.signedArea([l[$1.0], l[$1.1], l[$1.2]])) }
        #expect(abs(area - 12) < 1e-9)
    }

    @Test func selfCrossingOutlineIsNotSimple() {
        #expect(!Polygon2.isSimple([Vec2(0, 0), Vec2(2, 2), Vec2(2, 0), Vec2(0, 2)]))
        #expect(Polygon2.isSimple([Vec2(0, 0), Vec2(2, 0), Vec2(2, 2), Vec2(0, 2)]))
    }

    @Test func splitAndRemoveKeepOneMaterialPerWall() {
        var room = Room.shoebox(width: 4, depth: 4, height: 3, walls: "brick", floor: "carpet", ceiling: "plaster")
        room.splitWall(1, at: Vec2(2, 0))
        #expect(room.outline.count == 5 && room.walls.count == 5)
        let removed = room.removeVertex(2)
        #expect(removed)
        #expect(room.outline.count == 4 && room.walls.count == 4)
    }

    @Test func roomNormalsPointInwardForEitherWinding() {
        for reversed in [false, true] {
            var room = Room.shoebox(width: 4, depth: 6, height: 3, walls: "brick", floor: "carpet", ceiling: "plaster")
            if reversed { room.outline.reverse() }
            for s in room.surfaces() {
                let toCenter = Vec3(0, 1.5, 0) - s.vertices[0]
                #expect((toCenter * s.normal).sum() > 0)
            }
        }
    }

    @Test func meshTrianglesFaceIntoTheRoom() {
        let room = Room.shoebox(width: 4, depth: 6, height: 3, walls: "brick", floor: "carpet", ceiling: "plaster")
        let mesh = ProjectCompiler.mesh(room: room, materials: MaterialLibrary.acoustic)
        #expect(mesh.triangleCount == 12)
        #expect(mesh.triangleMaterials.count == 12)
        func vertex(_ i: UInt32) -> Vec3 {
            let k = Int(i) * 3
            return Vec3(Double(mesh.positions[k]), Double(mesh.positions[k + 1]), Double(mesh.positions[k + 2]))
        }
        for t in 0..<mesh.triangleCount {
            let a = vertex(mesh.indices[3 * t]), b = vertex(mesh.indices[3 * t + 1]), c = vertex(mesh.indices[3 * t + 2])
            let u = b - a, v = c - a
            let n = Vec3(u.y * v.z - u.z * v.y, u.z * v.x - u.x * v.z, u.x * v.y - u.y * v.x)
            #expect(((Vec3(0, 1.5, 0) - a) * n).sum() > 0)
        }
    }
}

@Suite struct AcousticsTests {
    @Test func sabineMatchesHandCalculationWithoutAir() {
        // 10 x 10 x 3: V 300, S 320. All surfaces alpha 0.2 gives A = 64, so T = 0.161 * 300 / 64 = 0.7547 before air.
        var m = MaterialLibrary.acoustic[0]
        m.absorption = Array(repeating: 0.2, count: AcousticMaterial.bandsHz.count)
        let room = Room.shoebox(width: 10, depth: 10, height: 3, walls: m.id, floor: m.id, ceiling: m.id)
        let r = RoomAcoustics.reverb(room: room, materials: [m])
        #expect(abs(r.sabineSeconds[0] - 0.7547) < 0.005)
        #expect(r.eyringSeconds[0] < r.sabineSeconds[0])
        #expect(r.sabineSeconds[7] < r.sabineSeconds[0]) // air takes the top band
    }

    @Test func firstOrderShoeboxHasExactlySixPathsOfImageLength() throws {
        let room = Room.shoebox(width: 4, depth: 6, height: 3, walls: "concrete", floor: "concrete", ceiling: "concrete")
        let source = SoundSource(name: "s", position: Vec3(1, 1, 1))
        let listener = Vec3(-1, 2, -1)
        let paths = ImageSources.paths(room: room, materials: MaterialLibrary.acoustic, source: source, sourcePosition: source.position,
                                       listener: listener, maxOrder: 1)
        #expect(paths.count == 6)
        // Floor bounce: image at y = -1, so length is |(-2, 3, -2)|.
        let floor = try #require(paths.first { $0.surfaces == [.floor] })
        #expect(abs(floor.lengthMetres - (4.0 + 9 + 4).squareRoot()) < 1e-9)
        #expect(abs(floor.points[1].y) < 1e-9)
    }

    @Test func bedOrderComesFromChannelCount() {
        #expect(AmbisonicBed.order(channels: 4) == 1)
        #expect(AmbisonicBed.order(channels: 16) == 3)
        #expect(AmbisonicBed.order(channels: 1) == nil)
        #expect(AmbisonicBed.order(channels: 5) == nil)
    }

    @Test func aBlockOccludesThePathsThroughIt() {
        var room = Room.shoebox(width: 8, depth: 8, height: 3, walls: "concrete", floor: "concrete", ceiling: "concrete")
        #expect(!ImageSources.directIsOccluded(room: room, from: Vec3(-3, 1, 0), to: Vec3(3, 1, 0)))
        room.blocks.append(Block(name: "b", center: Vec3(0, 1, 0), size: Vec3(1, 2.5, 1), material: "concrete"))
        #expect(ImageSources.directIsOccluded(room: room, from: Vec3(-3, 1, 0), to: Vec3(3, 1, 0)))
    }

    @Test func listenerBasisIsOrthonormalAndFacesMinusZ() {
        var l = Listener(position: .zero)
        #expect(abs(l.forward.z + 1) < 1e-12)
        l.yawDegrees = 37; l.pitchDegrees = 21
        #expect(abs((l.forward * l.up).sum()) < 1e-12)
        #expect(abs((l.forward * l.forward).sum() - 1) < 1e-12)
        #expect(abs((l.up * l.up).sum() - 1) < 1e-12)
    }

    @Test func orbitVelocityMatchesCircumference() {
        var p = SpatialProject.starter()
        p.sources[0].motion = .orbit(center: Vec2(0, 0), radius: 2, periodSeconds: 4)
        let v = ProjectCompiler.frame(project: p, time: 1).sources[0].velocity
        #expect(abs((v * v).sum().squareRoot() - 2 * Double.pi * 2 / 4) < 1e-3)
    }
}

@Suite struct DocumentTests {
    @Test func projectRoundTripsThroughJSON() throws {
        var p = SpatialProject.starter()
        p.sources[0].motion = .shuttle(offset: Vec3(1, 2, 3), periodSeconds: 5)
        p.bodies.append(SoundingBody(name: "b", material: "steel", position: Vec3(1, 1, 1), size: Vec3(0.5, 0.3, 0.002)))
        p.strikes.append(Strike(body: p.bodies[0].id, timeSeconds: 1.5, u: 0.2, v: 0.8, impulseNs: 0.1, contactSeconds: 0.002))
        #expect(try ProjectPackage.decode(ProjectPackage.encode(p)) == p)
    }

    @Test func newerFormatIsRefusedRatherThanMisread() throws {
        var p = SpatialProject.starter()
        p.version = SpatialProject.formatVersion + 1
        #expect(throws: ProjectPackage.Failure.self) { try ProjectPackage.decode(ProjectPackage.encode(p)) }
    }

    @Test func packageCarriesItsAudio() async throws {
        let dirA = tempURL("a"), dirB = tempURL("b")
        defer { try? FileManager.default.removeItem(at: dirA); try? FileManager.default.removeItem(at: dirB) }
        let source = ClipStore(directory: dirA)
        let wav = tempURL("clip.wav")
        defer { try? FileManager.default.removeItem(at: wav) }
        let writer = try WAVWriter(url: wav, format: .pcm16, sampleRate: 48000, channels: 1)
        try TestSignal.tone440.generate(sampleRate: 48000, seconds: 0.1).withUnsafeBufferPointer { try writer.append($0) }
        _ = try writer.finish()
        let clip = try await source.importFile(wav)
        #expect(clip.frames == 4800 && clip.channels == 1 && clip.sampleRate == 48000)

        var p = SpatialProject.starter()
        p.clips = [clip]
        let wrapper = try ProjectPackage.write(p, store: source)
        let dest = ClipStore(directory: dirB)
        let read = try ProjectPackage.read(wrapper, store: dest)
        #expect(read.clips == [clip])
        #expect(try await dest.decode(clip).frames == 4800)
    }
}

@Suite struct BodyTests {
    @Test func modelNameFollowsTheLibraryRule() {
        func name(_ s: Vec3, wall: Double = 0) -> String {
            SoundingBody(name: "b", material: "steel", position: .zero, size: s, wallThicknessM: wall).modelName
        }
        #expect(name(Vec3(0.6, 0.4, 0.003)) == "Plate")
        #expect(name(Vec3(0.3, 0.02, 0.02)) == "Bar")
        #expect(name(Vec3(0.4, 0.2, 0.2)) == "Block")
        #expect(name(Vec3(0.4, 0.2, 0.2), wall: 0.002) == "Panel")
    }

    @Test func thinWallCountsForMass() {
        let solid = SoundingBody(name: "b", material: "steel", position: .zero, size: Vec3(1, 1, 0.5))
        let panel = SoundingBody(name: "b", material: "steel", position: .zero, size: Vec3(1, 1, 0.5), wallThicknessM: 0.01)
        #expect(solid.massKg(density: 1000) == 500)
        #expect(panel.massKg(density: 1000) == 10)
    }
}

@Suite struct AudioFileTests {
    /// Read fields straight from the bytes, so the check is independent of the writer.
    private func field32(_ d: Data, _ at: Int) -> UInt32 { d[at..<(at + 4)].enumerated().reduce(0) { $0 | UInt32($1.element) << (8 * UInt32($1.offset)) } }

    @Test func headerAndSizesAreCorrect() throws {
        let url = tempURL("h.wav")
        defer { try? FileManager.default.removeItem(at: url) }
        let w = try WAVWriter(url: url, format: .pcm24, sampleRate: 44100, channels: 2)
        try [Float](repeating: 0.25, count: 200).withUnsafeBufferPointer { try w.append($0) }
        let s = try w.finish()
        let d = try Data(contentsOf: url)
        #expect(s.frames == 100 && s.peak == 0.25 && s.clippedSamples == 0)
        #expect(d.count == 44 + 600 && field32(d, 4) == 36 + 600 && field32(d, 40) == 600 && field32(d, 24) == 44100)
    }

    @Test func fixedPointCountsClipsAndFloatKeepsThem() throws {
        for (format, clipped) in [(WAVFormat.pcm16, Int64(1)), (.float32, 0)] {
            let url = tempURL("c.wav")
            defer { try? FileManager.default.removeItem(at: url) }
            let w = try WAVWriter(url: url, format: format, sampleRate: 48000, channels: 2)
            try [Float](arrayLiteral: 2, 0).withUnsafeBufferPointer { try w.append($0) }
            let s = try w.finish()
            #expect(s.clippedSamples == clipped && s.peak == 2)
        }
    }

    @Test func signalsAreDeterministicAndBounded() {
        for signal in TestSignal.allCases {
            let a = signal.generate(sampleRate: 48000, seconds: 1), b = signal.generate(sampleRate: 48000, seconds: 1)
            #expect(a == b)
            #expect(a.count == 48000)
            #expect(a.allSatisfy { abs($0) <= 1 })
        }
    }
}

/// Records what the exporter asks for. Its output is a ramp, to prove the samples reach the file unchanged; it is a test
/// double for the boundary and says nothing about the sound of the library.
private final class RecordingEngine: SpatialEngine, @unchecked Sendable {
    let status = EngineStatus(name: "Test", version: "0", canPlay: false, canRender: true, canStrike: true, hasLateReverb: false, limitation: nil)
    let preferredSampleRate = 48000.0
    var configured: EngineConfiguration?
    var request: RenderRequest?
    var failAfterBlocks: Int?

    func hrtfSets() async -> [HRTFSetInfo] { [] }
    func importHRTF(_ url: URL) async throws -> HRTFSetInfo { throw CocoaError(.featureUnsupported) }
    func configure(_ c: EngineConfiguration) async throws {}
    func play(from seconds: Double) async throws {}
    func stop() async {}
    func update(_ frame: EngineFrame) {}
    func strike(_ strike: Strike) {}
    func scrape(body: UUID, u: Double, v: Double, speedMetresPerSecond: Double, normalForceN: Double) {}
    func endScrape(body: UUID) {}
    func meters() -> MeterFrame { MeterFrame() }
    func setHeadOrientation(_ quaternion: SIMD4<Float>) {}
    func modes(of body: UUID) async -> [BodyMode] { [] }
    func measuredRoom() async -> MeasuredRoom? { nil }
    func render(_ configuration: EngineConfiguration, _ request: RenderRequest) -> AsyncThrowingStream<RenderBlock, Error> {
        configured = configuration
        self.request = request
        let fail = failAfterBlocks
        return AsyncThrowingStream { c in
            let total = Int64(request.durationSeconds * 48000)
            var done: Int64 = 0, blocks = 0
            while done < total {
                let n = min(Int64(request.blockFrames), total - done)
                c.yield(RenderBlock(samples: [Float](repeating: 0.5, count: Int(n) * 2), framesDone: done + n, framesTotal: total))
                done += n; blocks += 1
                if let fail, blocks >= fail {
                    c.finish(throwing: CocoaError(.fileWriteUnknown))
                    return
                }
            }
            c.finish()
        }
    }
}

@Suite struct ExportTests {
    @Test func exportCarriesTailAndOnlyStrikesInRange() async throws {
        var p = SpatialProject.starter()
        let b = SoundingBody(name: "b", material: "steel", position: .zero, size: Vec3(0.3, 0.02, 0.02))
        p.bodies = [b]
        p.strikes = [Strike(body: b.id, timeSeconds: 0.5, u: 0.5, v: 0.5, impulseNs: 1, contactSeconds: 0.001),
                     Strike(body: b.id, timeSeconds: 5, u: 0.5, v: 0.5, impulseNs: 1, contactSeconds: 0.001)]
        var o = ExportOptions()
        o.durationSeconds = 1; o.tailSeconds = 0.5
        let url = tempURL("x.wav")
        defer { try? FileManager.default.removeItem(at: url) }
        let engine = RecordingEngine()
        let r = try await Exporter.export(project: p, engine: engine, store: ClipStore(directory: tempURL("s")), options: o, to: url) { _ in }
        #expect(r.summary.frames == 72000 && r.summary.peak == 0.5)
        #expect(engine.request?.strikes.count == 1)
        #expect(engine.configured?.bodies.count == 1)
    }

    @Test func failedExportLeavesNoFile() async throws {
        let url = tempURL("f.wav")
        let engine = RecordingEngine()
        engine.failAfterBlocks = 2
        await #expect(throws: (any Error).self) {
            try await Exporter.export(project: .starter(), engine: engine, store: ClipStore(directory: tempURL("s")), options: ExportOptions(), to: url) { _ in }
        }
        #expect(!FileManager.default.fileExists(atPath: url.path))
    }

    @Test func unavailableEngineRefusesToRender() async {
        let url = tempURL("u.wav")
        await #expect(throws: UnavailableEngine.Missing.self) {
            try await Exporter.export(project: .starter(), engine: UnavailableEngine(reason: "off"), store: ClipStore(directory: tempURL("s")), options: ExportOptions(), to: url) { _ in }
        }
        #expect(!FileManager.default.fileExists(atPath: url.path))
    }
}
