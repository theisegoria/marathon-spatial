import Foundation

/// The studio's whole view of the Marathon Spatial library. All sound crosses this boundary: nothing else in the app
/// computes audio, so what is auditioned and exported is what the library does.
public protocol SpatialEngine: AnyObject, Sendable {
    var status: EngineStatus { get }
    /// The rate to configure at; for live output this is the output device's rate.
    var preferredSampleRate: Double { get }
    func hrtfSets() async -> [HRTFSetInfo]
    /// Convert a user's SOFA file into a set stored with the app.
    func importHRTF(_ url: URL) async throws -> HRTFSetInfo
    /// Apply what is not per-block state. May allocate and take time; never called from a render thread.
    func configure(_ configuration: EngineConfiguration) async throws
    /// Start live output on the current device at a timeline position.
    func play(from seconds: Double) async throws
    func stop() async
    /// Publish the newest pose. Must not block.
    func update(_ frame: EngineFrame)
    /// Head orientation relative to the listener pose, as x, y, z, w.
    func setHeadOrientation(_ quaternion: SIMD4<Float>)
    func strike(_ strike: Strike)
    /// A sustained rub across a body, updated as the pointer moves. Live only.
    func scrape(body: UUID, u: Double, v: Double, speedMetresPerSecond: Double, normalForceN: Double)
    func endScrape(body: UUID)
    func meters() -> MeterFrame
    /// The modes the library computed for a body, ascending frequency.
    func modes(of body: UUID) async -> [BodyMode]
    /// The room as the library's simulator sees it, or nil before a scene is configured.
    func measuredRoom() async -> MeasuredRoom?
    /// Render faster than real time into memory. No device is opened and live state is not touched.
    func render(_ configuration: EngineConfiguration, _ request: RenderRequest) -> AsyncThrowingStream<RenderBlock, Error>
}

public struct EngineStatus: Sendable, Equatable {
    public var name: String
    public var version: String
    public var canPlay: Bool
    public var canRender: Bool
    public var canStrike: Bool
    /// The library's late reverb is reachable.
    public var hasLateReverb: Bool
    /// Why the engine cannot do its work, in a sentence for the interface.
    public var limitation: String?

    public init(name: String, version: String, canPlay: Bool, canRender: Bool, canStrike: Bool, hasLateReverb: Bool, limitation: String?) {
        self.name = name; self.version = version; self.canPlay = canPlay
        self.canRender = canRender; self.canStrike = canStrike; self.hasLateReverb = hasLateReverb; self.limitation = limitation
    }
}

public struct HRTFSetInfo: Sendable, Hashable, Identifiable {
    public var id: String
    public var name: String
    public var detail: String
    public var licence: String
    public init(id: String, name: String, detail: String, licence: String) {
        self.id = id; self.name = name; self.detail = detail; self.licence = licence
    }
}

public struct BodyMode: Sendable, Hashable, Identifiable {
    public var frequencyHz: Double
    public var dampingRatio: Double
    public var id: Double { frequencyHz }
    public var decaySeconds: Double {
        let rate = dampingRatio * 2 * .pi * frequencyHz
        return rate > 0 ? log(1000.0) / rate : .infinity
    }
    public init(frequencyHz: Double, dampingRatio: Double) { self.frequencyHz = frequencyHz; self.dampingRatio = dampingRatio }
}

/// What the library's simulator reports for the room.
public struct MeasuredRoom: Sendable, Equatable {
    public var volume: Double
    public var surfaceArea: Double
    /// One per `AcousticMaterial.bandsHz`.
    public var t60Seconds: [Double]
    public init(volume: Double, surfaceArea: Double, t60Seconds: [Double]) {
        self.volume = volume; self.surfaceArea = surfaceArea; self.t60Seconds = t60Seconds
    }
    public var midSeconds: Double { (t60Seconds[3] + t60Seconds[4]) / 2 }
}

/// Triangles and materials as flat arrays, ready to pass through a C boundary.
public struct EngineMesh: Sendable, Equatable {
    public var positions: [Float] = []
    public var indices: [UInt32] = []
    /// One per triangle, indexing `EngineConfiguration.materials`.
    public var triangleMaterials: [UInt32] = []
    public var triangleCount: Int { indices.count / 3 }
}

public struct EngineSourceSetup: Sendable {
    public var source: SoundSource
    /// PCM to play. nil for a source with no content, which stays silent.
    public var clip: DecodedClip?
}

public struct EngineBedSetup: Sendable {
    public var bed: AmbisonicBed
    /// All channels, in ACN order.
    public var clip: DecodedClip
}

public struct EngineBodySetup: Sendable {
    public var body: SoundingBody
    public var material: PhysicalMaterial
}

public struct EngineConfiguration: Sendable {
    public var room: Room
    public var mesh: EngineMesh
    public var materials: [AcousticMaterial]
    public var sources: [EngineSourceSetup]
    public var bodies: [EngineBodySetup]
    public var beds: [EngineBedSetup]
    public var render: RenderSettings
    public var sampleRate: Double
}

public struct EngineSourcePose: Sendable, Equatable {
    public var id: UUID
    public var position: Vec3
    public var velocity: Vec3
}

/// Everything that moves, at one instant.
public struct EngineFrame: Sendable, Equatable {
    public var timeSeconds: Double
    public var listenerPosition: Vec3
    public var listenerVelocity: Vec3
    public var listenerForward: Vec3
    public var listenerUp: Vec3
    public var sources: [EngineSourcePose]
    public var bodies: [EngineSourcePose]
}

public struct MeterFrame: Sendable, Equatable {
    public var peak: [Float] = [0, 0]
    public var rms: [Float] = [0, 0]
    /// Linear level per source after the mix.
    public var sources: [UUID: Float] = [:]
    /// Output samples over full scale since playback began.
    public var overs: Int = 0
    public init() {}
}

public struct RenderRequest: Sendable {
    public var startSeconds: Double
    public var durationSeconds: Double
    public var blockFrames: Int
    /// Pose at any time in the render; deterministic, so an export is repeatable.
    public var frame: @Sendable (Double) -> EngineFrame
    /// Strikes inside the range, in time order.
    public var strikes: [Strike]

    public init(startSeconds: Double, durationSeconds: Double, blockFrames: Int, frame: @escaping @Sendable (Double) -> EngineFrame, strikes: [Strike]) {
        self.startSeconds = startSeconds; self.durationSeconds = durationSeconds; self.blockFrames = blockFrames
        self.frame = frame; self.strikes = strikes
    }
}

public struct RenderBlock: Sendable {
    /// Interleaved stereo.
    public var samples: [Float]
    public var framesDone: Int64
    public var framesTotal: Int64

    public init(samples: [Float], framesDone: Int64, framesTotal: Int64) {
        self.samples = samples; self.framesDone = framesDone; self.framesTotal = framesTotal
    }
}

/// Stands in when the library is not linked. It makes no sound, and says so.
public final class UnavailableEngine: SpatialEngine {
    public struct Missing: Error, LocalizedError {
        public var reason: String
        public var errorDescription: String? { reason }
        public init(reason: String) { self.reason = reason }
    }

    public let status: EngineStatus
    private let reason: String

    public init(reason: String) {
        self.reason = reason
        status = EngineStatus(name: "Not linked", version: "", canPlay: false, canRender: false, canStrike: false, hasLateReverb: false, limitation: reason)
    }

    public var preferredSampleRate: Double { 48000 }
    public func hrtfSets() async -> [HRTFSetInfo] { [] }
    public func importHRTF(_ url: URL) async throws -> HRTFSetInfo { throw Missing(reason: reason) }
    public func configure(_ configuration: EngineConfiguration) async throws {}
    public func play(from seconds: Double) async throws { throw Missing(reason: reason) }
    public func stop() async {}
    public func update(_ frame: EngineFrame) {}
    public func setHeadOrientation(_ quaternion: SIMD4<Float>) {}
    public func strike(_ strike: Strike) {}
    public func scrape(body: UUID, u: Double, v: Double, speedMetresPerSecond: Double, normalForceN: Double) {}
    public func endScrape(body: UUID) {}
    public func meters() -> MeterFrame { MeterFrame() }
    public func modes(of body: UUID) async -> [BodyMode] { [] }
    public func measuredRoom() async -> MeasuredRoom? { nil }
    public func render(_ configuration: EngineConfiguration, _ request: RenderRequest) -> AsyncThrowingStream<RenderBlock, Error> {
        AsyncThrowingStream { $0.finish(throwing: Missing(reason: reason)) }
    }
}
