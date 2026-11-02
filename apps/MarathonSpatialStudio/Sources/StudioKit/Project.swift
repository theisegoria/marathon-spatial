import Foundation

/// Identifies a material in `SpatialProject.materials`.
public typealias MaterialID = String

/// Frequency-dependent acoustic data, passed to the library as plain data.
public struct AcousticMaterial: Codable, Hashable, Sendable, Identifiable {
    /// Octave band centres of `absorption`, in Hz.
    public static let bandsHz: [Double] = [62.5, 125, 250, 500, 1000, 2000, 4000, 8000]

    public var id: MaterialID
    public var name: String
    public var absorption: [Double]
    public var scattering: Double
    /// Energy through one crossing of the surface, in dB; the library takes it as a linear energy fraction per band.
    public var transmissionDB: Double
    /// Palette index, so that the swatch colour lives in the UI and the file stays plain data.
    public var swatch: Int

    public init(id: MaterialID, name: String, absorption: [Double], scattering: Double, transmissionDB: Double, swatch: Int) {
        self.id = id; self.name = name; self.absorption = absorption
        self.scattering = scattering; self.transmissionDB = transmissionDB; self.swatch = swatch
    }

    /// Mean of the mid bands (500 Hz to 2 kHz), for captions and quick comparison.
    public var midAbsorption: Double { (absorption[3] + absorption[4] + absorption[5]) / 3 }
}

public struct Block: Codable, Hashable, Sendable, Identifiable {
    public var id = UUID()
    public var name: String
    public var center: Vec3
    public var size: Vec3
    public var material: MaterialID
    public init(name: String, center: Vec3, size: Vec3, material: MaterialID) {
        self.name = name; self.center = center; self.size = size; self.material = material
    }
}

/// An extruded polygon room: `outline` in plan, floor at y = 0, ceiling at y = height.
public struct Room: Codable, Hashable, Sendable {
    public var outline: [Vec2]
    public var height: Double
    /// One material per outline edge; edge i runs from outline[i] to outline[i + 1].
    public var walls: [MaterialID]
    public var floor: MaterialID
    public var ceiling: MaterialID
    public var blocks: [Block]

    public init(outline: [Vec2], height: Double, walls: [MaterialID], floor: MaterialID, ceiling: MaterialID, blocks: [Block] = []) {
        self.outline = outline; self.height = height; self.walls = walls
        self.floor = floor; self.ceiling = ceiling; self.blocks = blocks
    }

    public static func shoebox(width: Double, depth: Double, height: Double, walls: MaterialID, floor: MaterialID, ceiling: MaterialID) -> Room {
        let w = width / 2, d = depth / 2
        return Room(outline: [Vec2(-w, -d), Vec2(w, -d), Vec2(w, d), Vec2(-w, d)], height: height,
                    walls: Array(repeating: walls, count: 4), floor: floor, ceiling: ceiling)
    }

    public var centroid: Vec2 {
        guard !outline.isEmpty else { return Vec2(0, 0) }
        var s = Vec2(0, 0)
        for p in outline { s = s + p }
        return s * (1 / Double(outline.count))
    }

    public func contains(_ p: Vec3) -> Bool {
        p.y >= 0 && p.y <= height && Polygon2.contains(outline, p.plan)
    }

    /// Split edge `index` at `point`; both halves keep the material. Returns the new vertex index.
    @discardableResult public mutating func splitWall(_ index: Int, at point: Vec2) -> Int {
        outline.insert(point, at: index + 1)
        walls.insert(walls[index], at: index + 1)
        return index + 1
    }

    /// Remove a vertex, merging its two walls, if the room stays a simple polygon.
    @discardableResult public mutating func removeVertex(_ index: Int) -> Bool {
        guard outline.count > 3, outline.indices.contains(index) else { return false }
        var next = outline
        next.remove(at: index)
        guard Polygon2.isSimple(next) else { return false }
        outline = next
        walls.remove(at: index)
        return true
    }
}

public struct Directivity: Codable, Hashable, Sendable {
    /// 0 is omnidirectional, 0.5 cardioid, 1 figure of eight.
    public var dipoleWeight: Double = 0
    /// Sharpens the pattern; 1 is the plain pattern.
    public var power: Double = 1
    public var yawDegrees: Double = 0
    public var pitchDegrees: Double = 0
    public init() {}

    /// Linear gain toward a direction given in degrees off the source's axis.
    public func gain(offAxisDegrees a: Double) -> Double {
        let c = cos(a * .pi / 180)
        return pow(abs((1 - dipoleWeight) + dipoleWeight * c), max(power, 0.01))
    }
}

/// Scripted movement, evaluated identically for live playback and offline export.
public enum Motion: Codable, Hashable, Sendable {
    case still
    case orbit(center: Vec2, radius: Double, periodSeconds: Double)
    /// Moves to `offset` and back, once per period.
    case shuttle(offset: Vec3, periodSeconds: Double)

    /// Position at `time` for an object at rest at `home`.
    public func position(home: Vec3, time: Double) -> Vec3 {
        switch self {
        case .still: return home
        case let .orbit(center, radius, period):
            guard period != 0 else { return home }
            let a = 2 * Double.pi * time / period
            return Vec3(center.x + radius * cos(a), home.y, center.z + radius * sin(a))
        case let .shuttle(offset, period):
            guard period != 0 else { return home }
            return home + offset * (0.5 - 0.5 * cos(2 * Double.pi * time / period))
        }
    }
}

/// Procedural signals, generated locally, so a project needs no imported audio to be auditioned.
public enum TestSignal: String, Codable, Hashable, Sendable, CaseIterable, Identifiable {
    case pinkNoise, whiteNoise, noiseBursts, clickTrain, sineSweep, tone440
    public var id: String { rawValue }
    public var title: String {
        switch self {
        case .pinkNoise: "Pink noise"
        case .whiteNoise: "White noise"
        case .noiseBursts: "Noise bursts"
        case .clickTrain: "Click train"
        case .sineSweep: "Sine sweep"
        case .tone440: "Tone 440 Hz"
        }
    }
}

public enum SourceContent: Codable, Hashable, Sendable {
    case none
    case clip(UUID)
    case signal(TestSignal)
}

public struct SoundSource: Codable, Hashable, Sendable, Identifiable {
    public var id = UUID()
    public var name: String
    public var position: Vec3
    public var content: SourceContent = .none
    public var gainDB: Double = 0
    public var startSeconds: Double = 0
    public var loops = true
    public var muted = false
    public var solo = false
    public var directivity = Directivity()
    public var motion: Motion = .still
    public var doppler = true
    public var airAbsorption = true
    public var nearField = true
    public var reflections = true

    public init(name: String, position: Vec3, content: SourceContent = .none) {
        self.name = name; self.position = position; self.content = content
    }
}

public struct Listener: Codable, Hashable, Sendable {
    public var position: Vec3
    /// Heading in degrees: 0 faces -Z, positive turns toward +X.
    public var yawDegrees: Double = 0
    public var pitchDegrees: Double = 0
    public var motion: Motion = .still

    public init(position: Vec3) { self.position = position }

    public var forward: Vec3 {
        let y = yawDegrees * .pi / 180, p = pitchDegrees * .pi / 180
        return Vec3(sin(y) * cos(p), sin(p), -cos(y) * cos(p))
    }

    public var up: Vec3 {
        let y = yawDegrees * .pi / 180, p = pitchDegrees * .pi / 180
        return Vec3(-sin(y) * sin(p), cos(p), cos(y) * sin(p))
    }
}

/// Elastic constants for modal synthesis: plain data for the library.
public struct PhysicalMaterial: Codable, Hashable, Sendable, Identifiable {
    /// A library preset name. The library fills in contact properties from it; these fields override the elastic ones.
    public var id: String
    public var name: String
    public var densityKgM3: Double
    public var youngsModulusPa: Double
    public var poissonRatio: Double
    /// Structural loss factor: eta = 2 * xi. Larger rings shorter.
    public var lossFactor: Double

    public init(id: String, name: String, densityKgM3: Double, youngsModulusPa: Double, poissonRatio: Double, lossFactor: Double) {
        self.id = id; self.name = name; self.densityKgM3 = densityKgM3; self.youngsModulusPa = youngsModulusPa
        self.poissonRatio = poissonRatio; self.lossFactor = lossFactor
    }
}

public struct SoundingBody: Codable, Hashable, Sendable, Identifiable {
    public var id = UUID()
    public var name: String
    public var material: String
    public var position: Vec3
    /// Longest, middle and thickness in metres, in that order.
    public var size: Vec3
    /// Above 0, the box is a thin-walled panel of this thickness.
    public var wallThicknessM: Double = 0
    public var gainDB: Double = 0

    public init(name: String, material: String, position: Vec3, size: Vec3, wallThicknessM: Double = 0) {
        self.name = name; self.material = material; self.position = position; self.size = size; self.wallThicknessM = wallThicknessM
    }

    /// What the library models this box as, by the rule of its proportions.
    public var modelName: String {
        let s = sortedSize
        if wallThicknessM > 0, wallThicknessM < s.z { return "Panel" }
        if s.z < 0.2 * s.y { return "Plate" }
        return s.y < 0.2 * s.x ? "Bar" : "Block"
    }

    /// Size with its axes sorted longest first, as the library requires.
    public var sortedSize: Vec3 {
        let v = [size.x, size.y, size.z].sorted(by: >)
        return Vec3(v[0], v[1], v[2])
    }

    public func massKg(density: Double) -> Double {
        let s = sortedSize
        return density * s.x * s.y * (wallThicknessM > 0 ? min(wallThicknessM, s.z) : s.z)
    }

    /// Local strike position in metres from u, v in 0...1 across the broad face.
    public func localPosition(u: Double, v: Double) -> Vec3 {
        let s = sortedSize
        return Vec3((u - 0.5) * s.x, (v - 0.5) * s.y, s.z / 2)
    }
}

/// One hit on a body. `u` and `v` are in 0...1 across its striking face.
public struct Strike: Codable, Hashable, Sendable, Identifiable {
    public static let contactRange = 1e-5...0.02

    public var id = UUID()
    public var body: UUID
    public var timeSeconds: Double
    public var u: Double
    public var v: Double
    public var impulseNs: Double
    /// Half-sine contact duration: short is a hard mallet and excites high modes, long is a soft one.
    public var contactSeconds: Double

    public init(body: UUID, timeSeconds: Double, u: Double, v: Double, impulseNs: Double, contactSeconds: Double) {
        self.body = body; self.timeSeconds = timeSeconds; self.u = u; self.v = v
        self.impulseNs = impulseNs; self.contactSeconds = contactSeconds
    }
}

/// Higher-order ambisonic background, decoded through the same HRTFs. The clip holds (order + 1)^2 ACN/SN3D channels.
public struct AmbisonicBed: Codable, Hashable, Sendable, Identifiable {
    public var id = UUID()
    public var name: String
    public var clip: UUID
    public var gainDB: Double = 0
    public var yawDegrees: Double = 0
    public var loops = true
    /// Fixed to the head rather than the world.
    public var headLocked = false

    public init(name: String, clip: UUID) { self.name = name; self.clip = clip }

    /// The order a channel count represents, when it is a whole one.
    public static func order(channels: Int) -> Int? {
        let root = Int(Double(channels).squareRoot().rounded())
        return root >= 2 && root * root == channels ? root - 1 : nil
    }
}

public struct ClipReference: Codable, Hashable, Sendable, Identifiable {
    public var id = UUID()
    public var name: String
    /// SHA-256 of the file bytes: its name in the document and in the clip cache.
    public var sha256: String
    public var fileExtension: String
    public var frames: Int64
    public var sampleRate: Double
    public var channels: Int

    public init(name: String, sha256: String, fileExtension: String, frames: Int64, sampleRate: Double, channels: Int) {
        self.name = name; self.sha256 = sha256; self.fileExtension = fileExtension
        self.frames = frames; self.sampleRate = sampleRate; self.channels = channels
    }

    public var durationSeconds: Double { sampleRate > 0 ? Double(frames) / sampleRate : 0 }
    public var storedName: String { fileExtension.isEmpty ? sha256 : "\(sha256).\(fileExtension)" }
}

public enum OutputMode: String, Codable, Hashable, Sendable, CaseIterable, Identifiable {
    case binaural, stereo
    public var id: String { rawValue }
    public var title: String { self == .binaural ? "Binaural (headphones)" : "Stereo (speakers)" }
}

public struct RenderSettings: Codable, Hashable, Sendable {
    public var output: OutputMode = .binaural
    /// nil selects the library's default set.
    public var hrtfSet: String?
    /// 1 to 3.
    public var reflectionOrder = 3
    public var earlyReflections = true
    public var lateReverb = true
    public var reverbTrimDB = 0.0
    public var masterGainDB = 0.0
    /// A master peak limiter at 0 dBFS. Off keeps levels linear for measurement.
    public var limiter = true
    public init() {}
}

public struct SpatialProject: Codable, Hashable, Sendable {
    public static let formatVersion = 1

    public var version = SpatialProject.formatVersion
    public var room: Room
    public var materials: [AcousticMaterial]
    public var physicalMaterials: [PhysicalMaterial]
    public var listener: Listener
    public var sources: [SoundSource] = []
    public var bodies: [SoundingBody] = []
    /// Strikes recorded into the timeline. Export plays them.
    public var strikes: [Strike] = []
    public var clips: [ClipReference] = []
    public var beds: [AmbisonicBed] = []
    public var render = RenderSettings()
    public var timelineSeconds = 30.0

    public init(room: Room, materials: [AcousticMaterial] = MaterialLibrary.acoustic,
                physicalMaterials: [PhysicalMaterial] = MaterialLibrary.physical, listener: Listener) {
        self.room = room; self.materials = materials; self.physicalMaterials = physicalMaterials; self.listener = listener
    }

    public func material(_ id: MaterialID) -> AcousticMaterial? { materials.first { $0.id == id } }
    public func clip(_ id: UUID) -> ClipReference? { clips.first { $0.id == id } }

    /// Sources that sound: muted ones are out and, if any is soloed, only soloed ones are in.
    public var audibleSources: [SoundSource] {
        let anySolo = sources.contains { $0.solo }
        return sources.filter { !$0.muted && (!anySolo || $0.solo) }
    }

    /// A starter project: one room, one centred listener, and a source to hear immediately.
    public static func starter() -> SpatialProject {
        var p = SpatialProject(room: RoomPreset.builtIn[0].room, listener: Listener(position: Vec3(0, 1.6, 1)))
        var s = SoundSource(name: "Source 1", position: Vec3(1.5, 1.5, -1.5), content: .signal(.noiseBursts))
        s.directivity.dipoleWeight = 0.3
        p.sources.append(s)
        return p
    }
}
