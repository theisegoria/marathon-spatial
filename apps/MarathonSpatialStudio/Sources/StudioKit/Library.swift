import Foundation

/// Editable starting points. Absorption is typical published octave-band data for each construction, rounded. Published
/// tables start at 125 Hz, so the 62.5 Hz band repeats that value. A project keeps its own copy, so a user can enter measured values.
public enum MaterialLibrary {
    public static let acoustic: [AcousticMaterial] = [
        AcousticMaterial(id: "concrete", name: "Painted concrete", absorption: [0.01, 0.01, 0.01, 0.02, 0.02, 0.02, 0.03, 0.04], scattering: 0.10, transmissionDB: -50, swatch: 0),
        AcousticMaterial(id: "brick", name: "Bare brick", absorption: [0.03, 0.03, 0.03, 0.03, 0.04, 0.05, 0.07, 0.08], scattering: 0.25, transmissionDB: -45, swatch: 1),
        AcousticMaterial(id: "plaster", name: "Plaster on masonry", absorption: [0.02, 0.02, 0.03, 0.04, 0.05, 0.05, 0.05, 0.06], scattering: 0.10, transmissionDB: -45, swatch: 2),
        AcousticMaterial(id: "gypsum", name: "Gypsum board", absorption: [0.29, 0.29, 0.10, 0.05, 0.04, 0.07, 0.09, 0.09], scattering: 0.10, transmissionDB: -30, swatch: 3),
        AcousticMaterial(id: "glass", name: "Window glass", absorption: [0.18, 0.18, 0.06, 0.04, 0.03, 0.02, 0.02, 0.02], scattering: 0.05, transmissionDB: -28, swatch: 4),
        AcousticMaterial(id: "wood", name: "Wood panelling", absorption: [0.28, 0.28, 0.22, 0.17, 0.09, 0.10, 0.11, 0.11], scattering: 0.15, transmissionDB: -25, swatch: 5),
        AcousticMaterial(id: "floorboards", name: "Wood floor", absorption: [0.15, 0.15, 0.11, 0.10, 0.07, 0.06, 0.07, 0.07], scattering: 0.10, transmissionDB: -30, swatch: 6),
        AcousticMaterial(id: "carpet", name: "Heavy carpet", absorption: [0.08, 0.08, 0.24, 0.57, 0.69, 0.71, 0.73, 0.75], scattering: 0.20, transmissionDB: -40, swatch: 7),
        AcousticMaterial(id: "curtain", name: "Heavy curtain", absorption: [0.14, 0.14, 0.35, 0.55, 0.72, 0.70, 0.65, 0.65], scattering: 0.40, transmissionDB: -10, swatch: 8),
        AcousticMaterial(id: "acoustictile", name: "Acoustic ceiling tile", absorption: [0.50, 0.50, 0.70, 0.60, 0.70, 0.70, 0.50, 0.50], scattering: 0.15, transmissionDB: -20, swatch: 9),
        AcousticMaterial(id: "steel", name: "Steel plate", absorption: [0.04, 0.04, 0.04, 0.03, 0.03, 0.02, 0.02, 0.02], scattering: 0.10, transmissionDB: -35, swatch: 10),
        AcousticMaterial(id: "marble", name: "Polished stone", absorption: [0.01, 0.01, 0.01, 0.01, 0.01, 0.02, 0.02, 0.02], scattering: 0.05, transmissionDB: -55, swatch: 11),
        AcousticMaterial(id: "water", name: "Water surface", absorption: [0.01, 0.01, 0.01, 0.01, 0.015, 0.02, 0.025, 0.03], scattering: 0.05, transmissionDB: -60, swatch: 12),
    ]

    /// The library's presets for these names: handbook values, starting points for tuning by ear.
    public static let physical: [PhysicalMaterial] = [
        PhysicalMaterial(id: "steel", name: "Steel", densityKgM3: 7850, youngsModulusPa: 200e9, poissonRatio: 0.29, lossFactor: 6e-4),
        PhysicalMaterial(id: "sheet-metal", name: "Sheet metal", densityKgM3: 7850, youngsModulusPa: 200e9, poissonRatio: 0.29, lossFactor: 1.5e-3),
        PhysicalMaterial(id: "aluminium", name: "Aluminium", densityKgM3: 2700, youngsModulusPa: 69e9, poissonRatio: 0.33, lossFactor: 4e-4),
        PhysicalMaterial(id: "cast-iron", name: "Cast iron", densityKgM3: 7200, youngsModulusPa: 120e9, poissonRatio: 0.26, lossFactor: 2e-3),
        PhysicalMaterial(id: "brass", name: "Brass", densityKgM3: 8500, youngsModulusPa: 100e9, poissonRatio: 0.34, lossFactor: 1e-3),
        PhysicalMaterial(id: "copper", name: "Copper", densityKgM3: 8960, youngsModulusPa: 120e9, poissonRatio: 0.34, lossFactor: 2e-3),
        PhysicalMaterial(id: "glass", name: "Glass", densityKgM3: 2500, youngsModulusPa: 70e9, poissonRatio: 0.22, lossFactor: 1e-3),
        PhysicalMaterial(id: "ceramic", name: "Ceramic", densityKgM3: 2400, youngsModulusPa: 70e9, poissonRatio: 0.20, lossFactor: 1.5e-3),
        PhysicalMaterial(id: "concrete", name: "Concrete", densityKgM3: 2400, youngsModulusPa: 30e9, poissonRatio: 0.20, lossFactor: 2e-2),
        PhysicalMaterial(id: "granite", name: "Granite", densityKgM3: 2700, youngsModulusPa: 50e9, poissonRatio: 0.25, lossFactor: 3e-3),
        PhysicalMaterial(id: "hardwood", name: "Hardwood", densityKgM3: 720, youngsModulusPa: 12e9, poissonRatio: 0.35, lossFactor: 1e-2),
        PhysicalMaterial(id: "plywood", name: "Plywood", densityKgM3: 600, youngsModulusPa: 8e9, poissonRatio: 0.30, lossFactor: 1.5e-2),
        PhysicalMaterial(id: "acrylic", name: "Acrylic", densityKgM3: 1190, youngsModulusPa: 3.2e9, poissonRatio: 0.37, lossFactor: 3e-2),
        PhysicalMaterial(id: "rubber", name: "Rubber", densityKgM3: 1100, youngsModulusPa: 5e6, poissonRatio: 0.49, lossFactor: 0.2),
    ]
}

public struct RoomPreset: Codable, Hashable, Sendable, Identifiable {
    public var id: String
    public var name: String
    public var room: Room
    public var isUser = false

    public init(id: String, name: String, room: Room, isUser: Bool = false) {
        self.id = id; self.name = name; self.room = room; self.isUser = isUser
    }

    public static let builtIn: [RoomPreset] = [
        RoomPreset(id: "living", name: "Living room", room: .shoebox(width: 5, depth: 6.5, height: 2.6, walls: "plaster", floor: "carpet", ceiling: "gypsum")),
        RoomPreset(id: "studio", name: "Treated studio", room: .shoebox(width: 6, depth: 8, height: 3, walls: "curtain", floor: "carpet", ceiling: "acoustictile")),
        RoomPreset(id: "corridor", name: "Concrete corridor", room: .shoebox(width: 2.4, depth: 24, height: 2.7, walls: "concrete", floor: "concrete", ceiling: "concrete")),
        RoomPreset(id: "hall", name: "Wooden hall", room: .shoebox(width: 20, depth: 34, height: 12, walls: "wood", floor: "floorboards", ceiling: "plaster")),
        RoomPreset(id: "cathedral", name: "Stone nave", room: .shoebox(width: 18, depth: 60, height: 26, walls: "marble", floor: "marble", ceiling: "plaster")),
        RoomPreset(id: "hangar", name: "Steel hangar", room: .shoebox(width: 40, depth: 60, height: 16, walls: "steel", floor: "concrete", ceiling: "steel")),
        RoomPreset(id: "bathroom", name: "Tiled bathroom", room: .shoebox(width: 2.2, depth: 3, height: 2.4, walls: "glass", floor: "marble", ceiling: "plaster")),
        RoomPreset(id: "cistern", name: "Flooded cistern", room: .shoebox(width: 14, depth: 22, height: 7, walls: "concrete", floor: "water", ceiling: "concrete")),
    ]
}

/// A saved view of how to render: the part of a project worth reusing across rooms.
public struct RenderPreset: Codable, Hashable, Sendable, Identifiable {
    public var id: String
    public var name: String
    public var settings: RenderSettings
    public var isUser = false

    public init(id: String, name: String, settings: RenderSettings, isUser: Bool = false) {
        self.id = id; self.name = name; self.settings = settings; self.isUser = isUser
    }

    public static let builtIn: [RenderPreset] = {
        var dry = RenderSettings(); dry.earlyReflections = false; dry.lateReverb = false
        var early = RenderSettings(); early.lateReverb = false
        var full = RenderSettings()
        var speakers = RenderSettings(); speakers.output = .stereo
        var laptop = RenderSettings(); laptop.output = .speakers3D
        var desk = RenderSettings(); desk.output = .speakers3D
        desk.crosstalk = { var c = CrosstalkSettings(); c.spanDegrees = 30; c.distanceMetres = 0.8; return c }()
        return [
            RenderPreset(id: "full", name: "Full binaural", settings: full),
            RenderPreset(id: "early", name: "Direct and early reflections", settings: early),
            RenderPreset(id: "dry", name: "Direct sound only", settings: dry),
            RenderPreset(id: "speakers", name: "Stereo speakers", settings: speakers),
            RenderPreset(id: "laptop3d", name: "Laptop speakers, 3D", settings: laptop),
            RenderPreset(id: "desk3d", name: "Desk speakers, 3D", settings: desk),
        ]
    }()
}

/// User presets, as JSON in Application Support.
public actor PresetStore {
    public struct Contents: Codable, Sendable {
        public var rooms: [RoomPreset] = []
        public var renders: [RenderPreset] = []
        public init() {}
    }

    private let url: URL

    public init(directory: URL) { url = directory.appendingPathComponent("Presets.json") }

    public func load() -> Contents {
        guard let data = try? Data(contentsOf: url), let c = try? JSONDecoder().decode(Contents.self, from: data) else { return Contents() }
        return c
    }

    public func save(_ c: Contents) throws {
        try FileManager.default.createDirectory(at: url.deletingLastPathComponent(), withIntermediateDirectories: true)
        let e = JSONEncoder()
        e.outputFormatting = [.prettyPrinted, .sortedKeys]
        try e.encode(c).write(to: url, options: .atomic)
    }

    public func add(room: RoomPreset) throws {
        var c = load(); c.rooms.removeAll { $0.id == room.id }; c.rooms.append(room); try save(c)
    }

    public func add(render: RenderPreset) throws {
        var c = load(); c.renders.removeAll { $0.id == render.id }; c.renders.append(render); try save(c)
    }

    public func remove(id: String) throws {
        var c = load(); c.rooms.removeAll { $0.id == id }; c.renders.removeAll { $0.id == id }; try save(c)
    }
}
