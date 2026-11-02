import Foundation
import simd

/// What the user paints, and what a reflection bounces off.
public enum SurfaceRef: Hashable, Sendable {
    case wall(Int)
    case floor
    case ceiling
    case block(UUID)
}

/// A planar convex-or-concave polygon with a material.
public struct Surface: Sendable {
    public var ref: SurfaceRef
    public var material: MaterialID
    public var vertices: [Vec3]
    public var normal: Vec3
    public var area: Double
    /// The dominant axis dropped to project containment tests into 2D.
    fileprivate var dropAxis: Int

    init(ref: SurfaceRef, material: MaterialID, vertices: [Vec3], normal: Vec3, area: Double) {
        self.ref = ref; self.material = material; self.vertices = vertices; self.normal = normal; self.area = area
        let n = abs(normal)
        dropAxis = n.x >= n.y && n.x >= n.z ? 0 : (n.y >= n.z ? 1 : 2)
    }

    fileprivate func project(_ p: Vec3) -> Vec2 {
        switch dropAxis {
        case 0: Vec2(p.y, p.z)
        case 1: Vec2(p.x, p.z)
        default: Vec2(p.x, p.y)
        }
    }

    func contains(_ p: Vec3) -> Bool { Polygon2.contains(vertices.map(project), project(p)) }

    /// Mirror image of a point across the surface's plane.
    func mirror(_ p: Vec3) -> Vec3 { p - 2 * dot(p - vertices[0], normal) * normal }

    /// Where segment a-b crosses the surface, as a parameter along the segment, if it does.
    func intersect(_ a: Vec3, _ b: Vec3) -> Double? {
        let d = b - a
        let denom = dot(d, normal)
        guard abs(denom) > 1e-12 else { return nil }
        let t = dot(vertices[0] - a, normal) / denom
        guard t > 1e-6, t < 1 - 1e-6 else { return nil }
        return contains(a + d * t) ? t : nil
    }
}

public extension Room {
    /// Every sound-reflecting surface. Room normals point inward; block normals point outward.
    func surfaces() -> [Surface] {
        var out: [Surface] = []
        let ccw = Polygon2.signedArea(outline) > 0
        let n = outline.count
        for i in 0..<n {
            let a = outline[i], b = outline[(i + 1) % n]
            let e = b - a
            let len = e.length
            guard len > 0 else { continue }
            // For a polygon counter-clockwise in (x, z), the interior lies to the side (-e.z, e.x).
            let side = ccw ? 1.0 : -1.0
            let normal = Vec3(-e.z / len * side, 0, e.x / len * side)
            let verts = [Vec3(plan: a, y: 0), Vec3(plan: b, y: 0), Vec3(plan: b, y: height), Vec3(plan: a, y: height)]
            out.append(Surface(ref: .wall(i), material: walls[i], vertices: verts, normal: normal, area: len * height))
        }
        let area = abs(Polygon2.signedArea(outline))
        out.append(Surface(ref: .floor, material: floor, vertices: outline.map { Vec3(plan: $0, y: 0) }, normal: Vec3(0, 1, 0), area: area))
        out.append(Surface(ref: .ceiling, material: ceiling, vertices: outline.map { Vec3(plan: $0, y: height) }, normal: Vec3(0, -1, 0), area: area))
        for b in blocks {
            let lo = b.center - b.size / 2, hi = b.center + b.size / 2
            func quad(_ axis: Int, _ value: Double, _ sign: Double) {
                let u = (axis + 1) % 3, v = (axis + 2) % 3
                func at(_ s: Double, _ t: Double) -> Vec3 {
                    var p = Vec3.zero
                    p[axis] = value; p[u] = s; p[v] = t
                    return p
                }
                var normal = Vec3.zero
                normal[axis] = sign
                out.append(Surface(ref: .block(b.id), material: b.material,
                                   vertices: [at(lo[u], lo[v]), at(hi[u], lo[v]), at(hi[u], hi[v]), at(lo[u], hi[v])],
                                   normal: normal, area: (hi[u] - lo[u]) * (hi[v] - lo[v])))
            }
            for axis in 0..<3 {
                quad(axis, lo[axis], -1)
                quad(axis, hi[axis], 1)
            }
        }
        return out
    }

    /// Air volume in cubic metres, blocks excluded. Overlapping blocks are subtracted twice.
    public var volume: Double {
        let gross = abs(Polygon2.signedArea(outline)) * height
        return max(0, gross - blocks.reduce(0) { $0 + $1.size.x * $1.size.y * $1.size.z })
    }
}

public struct ReverbEstimate: Sendable, Equatable {
    public var bandsHz: [Double]
    public var sabineSeconds: [Double]
    public var eyringSeconds: [Double]
    public var midSeconds: Double { (eyringSeconds[3] + eyringSeconds[4]) / 2 }
}

public enum RoomAcoustics {
    /// Air attenuation in dB per km at 20 C and 50% humidity for `AcousticMaterial.bandsHz` (ISO 9613-1, rounded).
    static let airDBPerKm: [Double] = [0.12, 0.45, 1.3, 2.7, 4.7, 9.9, 29, 104]

    /// Statistical reverberation time in each band, Sabine and Eyring, with air absorption.
    public static func reverb(room: Room, materials: [AcousticMaterial]) -> ReverbEstimate {
        let surfaces = room.surfaces()
        let total = surfaces.reduce(0) { $0 + $1.area }
        let volume = room.volume
        var sabine: [Double] = [], eyring: [Double] = []
        for band in AcousticMaterial.bandsHz.indices {
            var a = 0.0
            for s in surfaces {
                let alpha = materials.first { $0.id == s.material }?.absorption[band] ?? 0.1
                a += s.area * alpha
            }
            let mean = total > 0 ? min(0.99, a / total) : 0
            let m = airDBPerKm[band] / 1000 / (10 / log(10.0))
            sabine.append(volume > 0 && a > 0 ? 0.161 * volume / (a + 4 * m * volume) : 0)
            let e = -total * log(1 - mean) + 4 * m * volume
            eyring.append(volume > 0 && e > 0 ? 0.161 * volume / e : 0)
        }
        return ReverbEstimate(bandsHz: AcousticMaterial.bandsHz, sabineSeconds: sabine, eyringSeconds: eyring)
    }
}

/// One specular path from a source to the listener.
public struct ReflectionPath: Sendable, Hashable, Identifiable {
    public var id: Int
    public var source: UUID
    /// Source, each bounce, listener.
    public var points: [Vec3]
    public var surfaces: [SurfaceRef]
    public var lengthMetres: Double
    /// Mid-band level relative to the same source at 1 m, including distance and surface loss.
    public var levelDB: Double
    /// Amplitude per `AcousticMaterial.bandsHz`, relative to the same source at 1 m.
    public var bandGains: [Double]
    public var order: Int { surfaces.count }
    public var delaySeconds: Double { lengthMetres / ImageSources.speedOfSound }
}

/// Early reflections by exact mirror images against this room, for drawing. What is heard comes from the library's simulator.
public enum ImageSources {
    public static let speedOfSound = 343.0

    public static func paths(room: Room, materials: [AcousticMaterial], source: SoundSource, sourcePosition: Vec3,
                             listener: Vec3, maxOrder: Int, limit: Int = 256) -> [ReflectionPath] {
        let surfaces = room.surfaces()
        let bandCount = AcousticMaterial.bandsHz.count
        let reflection = surfaces.map { s -> [Double] in
            let m = materials.first { $0.id == s.material }
            return (0..<bandCount).map { (1 - min(0.99, m?.absorption[$0] ?? 0.1)).squareRoot() }
        }
        var out: [ReflectionPath] = []
        var nextID = 0

        func occluded(_ a: Vec3, _ b: Vec3, ignoring: [Int]) -> Bool {
            for i in surfaces.indices where !ignoring.contains(i) {
                if surfaces[i].intersect(a, b) != nil { return true }
            }
            return false
        }

        func walk(sequence: [Int], image: Vec3) {
            if sequence.count <= 1, Task.isCancelled { return }
            if !sequence.isEmpty, let path = resolve(sequence: sequence) { out.append(path) }
            guard sequence.count < maxOrder else { return }
            for i in surfaces.indices where i != sequence.last {
                // An image behind the surface cannot reflect off its front.
                guard dot(image - surfaces[i].vertices[0], surfaces[i].normal) > 1e-9 else { continue }
                walk(sequence: sequence + [i], image: surfaces[i].mirror(image))
            }
        }

        /// Walk back from the listener through the mirrored images to find the bounce points, then check each leg.
        func resolve(sequence: [Int]) -> ReflectionPath? {
            var images = [sourcePosition]
            for i in sequence { images.append(surfaces[i].mirror(images.last!)) }
            var bounces: [Vec3] = []
            var target = listener
            for k in stride(from: sequence.count - 1, through: 0, by: -1) {
                let s = surfaces[sequence[k]]
                guard let t = s.intersect(target, images[k + 1]) else { return nil }
                let hit = target + (images[k + 1] - target) * t
                bounces.insert(hit, at: 0)
                target = hit
            }
            let points = [sourcePosition] + bounces + [listener]
            var length = 0.0
            var gains = [Double](repeating: 1, count: bandCount)
            for leg in 0..<(points.count - 1) {
                let ignore = [leg > 0 ? sequence[leg - 1] : -1, leg < sequence.count ? sequence[leg] : -1]
                if occluded(points[leg], points[leg + 1], ignoring: ignore) { return nil }
                length += simd_distance(points[leg], points[leg + 1])
            }
            let spreading = 1 / max(length, 1e-3)
            for b in 0..<bandCount {
                for i in sequence { gains[b] *= reflection[i][b] }
                gains[b] *= spreading
            }
            let level = 20 * log10(max((gains[3] + gains[4] + gains[5]) / 3, 1e-9))
            defer { nextID += 1 }
            return ReflectionPath(id: nextID, source: source.id, points: points, surfaces: sequence.map { surfaces[$0].ref },
                                  lengthMetres: length, levelDB: level, bandGains: gains)
        }

        walk(sequence: [], image: sourcePosition)
        out.sort { $0.levelDB > $1.levelDB }
        return Array(out.prefix(limit))
    }

    /// True when a straight line from source to listener is blocked.
    public static func directIsOccluded(room: Room, from a: Vec3, to b: Vec3) -> Bool {
        room.surfaces().contains { $0.intersect(a, b) != nil }
    }
}
