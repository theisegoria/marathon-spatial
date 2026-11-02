import simd

/// World space: metres, right-handed, Y up. The plan view looks down -Y at the X-Z plane.
public typealias Vec3 = SIMD3<Double>

public struct Vec2: Codable, Hashable, Sendable {
    public var x: Double
    public var z: Double
    public init(_ x: Double, _ z: Double) { self.x = x; self.z = z }
    public static func + (a: Vec2, b: Vec2) -> Vec2 { Vec2(a.x + b.x, a.z + b.z) }
    public static func - (a: Vec2, b: Vec2) -> Vec2 { Vec2(a.x - b.x, a.z - b.z) }
    public static func * (a: Vec2, s: Double) -> Vec2 { Vec2(a.x * s, a.z * s) }
    public var length: Double { (x * x + z * z).squareRoot() }
    public func dot(_ o: Vec2) -> Double { x * o.x + z * o.z }
    public func cross(_ o: Vec2) -> Double { x * o.z - z * o.x }
}

public extension SIMD3 where Scalar == Double {
    var plan: Vec2 { Vec2(x, z) }
    init(plan p: Vec2, y: Double) { self.init(p.x, y, p.z) }
}

public enum Polygon2 {
    /// Shoelace area, positive when the vertices wind counter-clockwise in (x, z).
    public static func signedArea(_ p: [Vec2]) -> Double {
        guard p.count >= 3 else { return 0 }
        var a = 0.0
        for i in p.indices { a += p[i].cross(p[(i + 1) % p.count]) }
        return a / 2
    }

    public static func contains(_ p: [Vec2], _ q: Vec2) -> Bool {
        guard p.count >= 3 else { return false }
        var inside = false
        var j = p.count - 1
        for i in p.indices {
            if (p[i].z > q.z) != (p[j].z > q.z) {
                let t = (q.z - p[i].z) / (p[j].z - p[i].z)
                if q.x < p[i].x + t * (p[j].x - p[i].x) { inside.toggle() }
            }
            j = i
        }
        return inside
    }

    /// True when no two non-adjacent edges cross and there is no zero-length edge.
    public static func isSimple(_ p: [Vec2]) -> Bool {
        let n = p.count
        guard n >= 3 else { return false }
        for i in 0..<n where (p[(i + 1) % n] - p[i]).length < 1e-9 { return false }
        guard n > 3 else { return abs(signedArea(p)) > 1e-9 }
        for i in 0..<n {
            for j in (i + 1)..<n where j != i && (j + 1) % n != i && (i + 1) % n != j {
                if segmentsCross(p[i], p[(i + 1) % n], p[j], p[(j + 1) % n]) { return false }
            }
        }
        return abs(signedArea(p)) > 1e-9
    }

    static func segmentsCross(_ a: Vec2, _ b: Vec2, _ c: Vec2, _ d: Vec2) -> Bool {
        func orient(_ p: Vec2, _ q: Vec2, _ r: Vec2) -> Double { (q - p).cross(r - p) }
        let o1 = orient(a, b, c), o2 = orient(a, b, d), o3 = orient(c, d, a), o4 = orient(c, d, b)
        return (o1 > 0) != (o2 > 0) && (o3 > 0) != (o4 > 0)
    }

    /// Ear clipping for a simple polygon. Returns index triples, wound counter-clockwise in (x, z).
    public static func triangulate(_ input: [Vec2]) -> [(Int, Int, Int)] {
        guard input.count >= 3 else { return [] }
        var order = Array(input.indices)
        if signedArea(input) < 0 { order.reverse() }
        var out: [(Int, Int, Int)] = []
        var guardCount = order.count * order.count + 8
        while order.count > 3, guardCount > 0 {
            guardCount -= 1
            var clipped = false
            for i in order.indices {
                let ia = order[(i + order.count - 1) % order.count], ib = order[i], ic = order[(i + 1) % order.count]
                let a = input[ia], b = input[ib], c = input[ic]
                guard (b - a).cross(c - b) > 1e-12 else { continue }
                let blocked = order.contains { k in
                    guard k != ia, k != ib, k != ic else { return false }
                    let p = input[k]
                    return (b - a).cross(p - a) >= 0 && (c - b).cross(p - b) >= 0 && (a - c).cross(p - c) >= 0
                }
                if blocked { continue }
                out.append((ia, ib, ic))
                order.remove(at: i)
                clipped = true
                break
            }
            if !clipped { break }
        }
        if order.count == 3 { out.append((order[0], order[1], order[2])) }
        return out
    }
}

/// Closest point on segment ab to q, and the parameter along it.
public func closestOnSegment(_ a: Vec2, _ b: Vec2, _ q: Vec2) -> (point: Vec2, t: Double) {
    let d = b - a
    let l = d.dot(d)
    guard l > 0 else { return (a, 0) }
    let t = min(1, max(0, (q - a).dot(d) / l))
    return (a + d * t, t)
}
