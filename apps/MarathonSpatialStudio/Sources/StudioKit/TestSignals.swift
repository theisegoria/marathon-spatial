import Foundation

/// A small seeded generator, for reproducible noise and dither.
public struct Xorshift: Sendable {
    private var state: UInt64
    public init(seed: UInt64) { state = seed == 0 ? 0x9E37_79B9_7F4A_7C15 : seed }
    public mutating func next() -> UInt64 {
        state ^= state << 13
        state ^= state >> 7
        state ^= state << 17
        return state
    }
    /// Uniform in 0..<1.
    public mutating func nextFloat() -> Float { Float(next() >> 40) / Float(1 << 24) }
}

public extension TestSignal {
    /// Mono samples at `sampleRate`, peak about -6 dBFS. Same input, same output.
    func generate(sampleRate: Double, seconds: Double, seed: UInt64 = 1) -> [Float] {
        let n = max(0, Int(sampleRate * seconds))
        var out = [Float](repeating: 0, count: n)
        var rng = Xorshift(seed: seed)
        func white() -> Float { rng.nextFloat() * 2 - 1 }
        switch self {
        case .whiteNoise:
            for i in 0..<n { out[i] = 0.5 * white() }
        case .pinkNoise:
            // Paul Kellet's filter.
            var b0: Float = 0, b1: Float = 0, b2: Float = 0
            for i in 0..<n {
                let w = white()
                b0 = 0.99765 * b0 + w * 0.0990460
                b1 = 0.96300 * b1 + w * 0.2965164
                b2 = 0.57000 * b2 + w * 1.0526913
                out[i] = 0.125 * (b0 + b1 + b2 + w * 0.1848)
            }
        case .noiseBursts:
            // 80 ms of noise every 500 ms, with 2 ms ramps.
            let period = max(1, Int(sampleRate * 0.5)), burst = Int(sampleRate * 0.08), ramp = max(1, Int(sampleRate * 0.002))
            for i in 0..<n {
                let p = i % period
                guard p < burst else { continue }
                let env = Float(min(1, min(Double(p) / Double(ramp), Double(burst - 1 - p) / Double(ramp))))
                out[i] = 0.5 * env * white()
            }
        case .clickTrain:
            let period = max(1, Int(sampleRate * 0.25))
            var i = 0
            while i < n { out[i] = 0.5; i += period }
        case .sineSweep:
            // Exponential sweep, 20 Hz to 20 kHz, over the whole signal.
            let f0 = 20.0, f1 = min(20000.0, sampleRate * 0.45), k = log(f1 / f0)
            let duration = max(seconds, 1e-6)
            for i in 0..<n {
                let t = Double(i) / sampleRate
                out[i] = Float(0.5 * sin(2 * .pi * f0 * duration / k * (exp(t / duration * k) - 1)))
            }
        case .tone440:
            for i in 0..<n { out[i] = Float(0.5 * sin(2 * .pi * 440 * Double(i) / sampleRate)) }
        }
        return out
    }
}
