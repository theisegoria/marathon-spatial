import Foundation

public enum WAVFormat: String, Codable, Hashable, Sendable, CaseIterable, Identifiable {
    case pcm16, pcm24, float32
    public var id: String { rawValue }
    public var title: String {
        switch self {
        case .pcm16: "16-bit PCM"
        case .pcm24: "24-bit PCM"
        case .float32: "32-bit float"
        }
    }
    var bytesPerSample: Int { self == .pcm16 ? 2 : (self == .pcm24 ? 3 : 4) }
}

public struct WAVSummary: Sendable, Equatable {
    public var frames: Int64
    public var peak: Float
    /// Samples that exceeded full scale and were clipped; zero for float files, which keep them.
    public var clippedSamples: Int64
    public var bytes: Int64
}

public enum WAVError: Error, LocalizedError {
    case tooLarge
    case channelCount
    public var errorDescription: String? {
        switch self {
        case .tooLarge: "The render would exceed the 4 GB WAV limit. Shorten it or lower the sample rate."
        case .channelCount: "Interleaved data does not match the channel count."
        }
    }
}

/// A streaming RIFF/WAVE writer. Fixed-point formats get TPDF dither from a seeded generator, so output is reproducible.
public final class WAVWriter {
    private let handle: FileHandle
    private let format: WAVFormat
    private let channels: Int
    private var frames: Int64 = 0
    private var peak: Float = 0
    private var clipped: Int64 = 0
    private var dither = Xorshift(seed: 0x4D53_5041)
    private var finished = false

    public init(url: URL, format: WAVFormat, sampleRate: Int, channels: Int) throws {
        guard channels > 0 else { throw WAVError.channelCount }
        guard FileManager.default.createFile(atPath: url.path, contents: nil) else { throw CocoaError(.fileWriteUnknown) }
        handle = try FileHandle(forWritingTo: url)
        self.format = format
        self.channels = channels
        try handle.write(contentsOf: Self.header(format: format, sampleRate: sampleRate, channels: channels, dataBytes: 0))
    }

    static func header(format: WAVFormat, sampleRate: Int, channels: Int, dataBytes: UInt32) -> Data {
        var d = Data()
        func put32(_ v: UInt32) { withUnsafeBytes(of: v.littleEndian) { d.append(contentsOf: $0) } }
        func put16(_ v: UInt16) { withUnsafeBytes(of: v.littleEndian) { d.append(contentsOf: $0) } }
        let block = channels * format.bytesPerSample
        d.append(contentsOf: Array("RIFF".utf8)); put32(36 &+ dataBytes)
        d.append(contentsOf: Array("WAVE".utf8))
        d.append(contentsOf: Array("fmt ".utf8)); put32(16)
        put16(format == .float32 ? 3 : 1)
        put16(UInt16(channels)); put32(UInt32(sampleRate)); put32(UInt32(sampleRate * block))
        put16(UInt16(block)); put16(UInt16(format.bytesPerSample * 8))
        d.append(contentsOf: Array("data".utf8)); put32(dataBytes)
        return d
    }

    /// Append interleaved samples.
    public func append(_ interleaved: UnsafeBufferPointer<Float>) throws {
        guard interleaved.count % channels == 0 else { throw WAVError.channelCount }
        let added = Int64(interleaved.count / channels)
        let bytes = (frames + added) * Int64(channels) * Int64(format.bytesPerSample)
        guard bytes + 36 <= Int64(UInt32.max) else { throw WAVError.tooLarge }

        var out = Data(count: interleaved.count * format.bytesPerSample)
        let scale: Float = format == .pcm16 ? 32767 : 8_388_607
        out.withUnsafeMutableBytes { (raw: UnsafeMutableRawBufferPointer) in
            for (i, sample) in interleaved.enumerated() {
                let finite = sample.isFinite ? sample : 0
                peak = max(peak, abs(finite))
                switch format {
                case .float32:
                    raw.storeBytes(of: finite.bitPattern.littleEndian, toByteOffset: i * 4, as: UInt32.self)
                case .pcm16, .pcm24:
                    if abs(finite) > 1 { clipped += 1 }
                    let lsb = dither.nextFloat() - dither.nextFloat()
                    let v = max(-scale - 1, min(scale, (max(-1, min(1, finite)) * scale + lsb).rounded()))
                    if format == .pcm16 {
                        raw.storeBytes(of: Int16(v).littleEndian, toByteOffset: i * 2, as: Int16.self)
                    } else {
                        let w = Int32(v)
                        raw[i * 3] = UInt8(truncatingIfNeeded: w)
                        raw[i * 3 + 1] = UInt8(truncatingIfNeeded: w >> 8)
                        raw[i * 3 + 2] = UInt8(truncatingIfNeeded: w >> 16)
                    }
                }
            }
        }
        try handle.write(contentsOf: out)
        frames += added
    }

    /// Patch the sizes into the header and close the file.
    public func finish() throws -> WAVSummary {
        let bytes = frames * Int64(channels) * Int64(format.bytesPerSample)
        if !finished {
            finished = true
            try handle.seek(toOffset: 4)
            try handle.write(contentsOf: withUnsafeBytes(of: UInt32(36 + bytes).littleEndian) { Data($0) })
            try handle.seek(toOffset: 40)
            try handle.write(contentsOf: withUnsafeBytes(of: UInt32(bytes).littleEndian) { Data($0) })
            try handle.close()
        }
        return WAVSummary(frames: frames, peak: peak, clippedSamples: clipped, bytes: bytes + 44)
    }
}
