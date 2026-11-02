import AVFoundation
import CryptoKit
import Foundation

public enum ClipError: Error, LocalizedError {
    case unreadable(String)
    case tooLong(String)
    public var errorDescription: String? {
        switch self {
        case .unreadable(let name): "\(name) could not be read as audio."
        case .tooLong(let name): "\(name) is too long to load into memory."
        }
    }
}

/// Decoded PCM, one array per channel. Reading a file through AVAudioFile touches no output device.
public struct DecodedClip: Sendable {
    public var channels: [[Float]]
    public var sampleRate: Double
    public var frames: Int { channels.first?.count ?? 0 }

    /// Mono, as the sum of the channels at equal power. A point source plays this.
    public var monoMix: [Float] {
        guard channels.count > 1 else { return channels.first ?? [] }
        let scale = 1 / Float(channels.count).squareRoot()
        var out = [Float](repeating: 0, count: frames)
        for channel in channels {
            for i in 0..<out.count { out[i] += channel[i] * scale }
        }
        return out
    }

    public var interleaved: [Float] {
        var out = [Float](repeating: 0, count: frames * channels.count)
        for (c, channel) in channels.enumerated() {
            for i in 0..<frames { out[i * channels.count + c] = channel[i] }
        }
        return out
    }
    public init(channels: [[Float]], sampleRate: Double) { self.channels = channels; self.sampleRate = sampleRate }
}

/// Content-addressed clip files shared by all open documents. A document package holds its own copies; opening one fills
/// this cache, which keeps the document value small enough to diff and undo.
public actor ClipStore {
    public static let maximumFrames: Int64 = 48_000 * 60 * 60

    public static let shared = ClipStore(directory: FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask)[0]
        .appendingPathComponent("Marathon Spatial Studio/Clips"))

    public nonisolated let directory: URL

    public init(directory: URL) { self.directory = directory }

    public nonisolated func url(for clip: ClipReference) -> URL { directory.appendingPathComponent(clip.storedName) }

    public nonisolated func url(storedName: String) -> URL { directory.appendingPathComponent(storedName) }

    private nonisolated func prepare() throws { try FileManager.default.createDirectory(at: directory, withIntermediateDirectories: true) }

    public nonisolated static func sha256(of url: URL) throws -> String {
        let handle = try FileHandle(forReadingFrom: url)
        defer { try? handle.close() }
        var hash = SHA256()
        while let chunk = try handle.read(upToCount: 1 << 20), !chunk.isEmpty { hash.update(data: chunk) }
        return hash.finalize().map { String(format: "%02x", $0) }.joined()
    }

    /// Copy a user's audio file into the cache and describe it.
    public func importFile(_ source: URL) throws -> ClipReference {
        let scoped = source.startAccessingSecurityScopedResource()
        defer { if scoped { source.stopAccessingSecurityScopedResource() } }
        let name = source.deletingPathExtension().lastPathComponent
        let file: AVAudioFile
        do { file = try AVAudioFile(forReading: source) } catch { throw ClipError.unreadable(source.lastPathComponent) }
        guard file.length <= Self.maximumFrames else { throw ClipError.tooLong(source.lastPathComponent) }
        let sha = try Self.sha256(of: source)
        let ref = ClipReference(name: name, sha256: sha, fileExtension: source.pathExtension.lowercased(), frames: file.length,
                                sampleRate: file.fileFormat.sampleRate, channels: Int(file.fileFormat.channelCount))
        try prepare()
        let dest = url(for: ref)
        if !FileManager.default.fileExists(atPath: dest.path) { try FileManager.default.copyItem(at: source, to: dest) }
        return ref
    }

    /// Store bytes read from a document package, if they are not cached already.
    public nonisolated func adopt(storedName: String, data: Data) throws {
        try prepare()
        let dest = url(storedName: storedName)
        if !FileManager.default.fileExists(atPath: dest.path) { try data.write(to: dest, options: .atomic) }
    }

    public func decode(_ clip: ClipReference) throws -> DecodedClip {
        let file: AVAudioFile
        do {
            file = try AVAudioFile(forReading: url(for: clip), commonFormat: .pcmFormatFloat32, interleaved: false)
        } catch { throw ClipError.unreadable(clip.name) }
        let capacity = AVAudioFrameCount(min(file.length, Self.maximumFrames))
        guard let buffer = AVAudioPCMBuffer(pcmFormat: file.processingFormat, frameCapacity: capacity), capacity > 0 else {
            throw ClipError.unreadable(clip.name)
        }
        try file.read(into: buffer)
        let frames = Int(buffer.frameLength)
        guard let data = buffer.floatChannelData else { throw ClipError.unreadable(clip.name) }
        let channels = (0..<Int(buffer.format.channelCount)).map { Array(UnsafeBufferPointer(start: data[$0], count: frames)) }
        return DecodedClip(channels: channels, sampleRate: file.processingFormat.sampleRate)
    }

    /// Peak per bucket across all channels, for a waveform drawing.
    public func overview(_ clip: ClipReference, buckets: Int) throws -> [Float] {
        let decoded = try decode(clip)
        return Self.overview(decoded, buckets: buckets)
    }

    public nonisolated static func overview(_ clip: DecodedClip, buckets: Int) -> [Float] {
        guard buckets > 0, clip.frames > 0 else { return [] }
        var out = [Float](repeating: 0, count: buckets)
        for channel in clip.channels {
            for i in 0..<clip.frames {
                let b = i * buckets / clip.frames
                out[b] = max(out[b], abs(channel[i]))
            }
        }
        return out
    }
}
