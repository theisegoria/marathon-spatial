import Foundation

/// A `.mspatial` project is a package: `project.json` plus `Audio/<sha256>.<ext>`. Everything that makes it sound travels
/// with it, except the HRTF sets, which belong to the library.
public enum ProjectPackage {
    public static let manifestName = "project.json"
    public static let audioDirectory = "Audio"

    public enum Failure: Error, LocalizedError {
        case missingManifest
        case newerVersion(Int)
        public var errorDescription: String? {
            switch self {
            case .missingManifest: "The project has no project.json."
            case .newerVersion(let v): "This project was saved by a newer version (format \(v))."
            }
        }
    }

    public static func encode(_ project: SpatialProject) throws -> Data {
        let e = JSONEncoder()
        e.outputFormatting = [.prettyPrinted, .sortedKeys]
        return try e.encode(project)
    }

    public static func decode(_ data: Data) throws -> SpatialProject {
        struct Header: Decodable { var version: Int }
        let header = try JSONDecoder().decode(Header.self, from: data)
        guard header.version <= SpatialProject.formatVersion else { throw Failure.newerVersion(header.version) }
        return try JSONDecoder().decode(SpatialProject.self, from: data)
    }

    public static func read(_ wrapper: FileWrapper, store: ClipStore) throws -> SpatialProject {
        guard let manifest = wrapper.fileWrappers?[manifestName]?.regularFileContents else { throw Failure.missingManifest }
        let project = try decode(manifest)
        if let audio = wrapper.fileWrappers?[audioDirectory]?.fileWrappers {
            for (name, file) in audio {
                if let data = file.regularFileContents { try store.adopt(storedName: name, data: data) }
            }
        }
        return project
    }

    /// Build the package. Clips are referenced from the cache rather than held in memory.
    public static func write(_ project: SpatialProject, store: ClipStore) throws -> FileWrapper {
        var audio: [String: FileWrapper] = [:]
        for clip in project.clips where audio[clip.storedName] == nil {
            let file = try FileWrapper(url: store.url(for: clip), options: [])
            file.preferredFilename = clip.storedName
            audio[clip.storedName] = file
        }
        let manifest = FileWrapper(regularFileWithContents: try encode(project))
        manifest.preferredFilename = manifestName
        let audioWrapper = FileWrapper(directoryWithFileWrappers: audio)
        audioWrapper.preferredFilename = audioDirectory
        return FileWrapper(directoryWithFileWrappers: [manifestName: manifest, audioDirectory: audioWrapper])
    }
}
