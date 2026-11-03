// swift-tools-version: 6.2
// Marathon Spatial Studio: a native macOS 27 tool for spatial and physics-based audio on any headphones.
// StudioKit holds the document, geometry, acoustics and the engine boundary, with no SwiftUI. The app target holds the
// interface. Every sound the app makes goes through the Marathon Spatial library.
//
// The library is built by its own CMake project (see docs/BUILDING.md). This manifest looks for its static libraries in
// MARATHON_SPATIAL_LIB_DIR, or else through ../../build. With all of them found, the app links them and LibraryEngine replaces
// UnavailableEngine. Without them it still builds, and says in the window that it cannot make sound. Scripts/bundle.sh sets
// MARATHON_SPATIAL_REQUIRE_LINKED, so a bundle can never be the silent variant by accident.
import Foundation
import PackageDescription

let here = URL(fileURLWithPath: #filePath).deletingLastPathComponent().path
let libraryRoot = URL(fileURLWithPath: "\(here)/../..").standardized.path
// Everything but device output. StudioRender and its tests link only these, so they cannot open a device.
let renderLibraries = ["marathon_spatial_core", "marathon_spatial_mysofa", "marathon_spatial_physics", "marathon_spatial_acoustics"]
let spatialLibraries = renderLibraries + ["marathon_spatial_output"]
let environment = ProcessInfo.processInfo.environment
let buildDirectory = environment["MARATHON_SPATIAL_LIB_DIR"] ?? "\(libraryRoot)/build"
let searchDirectories = [buildDirectory, "\(buildDirectory)/lib", "\(buildDirectory)/src/physics"]
func directory(of name: String) -> String? { searchDirectories.first { FileManager.default.fileExists(atPath: "\($0)/lib\(name).a") } }
let missing = spatialLibraries.filter { directory(of: $0) == nil }
let linked = missing.isEmpty
if !linked, environment["MARATHON_SPATIAL_REQUIRE_LINKED"] != nil {
    fatalError("Not found in \(searchDirectories): \(missing.map { "lib\($0).a" }.joined(separator: ", ")). Build the library first, see docs/BUILDING.md.")
}

func linkFlags(_ names: [String]) -> [String] {
    Set(names.compactMap(directory(of:))).sorted().map { "-L\($0)" } + names.map { "-l\($0)" } + ["-lc++"]
}
let frameworks: [LinkerSetting] = [
    .linkedFramework("Accelerate"), .linkedFramework("AudioToolbox"), .linkedFramework("CoreFoundation"), .linkedFramework("Foundation"),
    .linkedFramework("Metal"), .linkedLibrary("z"),
]

var targets: [Target] = [
    .target(name: "StudioKit", swiftSettings: [.swiftLanguageMode(.v6)]),
    .testTarget(name: "StudioKitTests", dependencies: ["StudioKit"], swiftSettings: [.swiftLanguageMode(.v6)]),
]

if linked {
    targets.append(.target(name: "CStudioBridge", cxxSettings: [.unsafeFlags(["-std=c++23"])]))
    targets.append(.target(name: "StudioRender", dependencies: ["StudioKit", "CStudioBridge"], swiftSettings: [.swiftLanguageMode(.v6)]))
    targets.append(.executableTarget(
        name: "MarathonSpatialStudio", dependencies: ["StudioKit", "StudioRender", "CStudioBridge"],
        swiftSettings: [.swiftLanguageMode(.v6), .define("MARATHON_SPATIAL_LINKED")],
        linkerSettings: [.unsafeFlags(linkFlags(spatialLibraries)), .linkedFramework("CoreAudio"), .linkedFramework("AVFoundation"), .linkedFramework("Vision"),
                         .linkedFramework("CoreMedia"), .linkedFramework("CoreVideo"), .linkedFramework("QuartzCore")] + frameworks))
    targets.append(.testTarget(
        name: "StudioRenderTests", dependencies: ["StudioKit", "StudioRender"], swiftSettings: [.swiftLanguageMode(.v6)],
        linkerSettings: [.unsafeFlags(linkFlags(renderLibraries))] + frameworks))
} else {
    targets.append(.executableTarget(name: "MarathonSpatialStudio", dependencies: ["StudioKit"],
                                     swiftSettings: [.swiftLanguageMode(.v6)]))
}

let package = Package(
    name: "MarathonSpatialStudio",
    platforms: [.macOS("27.0")],
    products: [
        .library(name: "StudioKit", targets: ["StudioKit"]),
        .executable(name: "MarathonSpatialStudio", targets: ["MarathonSpatialStudio"]),
    ],
    targets: targets
)
