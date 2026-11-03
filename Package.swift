// swift-tools-version: 6.2
// Marathon Spatial: spatial and physics-based audio for any headphones. macOS 27+, Apple silicon M3+.
// CMarathonSpatial is the engine with no device code, CMarathonSpatialOutput the Core Audio output, and MarathonSpatial
// the Swift wrapper over both. MarathonSpatialHRTFData carries the vendored HRTF sets as a resource bundle.
import PackageDescription

// The engine's C targets sit at the package root, beside everything that is not theirs.
let notEngine = ["LICENSE", "LICENSES", "NOTICE", "README.md", "REUSE.toml", "SECURITY.md", "CONTRIBUTING.md", "THIRD_PARTY_NOTICES.md",
                 "CMakeLists.txt", "apps", "cmake", "data", "demo", "docs", "examples", "shaders", "swift", "tests", "third_party", "tools"]

let package = Package(
    name: "MarathonSpatial",
    platforms: [.macOS("27.0")],
    products: [
        .library(name: "MarathonSpatial", targets: ["MarathonSpatial"]),
        .library(name: "MarathonSpatialHRTFData", targets: ["MarathonSpatialHRTFData"]),
    ],
    targets: [
        .target(
            name: "CMarathonSpatialMysofa",
            path: "third_party/libmysofa",
            sources: ["src/hdf", "src/hrtf"],
            publicHeadersPath: "src/hrtf",
            cSettings: [.unsafeFlags(["-w"])],
            linkerSettings: [.linkedLibrary("z")]),
        .target(
            name: "CMarathonSpatial",
            dependencies: ["CMarathonSpatialMysofa"],
            path: ".",
            exclude: notEngine,
            sources: ["src/core/Api.cpp", "src/core/Ambisonics.cpp", "src/core/Dsp.cpp", "src/core/Engine.cpp", "src/core/Speakers.cpp", "src/core/DeviceChannels.cpp", "src/core/Crosstalk.cpp",
                      "src/hrtf/HrtfBuild.cpp", "src/hrtf/HrtfFile.cpp", "src/hrtf/SphereHead.cpp"],
            publicHeadersPath: "include/CMarathonSpatial",
            cxxSettings: [.headerSearchPath("src"), .headerSearchPath("include"), .unsafeFlags(["-O3", "-fno-math-errno"], .when(configuration: .release))],
            linkerSettings: [.linkedFramework("Accelerate"), .linkedFramework("AudioToolbox"), .linkedFramework("CoreFoundation")]),
        .target(
            name: "CMarathonSpatialOutput",
            dependencies: ["CMarathonSpatial"],
            path: ".",
            exclude: notEngine,
            sources: ["src/core/CoreAudioOutput.cpp"],
            publicHeadersPath: "include/CMarathonSpatialOutput",
            cxxSettings: [.headerSearchPath("include")],
            linkerSettings: [.linkedFramework("CoreAudio")]),
        .target(
            name: "MarathonSpatial",
            dependencies: ["CMarathonSpatial", "CMarathonSpatialOutput"],
            path: "swift/MarathonSpatial",
            swiftSettings: [.swiftLanguageMode(.v6)],
            linkerSettings: [.linkedFramework("CoreMotion")]),
        .target(
            name: "MarathonSpatialHRTFData",
            path: "data",
            sources: ["HRTFData.swift"],
            resources: [.copy("hrtf")],
            swiftSettings: [.swiftLanguageMode(.v6)]),
        .testTarget(
            name: "MarathonSpatialTests",
            dependencies: ["MarathonSpatial", "MarathonSpatialHRTFData"],
            path: "tests/swift",
            swiftSettings: [.swiftLanguageMode(.v6)]),
    ],
    cxxLanguageStandard: .cxx20
)
