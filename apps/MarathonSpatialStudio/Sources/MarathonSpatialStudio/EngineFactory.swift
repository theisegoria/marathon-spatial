import StudioKit

/// The one place that chooses the engine. The package defines MARATHON_SPATIAL_LINKED when it finds the built library.
enum EngineFactory {
    static func make() -> SpatialEngine {
        #if MARATHON_SPATIAL_LINKED
        LibraryEngine()
        #else
        UnavailableEngine(reason: "The Marathon Spatial library is not linked into this build, so playback, strikes and export are off.")
        #endif
    }
}
