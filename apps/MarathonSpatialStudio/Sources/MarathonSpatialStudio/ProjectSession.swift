import Observation
import StudioKit
import SwiftUI

enum Workspace: String, CaseIterable, Identifiable {
    case room, playground
    var id: String { rawValue }
    var title: String { self == .room ? "Room" : "Playground" }
    var symbol: String { self == .room ? "square.dashed" : "hammer" }
}

enum EditTool: String, CaseIterable, Identifiable {
    case select, addSource, addBlock, paint, outline, strike
    var id: String { rawValue }
    var title: String {
        switch self {
        case .select: "Select and Move"
        case .addSource: "Add Source"
        case .addBlock: "Add Block"
        case .paint: "Paint Material"
        case .outline: "Edit Outline"
        case .strike: "Strike Bodies"
        }
    }
    var symbol: String {
        switch self {
        case .select: "cursorarrow"
        case .addSource: "speaker.plus"
        case .addBlock: "cube"
        case .paint: "paintbrush.pointed"
        case .outline: "skew"
        case .strike: "hammer"
        }
    }
    /// Pressed with the room focused.
    var key: Character {
        switch self {
        case .select: "v"
        case .addSource: "s"
        case .addBlock: "b"
        case .paint: "p"
        case .outline: "o"
        case .strike: "k"
        }
    }
}

enum SceneSelection: Hashable {
    case room
    case listener
    case source(UUID)
    case body(UUID)
    case block(UUID)
    case wall(Int)
    case floor
    case ceiling
    case bed(UUID)
}

/// State of one window that is not part of the document: tools, transport, meters and the engine.
@MainActor @Observable
final class ProjectSession {
    let engine: SpatialEngine

    var workspace = Workspace.room
    var tool = EditTool.select
    var selection: SceneSelection? = .room
    var brush: MaterialID = "concrete"
    var zoom = 1.0
    var showReflections = true
    var maximumDisplayOrder = 3
    var requestDelete = false

    var isImporting = false
    var isExporting = false
    var isEditingMaterials = false
    var alert: String?

    var hrtfSets: [HRTFSetInfo] = []
    var isPlaying = false
    var isRecordingStrikes = false
    var playhead = 0.0
    var meters = MeterFrame()
    var isConfiguring = false
    var configurationError: String?
    /// Counts finished engine configurations. Views that read from the live engine (a body's modes) key their
    /// queries on it: the engine learns of a new body only once its configuration lands, after the view first asks.
    private(set) var configurationGeneration = 0

    var reflectionPaths: [ReflectionPath] = []
    var isImportingHRTF = false
    var headTracking = false {
        didSet {
            guard headTracking != oldValue else { return }
            if headTracking { headTracker.start(engine: engine) } else { headTracker.stop(engine: engine) }
        }
    }
    let headTracker = HeadTracker()
    /// Headphone head tracking belongs to the binaural output alone.
    var outputIsHeadphones: Bool { (project?.render.output ?? .binaural) == .binaural }
    let cameraTracker = CameraHeadTracker()
    private var cameraRunning = false

    /// The camera runs while the document asks for 3D speakers that follow the head, and only then.
    private func syncCamera(_ render: RenderSettings) {
        let want = render.output == .speakers3D && render.crosstalkSettings.headTracking
        guard want != cameraRunning else { return }
        cameraRunning = want
        if want { cameraTracker.start(engine: engine) } else { cameraTracker.stop(engine: engine) }
    }

    var strikeImpulseNs = 0.05
    var strikeContactSeconds = 0.001
    var scrapeForceN = 5.0
    var lastStruck: UUID?

    /// The newest document, for tick-driven updates and recording.
    private var project: SpatialProject?
    private var configuredProject: SpatialProject?
    private var configureTask: Task<Void, Never>?
    private var reflectionTask: Task<Void, Never>?
    private var playTask: Task<Void, Never>?
    private var meterTask: Task<Void, Never>?
    private var playbackOrigin = ContinuousClock.now
    private var playbackStart = 0.0
    private var pathLimit = 64
    private var seekTask: Task<Void, Never>?

    init(engine: SpatialEngine = EngineFactory.make()) {
        self.engine = engine
    }

    func start() async {
        hrtfSets = await engine.hrtfSets()
    }

    // MARK: Project changes

    /// Called on every document or drag-preview change.
    func projectChanged(_ project: SpatialProject, pathLimit: Int) {
        self.project = project
        self.pathLimit = pathLimit
        syncCamera(project.render)
        // Headphone head tracking turns the binaural image with the head; speakers stand still, so it stops for them.
        if project.render.output != .binaural, headTracking { headTracking = false }
        refreshReflections(project)
        if isPlaying { engine.update(ProjectCompiler.frame(project: project, time: playhead)) }
        guard needsReconfigure(project) else { return }
        // The limiter is fixed when an engine is created, so changing it rebuilds the engine and closes its output.
        if isPlaying, let c = configuredProject, c.render.limiter != project.render.limiter { stop() }
        scheduleConfigure(project)
    }

    /// Poses reach the engine through `update`; only the things below need a reload.
    private func needsReconfigure(_ p: SpatialProject) -> Bool {
        guard let c = configuredProject else { return true }
        return p.room != c.room || p.materials != c.materials || p.physicalMaterials != c.physicalMaterials || p.bodies != c.bodies
            || p.render != c.render || p.clips != c.clips || p.beds != c.beds
            || p.audibleSources.map(\.setupKey) != c.audibleSources.map(\.setupKey)
    }

    private func scheduleConfigure(_ project: SpatialProject) {
        configureTask?.cancel()
        let engine = engine
        configureTask = Task { [weak self] in
            try? await Task.sleep(for: .milliseconds(120))
            guard !Task.isCancelled else { return }
            self?.isConfiguring = true
            defer { self?.isConfiguring = false }
            do {
                let configuration = try await ProjectCompiler.configuration(project: project, store: .shared, sampleRate: engine.preferredSampleRate)
                guard !Task.isCancelled else { return }
                try await engine.configure(configuration)
                self?.configuredProject = project
                self?.configurationGeneration += 1
                self?.configurationError = nil
            } catch is CancellationError {
            } catch {
                self?.configurationError = error.localizedDescription
            }
        }
    }

    /// What the room draws. The library's simulator traces what is heard, on its own.
    private func refreshReflections(_ project: SpatialProject) {
        reflectionTask?.cancel()
        let limit = pathLimit
        guard project.render.earlyReflections else {
            reflectionPaths = []
            return
        }
        let order = min(project.render.reflectionOrder, 3)
        let time = playhead
        let draw = showReflections
        reflectionTask = Task.detached(priority: .userInitiated) { [weak self] in
            let frame = ProjectCompiler.frame(project: project, time: time)
            var all: [ReflectionPath] = []
            for source in project.audibleSources {
                guard source.reflections, let pose = frame.sources.first(where: { $0.id == source.id }) else { continue }
                if Task.isCancelled { return }
                let paths = ImageSources.paths(room: project.room, materials: project.materials, source: source, sourcePosition: pose.position,
                                               listener: frame.listenerPosition, maxOrder: order, limit: 256)
                all += paths
            }
            all.sort { $0.levelDB > $1.levelDB }
            let shown = draw ? Array(all.prefix(limit)) : []
            if Task.isCancelled { return }
            await MainActor.run { [weak self] in
                guard !Task.isCancelled else { return }
                self?.reflectionPaths = shown
            }
        }
    }

    // MARK: Transport

    func togglePlayback() {
        if isPlaying { stop() } else { play() }
    }

    func play() {
        guard engine.status.canPlay, let project else { return }
        let engine = engine
        let from = playhead
        playTask = Task { [weak self] in
            do {
                engine.update(ProjectCompiler.frame(project: project, time: from))
                try await engine.play(from: from)
            } catch {
                self?.alert = error.localizedDescription
                return
            }
            guard let self else { return }
            self.isPlaying = true
            self.playbackOrigin = .now
            self.playbackStart = from
            self.firedUpTo = from
            await self.runClock()
        }
    }

    private func runClock() async {
        let hz = max(15, UserDefaults.standard.double(forKey: "meterRefreshHz"))
        var tick = 0
        while isPlaying, !Task.isCancelled {
            let elapsed = playbackOrigin.duration(to: .now)
            let t = playbackStart + Double(elapsed.components.seconds) + Double(elapsed.components.attoseconds) * 1e-18
            playhead = t
            if let project {
                if t >= project.timelineSeconds {
                    stop()
                    playhead = 0
                    break
                }
                engine.update(ProjectCompiler.frame(project: project, time: t))
                fireRecordedStrikes(project, upTo: t)
                // Scripted motion changes the paths with no document edit. About six times a second is enough to draw.
                tick += 1
                if Int(hz / 6) > 0, tick % Int(hz / 6) == 0, project.hasMotion { refreshReflections(project) }
            }
            meters = engine.meters()
            try? await Task.sleep(for: .seconds(1 / hz))
        }
    }

    private var firedUpTo = 0.0
    private func fireRecordedStrikes(_ project: SpatialProject, upTo t: Double) {
        for s in project.strikes where s.timeSeconds > firedUpTo && s.timeSeconds <= t { engine.strike(s) }
        firedUpTo = t
    }

    func stop() {
        isPlaying = false
        playTask?.cancel()
        let engine = engine
        Task { await engine.stop() }
        meters = MeterFrame()
        firedUpTo = playhead
    }

    func seek(to t: Double) {
        let resume = isPlaying || seekTask != nil
        if isPlaying { stop() }
        playhead = max(0, t)
        firedUpTo = playhead
        if let project { projectChanged(project, pathLimit: pathLimit) }
        // A scrub is a stream of seeks: restart once, after it settles.
        seekTask?.cancel()
        seekTask = nil
        guard resume else { return }
        seekTask = Task { [weak self] in
            try? await Task.sleep(for: .milliseconds(150))
            guard !Task.isCancelled else { return }
            self?.seekTask = nil
            self?.play()
        }
    }

    // MARK: Editing helpers

    /// Strike a body now through the library. With recording on, the strike goes on the timeline too.
    func strike(_ body: SoundingBody, u: Double, v: Double, project: inout SpatialProject) {
        let s = Strike(body: body.id, timeSeconds: playhead, u: u, v: v, impulseNs: strikeImpulseNs, contactSeconds: strikeContactSeconds)
        engine.strike(s)
        lastStruck = body.id
        if isRecordingStrikes && isPlaying { project.strikes.append(s) }
    }

    func cycleSelection() {
        guard let project else { return }
        var order: [SceneSelection] = [.listener]
        order += project.sources.map { .source($0.id) }
        order += project.bodies.map { .body($0.id) }
        order += project.room.blocks.map { .block($0.id) }
        order += project.beds.map { .bed($0.id) }
        guard let current = selection, let i = order.firstIndex(of: current) else {
            selection = order.first
            return
        }
        selection = order[(i + 1) % order.count]
    }

    func importClips(_ urls: [URL], into document: Binding<SpatialProject>) {
        Task {
            for url in urls {
                do {
                    let clip = try await ClipStore.shared.importFile(url)
                    if let existing = document.wrappedValue.clips.first(where: { $0.sha256 == clip.sha256 }) {
                        attach(existing, to: document)
                    } else {
                        document.wrappedValue.clips.append(clip)
                        attach(clip, to: document)
                    }
                } catch {
                    alert = error.localizedDescription
                }
            }
        }
    }

    /// A newly imported clip goes straight to the selected source, or to a new one.
    private func attach(_ clip: ClipReference, to document: Binding<SpatialProject>) {
        if case let .source(id)? = selection, let i = document.wrappedValue.sources.firstIndex(where: { $0.id == id }),
           document.wrappedValue.sources[i].content == .none {
            document.wrappedValue.sources[i].content = .clip(clip.id)
            return
        }
        let center = document.wrappedValue.room.centroid
        var s = SoundSource(name: clip.name, position: Vec3(plan: center, y: 1.5), content: .clip(clip.id))
        s.position.x += Double(document.wrappedValue.sources.count) * 0.4
        document.wrappedValue.sources.append(s)
        selection = .source(s.id)
    }
}

private extension SpatialProject {
    var hasMotion: Bool { listener.motion != .still || sources.contains { $0.motion != .still } }
}

private extension SoundSource {
    /// What a source's setup depends on: everything except the pose fields sent per frame.
    var setupKey: SoundSource {
        var s = self
        s.position = .zero
        s.motion = .still
        return s
    }
}
