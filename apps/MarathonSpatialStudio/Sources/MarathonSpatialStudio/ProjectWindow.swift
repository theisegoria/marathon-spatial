import StudioKit
import SwiftUI
import UniformTypeIdentifiers

struct ProjectWindow: View {
    @Binding var document: ProjectDocument
    @State private var session = ProjectSession()
    @State private var showInspector = true
    @AppStorage("reflectionPathLimit") private var pathLimit = 64

    private var project: Binding<SpatialProject> { $document.project }

    var body: some View {
        NavigationSplitView {
            SceneSidebar(project: project, session: session)
                .navigationSplitViewColumnWidth(min: 220, ideal: 260)
        } detail: {
            VStack(spacing: 0) {
                switch session.workspace {
                case .room: RoomWorkspace(project: project, session: session)
                case .playground: PlaygroundView(project: project, session: session)
                }
            }
            .safeAreaInset(edge: .bottom, spacing: 0) {
                // GlassEffectContainer takes all the height it is offered: left alone, it squeezes the room plan to nothing.
                TransportBar(project: project, session: session)
                    .fixedSize(horizontal: false, vertical: true)
                    .padding(10)
            }
            .overlay(alignment: .top) {
                if let limitation = session.engine.status.limitation {
                    EngineBanner(text: limitation).padding(.top, 8)
                }
            }
        }
        .inspector(isPresented: $showInspector) {
            InspectorView(project: project, session: session)
                .inspectorColumnWidth(min: 280, ideal: 320, max: 420)
        }
        .toolbar { toolbar }
        .focusedSceneValue(\.projectSession, session)
        .task { await session.start() }
        .onChange(of: document.project, initial: true) { _, new in
            session.projectChanged(new, pathLimit: pathLimit)
        }
        .onChange(of: session.showReflections) { session.projectChanged(document.project, pathLimit: pathLimit) }
        .onChange(of: session.maximumDisplayOrder) { session.projectChanged(document.project, pathLimit: pathLimit) }
        .onChange(of: session.requestDelete) { _, requested in
            if requested {
                deleteSelection()
                session.requestDelete = false
            }
        }
        .fileImporter(isPresented: $session.isImporting, allowedContentTypes: [.audio], allowsMultipleSelection: true) { result in
            switch result {
            case .success(let urls): session.importClips(urls, into: project)
            case .failure(let error): session.alert = error.localizedDescription
            }
        }
        .sheet(isPresented: $session.isExporting) { ExportSheet(project: document.project, session: session) }
        .sheet(isPresented: $session.isEditingMaterials) { MaterialEditor(project: project) }
        .alert("Marathon Spatial Studio", isPresented: Binding(get: { session.alert != nil }, set: { if !$0 { session.alert = nil } })) {
            Button("OK", role: .cancel) {}
        } message: {
            Text(session.alert ?? "")
        }
    }

    @ToolbarContentBuilder private var toolbar: some ToolbarContent {
        ToolbarItem(placement: .principal) {
            Picker("Workspace", selection: $session.workspace) {
                ForEach(Workspace.allCases) { Label($0.title, systemImage: $0.symbol).tag($0) }
            }
            .pickerStyle(.segmented)
            .fixedSize()
        }
        ToolbarItemGroup(placement: .primaryAction) {
            Button("Import Audio", systemImage: "square.and.arrow.down") { session.isImporting = true }
                .help("Bring in your own audio")
            Button("Export", systemImage: "square.and.arrow.up") { session.isExporting = true }
                .help("Render offline to a WAV file")
            Button("Inspector", systemImage: "sidebar.trailing") { showInspector.toggle() }
                .help("Show or hide the inspector")
        }
    }

    private func deleteSelection() {
        switch session.selection {
        case .source(let id)?: document.project.sources.removeAll { $0.id == id }
        case .body(let id)?:
            document.project.bodies.removeAll { $0.id == id }
            document.project.strikes.removeAll { $0.body == id }
        case .block(let id)?: document.project.room.blocks.removeAll { $0.id == id }
        case .bed(let id)?: document.project.beds.removeAll { $0.id == id }
        default: return
        }
        session.selection = .room
    }
}

struct EngineBanner: View {
    var text: String

    var body: some View {
        Label(text, systemImage: "exclamationmark.triangle")
            .font(.callout)
            .padding(.horizontal, 14)
            .padding(.vertical, 8)
            .glassPanel(cornerRadius: 12)
            .accessibilityLabel("Engine notice: \(text)")
    }
}

struct SceneSidebar: View {
    @Binding var project: SpatialProject
    var session: ProjectSession

    var body: some View {
        @Bindable var session = session
        List(selection: $session.selection) {
            Section("Scene") {
                Label("Room", systemImage: "square.dashed").tag(SceneSelection.room)
                Label("Listener", systemImage: "person.fill").tag(SceneSelection.listener)
            }
            Section("Sources") {
                ForEach(project.sources) { source in
                    SourceRow(source: source, level: session.meters.sources[source.id] ?? 0)
                        .tag(SceneSelection.source(source.id))
                        .contextMenu {
                            Button("Delete", role: .destructive) { project.sources.removeAll { $0.id == source.id } }
                        }
                }
                Button("Add Source", systemImage: "plus") { addSource() }
                    .buttonStyle(.borderless)
            }
            Section("Bodies") {
                ForEach(project.bodies) { body in
                    Label(body.name, systemImage: "hammer").tag(SceneSelection.body(body.id))
                }
                Button("Add Body", systemImage: "plus") { addBody() }
                    .buttonStyle(.borderless)
            }
            Section("Blocks") {
                ForEach(project.room.blocks) { block in
                    Label(block.name, systemImage: "cube").tag(SceneSelection.block(block.id))
                }
            }
            Section("Beds") {
                ForEach(project.beds) { bed in
                    Label(bed.name, systemImage: "globe").tag(SceneSelection.bed(bed.id))
                }
            }
            Section("Audio") {
                ForEach(project.clips) { clip in
                    ClipRow(clip: clip)
                        .contextMenu {
                            if AmbisonicBed.order(channels: clip.channels) != nil {
                                Button("Use as Ambisonic Bed") { addBed(clip) }
                            }
                            Button("Remove from Project", role: .destructive) { remove(clip) }
                        }
                }
                if project.clips.isEmpty {
                    Text("Drop audio files here, or use Import Audio.")
                        .font(.caption)
                        .foregroundStyle(.secondary)
                }
            }
        }
        .dropDestination(for: URL.self) { urls, _ in
            session.importClips(urls, into: $project)
            return !urls.isEmpty
        }
    }

    private func addSource() {
        let s = SoundSource(name: "Source \(project.sources.count + 1)", position: Vec3(plan: project.room.centroid, y: 1.5))
        project.sources.append(s)
        session.selection = .source(s.id)
    }

    private func addBed(_ clip: ClipReference) {
        let bed = AmbisonicBed(name: clip.name, clip: clip.id)
        project.beds.append(bed)
        session.selection = .bed(bed.id)
    }

    private func remove(_ clip: ClipReference) {
        project.clips.removeAll { $0.id == clip.id }
        project.beds.removeAll { $0.clip == clip.id }
        for i in project.sources.indices where project.sources[i].content == .clip(clip.id) { project.sources[i].content = .none }
    }

    private func addBody() {
        let b = SoundingBody(name: "Body \(project.bodies.count + 1)", material: project.physicalMaterials.first?.id ?? "steel",
                             position: Vec3(plan: project.room.centroid, y: 1), size: Vec3(0.6, 0.4, 0.003))
        project.bodies.append(b)
        session.selection = .body(b.id)
    }
}

private struct SourceRow: View {
    var source: SoundSource
    var level: Float

    var body: some View {
        HStack {
            Label(source.name, systemImage: source.muted ? "speaker.slash" : "speaker.wave.2")
            Spacer()
            if source.solo {
                Text("S").font(.caption.bold()).foregroundStyle(.yellow).accessibilityLabel("Solo")
            }
            Capsule()
                .fill(Palette.sourceTint)
                .frame(width: 30 * CGFloat(max(0, min(1, (decibels(level) + 60) / 60))), height: 4)
                .accessibilityHidden(true)
        }
    }
}

private struct ClipRow: View {
    var clip: ClipReference
    @State private var overview: [Float] = []

    var body: some View {
        VStack(alignment: .leading, spacing: 2) {
            Text(clip.name).lineLimit(1)
            Text("\(clip.durationSeconds, format: .number.precision(.fractionLength(1))) s, \(clip.channels) ch, \(Int(clip.sampleRate)) Hz")
                .font(.caption)
                .foregroundStyle(.secondary)
            Canvas { context, size in
                guard !overview.isEmpty else { return }
                let w = size.width / CGFloat(overview.count)
                for (i, peak) in overview.enumerated() {
                    let h = max(1, CGFloat(peak) * size.height)
                    context.fill(Path(CGRect(x: CGFloat(i) * w, y: (size.height - h) / 2, width: max(1, w - 0.5), height: h)), with: .color(.secondary))
                }
            }
            .frame(height: 22)
            .accessibilityHidden(true)
        }
        .task(id: clip.id) { overview = (try? await ClipStore.shared.overview(clip, buckets: 120)) ?? [] }
        .accessibilityElement(children: .combine)
    }
}

struct TransportBar: View {
    @Binding var project: SpatialProject
    var session: ProjectSession

    var body: some View {
        @Bindable var session = session
        GlassEffectContainer(spacing: 10) {
            HStack(spacing: 12) {
                HStack(spacing: 4) {
                    Button("Return to Start", systemImage: "backward.end.fill") { session.seek(to: 0) }
                        .labelStyle(.iconOnly)
                    Button(session.isPlaying ? "Stop" : "Play", systemImage: session.isPlaying ? "stop.fill" : "play.fill") {
                        session.togglePlayback()
                    }
                    .labelStyle(.iconOnly)
                    .disabled(!session.engine.status.canPlay)
                    Toggle("Record strikes", systemImage: "record.circle", isOn: $session.isRecordingStrikes)
                        .toggleStyle(.button)
                        .labelStyle(.iconOnly)
                        .tint(.red)
                        .help("Record strikes onto the timeline while playing")
                }
                .buttonStyle(.glass)
                .controlSize(.large)

                Text(timecode(session.playhead))
                    .font(.body.monospacedDigit())
                    .frame(width: 70, alignment: .trailing)
                    .accessibilityLabel("Playhead")
                    .accessibilityValue("\(Int(session.playhead)) seconds")

                TimelineScrubber(project: $project, session: session)

                Text(timecode(project.timelineSeconds))
                    .font(.body.monospacedDigit())
                    .foregroundStyle(.secondary)
                    .accessibilityHidden(true)

                OutputMeters(meters: session.meters)
                    .frame(width: 180)
            }
            .padding(.horizontal, 12)
            .padding(.vertical, 8)
            .glassPanel(cornerRadius: 22)
        }
    }

    private func timecode(_ t: Double) -> String {
        let total = Int(max(t, 0))
        return String(format: "%d:%02d.%d", total / 60, total % 60, Int(max(t, 0) * 10) % 10)
    }
}

struct TimelineScrubber: View {
    @Binding var project: SpatialProject
    var session: ProjectSession

    var body: some View {
        ZStack {
            Slider(value: Binding(get: { min(session.playhead, project.timelineSeconds) }, set: { session.seek(to: $0) }),
                   in: 0...max(project.timelineSeconds, 0.1))
                .accessibilityLabel("Timeline position")
            GeometryReader { geo in
                ForEach(project.strikes) { s in
                    Rectangle()
                        .fill(Palette.bodyTint)
                        .frame(width: 2, height: 10)
                        .offset(x: geo.size.width * CGFloat(s.timeSeconds / max(project.timelineSeconds, 0.1)), y: 2)
                }
            }
            .allowsHitTesting(false)
            .accessibilityHidden(true)
        }
    }
}

/// Stereo peak and RMS bars. dB readouts, so level is never shown by colour alone.
struct OutputMeters: View {
    var meters: MeterFrame

    var body: some View {
        VStack(spacing: 3) {
            bar(0, "L")
            bar(1, "R")
        }
        .accessibilityElement(children: .ignore)
        .accessibilityLabel("Output level")
        .accessibilityValue(String(format: "Peak %.0f dB left, %.0f dB right%@", decibels(meters.peak[0]), decibels(meters.peak[1]),
                                   meters.overs > 0 ? ", clipped" : ""))
    }

    private func bar(_ channel: Int, _ name: String) -> some View {
        let peak = position(meters.peak[channel]), rms = position(meters.rms[channel])
        return HStack(spacing: 4) {
            Text(name).font(.caption2.monospaced()).foregroundStyle(.secondary)
            GeometryReader { geo in
                ZStack(alignment: .leading) {
                    Capsule().fill(.quaternary)
                    Capsule().fill(Color.green.opacity(0.8)).frame(width: geo.size.width * rms)
                    Rectangle().fill(meters.peak[channel] >= 1 ? Color.red : Color.primary)
                        .frame(width: 2)
                        .offset(x: geo.size.width * peak)
                }
            }
            .frame(height: 6)
            Text(String(format: "%4.0f", max(-99, decibels(meters.peak[channel]))))
                .font(.caption2.monospacedDigit())
                .frame(width: 28, alignment: .trailing)
        }
    }

    /// Maps -60...0 dBFS onto the bar.
    private func position(_ linear: Float) -> CGFloat { CGFloat(max(0, min(1, (decibels(linear) + 60) / 60))) }
}
