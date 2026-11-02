import Charts
import StudioKit
import SwiftUI

/// Hit physical objects and hear them through the room. Every strike goes to the library's modal synthesis.
struct PlaygroundView: View {
    @Binding var project: SpatialProject
    var session: ProjectSession
    @FocusState private var focused: Bool

    var body: some View {
        @Bindable var session = session
        VStack(spacing: 0) {
            GlassEffectContainer(spacing: 10) {
                HStack(spacing: 16) {
                    control("Impulse") {
                        Slider(value: Binding(get: { log10(session.strikeImpulseNs) }, set: { session.strikeImpulseNs = pow(10, $0) }), in: -3...1)
                            .accessibilityLabel("Strike impulse")
                            .accessibilityValue("\(session.strikeImpulseNs, format: .number.precision(.significantDigits(2))) newton seconds")
                    }
                    Text("\(session.strikeImpulseNs, format: .number.precision(.significantDigits(2))) N s")
                        .font(.caption.monospacedDigit())
                        .frame(width: 66, alignment: .leading)
                        .accessibilityHidden(true)
                    control("Mallet") {
                        Slider(value: Binding(get: { -log10(session.strikeContactSeconds) },
                                              set: { session.strikeContactSeconds = pow(10, -$0) }),
                               in: -log10(Strike.contactRange.upperBound)...(-log10(Strike.contactRange.lowerBound)))
                            .accessibilityLabel("Mallet hardness")
                            .accessibilityValue("\(session.strikeContactSeconds * 1000, format: .number.precision(.significantDigits(2))) millisecond contact")
                    }
                    Text("\(session.strikeContactSeconds * 1000, format: .number.precision(.significantDigits(2))) ms contact")
                        .font(.caption.monospacedDigit())
                        .frame(width: 96, alignment: .leading)
                        .accessibilityHidden(true)
                    control("Press") {
                        Slider(value: $session.scrapeForceN, in: 0.5...50)
                            .accessibilityLabel("Scrape pressure")
                            .accessibilityValue("\(Int(session.scrapeForceN)) newtons")
                    }
                    Menu("Add Body", systemImage: "plus") {
                        ForEach(bodyPresets, id: \.0) { name, material, size in
                            Button(name) { add(name: name, material: material, size: size) }
                        }
                    }
                    .fixedSize()
                }
                .padding(.horizontal, 14)
                .padding(.vertical, 8)
                .glassPanel(cornerRadius: 16)
            }
            .fixedSize(horizontal: false, vertical: true)
            .padding(10)

            if project.bodies.isEmpty {
                ContentUnavailableView("No Bodies", systemImage: "hammer", description: Text("Add a body. Click to strike it, drag to scrape it, or press keys 1 to 9."))
                    .frame(maxHeight: .infinity)
            } else {
                ScrollView {
                    LazyVGrid(columns: [GridItem(.adaptive(minimum: 300), spacing: 12)], spacing: 12) {
                        ForEach(Array(project.bodies.enumerated()), id: \.element.id) { index, body in
                            BodyCard(index: index, item: body, project: $project, session: session)
                        }
                    }
                    .padding(.horizontal, 12)
                }
            }
            if !project.strikes.isEmpty {
                StrikeList(project: $project).frame(height: 120).padding(.horizontal, 12)
            }
        }
        .focusable()
        .focusEffectDisabled()
        .focused($focused)
        .onAppear { focused = true }
        .onKeyPress(characters: .decimalDigits) { press in
            guard let n = Int(press.characters), n >= 1, n <= project.bodies.count else { return .ignored }
            session.strike(project.bodies[n - 1], u: 0.37, v: 0.41, project: &project)
            return .handled
        }
        .onKeyPress(.space) {
            session.togglePlayback()
            return .handled
        }
    }

    /// A caption that keeps its width beside a slider that gives way when the window narrows.
    private func control(_ title: String, @ViewBuilder slider: () -> some View) -> some View {
        HStack(spacing: 6) {
            Text(title).fixedSize()
            slider().labelsHidden().frame(minWidth: 70, maxWidth: 150)
        }
    }

    private var bodyPresets: [(String, String, Vec3)] {
        [("Steel plate", "steel", Vec3(0.6, 0.4, 0.003)),
         ("Thunder sheet", "sheet-metal", Vec3(1.8, 0.9, 0.001)),
         ("Glass pane", "glass", Vec3(0.5, 0.35, 0.004)),
         ("Hardwood plank", "hardwood", Vec3(1.0, 0.1, 0.02)),
         ("Aluminium bar", "aluminium", Vec3(0.8, 0.03, 0.02)),
         ("Brass chime", "brass", Vec3(0.3, 0.02, 0.02)),
         ("Concrete block", "concrete", Vec3(0.4, 0.2, 0.2)),
         ("Ceramic tile", "ceramic", Vec3(0.3, 0.3, 0.008))]
    }

    private func add(name: String, material: String, size: Vec3) {
        var b = SoundingBody(name: name, material: material, position: Vec3(plan: project.room.centroid, y: 1.2), size: size)
        b.position.x += Double(project.bodies.count) * 0.6
        project.bodies.append(b)
        session.selection = .body(b.id)
    }
}

private struct BodyCard: View {
    var index: Int
    var item: SoundingBody
    @Binding var project: SpatialProject
    var session: ProjectSession
    @State private var modes: [BodyMode] = []
    @State private var lastHit: CGPoint?
    @State private var pressed = false
    @State private var scraping = false

    private var materialName: String { project.physicalMaterials.first { $0.id == item.material }?.name ?? item.material }

    var body: some View {
        VStack(alignment: .leading, spacing: 8) {
            HStack {
                Text("\(index + 1)").font(.caption.monospaced()).foregroundStyle(.secondary).accessibilityHidden(true)
                Text(item.name).font(.headline)
                Spacer()
                Text(materialName).font(.caption).foregroundStyle(.secondary)
            }
            GeometryReader { geo in
                let face = faceRect(in: geo.size)
                Canvas { c, _ in
                    c.fill(Path(roundedRect: face, cornerRadius: 3), with: .color(Palette.bodyTint.opacity(session.lastStruck == item.id ? 0.55 : 0.3)))
                    c.stroke(Path(roundedRect: face, cornerRadius: 3), with: .color(.primary.opacity(0.5)), lineWidth: 1)
                    if let lastHit {
                        c.stroke(Path(ellipseIn: CGRect(x: lastHit.x - 6, y: lastHit.y - 6, width: 12, height: 12)), with: .color(.primary), lineWidth: 1.5)
                    }
                }
                .contentShape(Rectangle())
                .gesture(DragGesture(minimumDistance: 0).onChanged { value in
                    press(value, face: face)
                }.onEnded { _ in
                    if scraping { session.engine.endScrape(body: item.id) }
                    scraping = false
                    pressed = false
                })
            }
            .frame(height: 130)
            .accessibilityElement()
            .accessibilityLabel("Strike \(item.name)")
            .accessibilityValue("\(materialName) \(item.modelName)")
            .accessibilityAddTraits(.isButton)
            .accessibilityAction { session.strike(item, u: 0.37, v: 0.41, project: &project) }
            .accessibilityAction(named: "Strike centre") { session.strike(item, u: 0.5, v: 0.5, project: &project) }
            .accessibilityAction(named: "Strike edge") { session.strike(item, u: 0.05, v: 0.5, project: &project) }
            .accessibilityAction(named: "Show in inspector") { session.selection = .body(item.id) }
            if modes.isEmpty {
                Text(session.engine.status.canStrike ? "Computing modes" : "Modes appear when the library is linked.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            } else {
                ModeChart(modes: modes).frame(height: 70)
            }
        }
        .padding(12)
        .glassPanel(cornerRadius: 18)
        .task(id: item) { modes = await session.engine.modes(of: item.id) }
    }

    /// Pressing strikes straight away, for the least latency. Moving on from there becomes a scrape.
    private func press(_ value: DragGesture.Value, face: CGRect) {
        guard face.contains(value.location) else { return }
        let u = (value.location.x - face.minX) / face.width, v = (value.location.y - face.minY) / face.height
        lastHit = value.location
        if !pressed {
            pressed = true
            session.strike(item, u: u, v: v, project: &project)
        } else if hypot(value.translation.width, value.translation.height) > 4 {
            scraping = true
            let metresPerPoint = item.sortedSize.x / max(face.width, 1)
            session.engine.scrape(body: item.id, u: u, v: v, speedMetresPerSecond: hypot(value.velocity.width, value.velocity.height) * metresPerPoint,
                                  normalForceN: session.scrapeForceN)
        }
    }

    /// The broad face, drawn to scale inside the card.
    private func faceRect(in size: CGSize) -> CGRect {
        let s = item.sortedSize
        let aspect = max(s.y, 1e-3) / max(s.x, 1e-3)
        let w = min(size.width - 8, (size.height - 8) / max(aspect, 0.05))
        let h = max(10, w * aspect)
        return CGRect(x: (size.width - w) / 2, y: (size.height - h) / 2, width: w, height: h)
    }
}

/// Modes from the library: frequency across, decay as height.
struct ModeChart: View {
    var modes: [BodyMode]

    var body: some View {
        Chart(modes) { mode in
            RuleMark(x: .value("Frequency", mode.frequencyHz), yStart: .value("Decay", 0.0), yEnd: .value("Decay", min(mode.decaySeconds, 30)))
                .foregroundStyle(Palette.bodyTint)
        }
        .chartXScale(type: .log)
        .chartXAxisLabel("Hz")
        .chartYAxisLabel("s to -60 dB")
        .accessibilityLabel("\(modes.count) modes, lowest \(Int(modes.first?.frequencyHz ?? 0)) hertz")
    }
}

private struct StrikeList: View {
    @Binding var project: SpatialProject

    var body: some View {
        VStack(alignment: .leading, spacing: 4) {
            HStack {
                Text("Recorded strikes").font(.headline)
                Spacer()
                Button("Clear", role: .destructive) { project.strikes.removeAll() }
            }
            List {
                ForEach(project.strikes.sorted { $0.timeSeconds < $1.timeSeconds }) { s in
                    HStack {
                        Text(String(format: "%.3f s", s.timeSeconds)).monospacedDigit()
                        Text(project.bodies.first { $0.id == s.body }?.name ?? "Removed body")
                        Spacer()
                        Text("\(s.impulseNs, format: .number.precision(.significantDigits(2))) N s").foregroundStyle(.secondary)
                        Button("Delete Strike", systemImage: "xmark") { project.strikes.removeAll { $0.id == s.id } }
                            .labelStyle(.iconOnly)
                            .buttonStyle(.borderless)
                    }
                }
            }
            .listStyle(.plain)
        }
    }
}
