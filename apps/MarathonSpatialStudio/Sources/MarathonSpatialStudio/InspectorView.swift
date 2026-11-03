import Charts
import StudioKit
import SwiftUI
import UniformTypeIdentifiers

struct InspectorView: View {
    @Binding var project: SpatialProject
    var session: ProjectSession

    var body: some View {
        Form {
            switch session.selection {
            case .listener?: ListenerInspector(listener: $project.listener, session: session)
            case .bed(let id)?:
                if let i = project.beds.firstIndex(where: { $0.id == id }) {
                    BedInspector(bed: $project.beds[i], clip: project.clip(project.beds[i].clip))
                }
            case .source(let id)?:
                if let i = project.sources.firstIndex(where: { $0.id == id }) {
                    SourceInspector(source: $project.sources[i], project: project, session: session)
                }
            case .body(let id)?:
                if let i = project.bodies.firstIndex(where: { $0.id == id }) {
                    BodyInspector(item: $project.bodies[i], materials: project.physicalMaterials, session: session)
                }
            case .block(let id)?:
                if let i = project.room.blocks.firstIndex(where: { $0.id == id }) {
                    BlockInspector(block: $project.room.blocks[i], materials: project.materials)
                }
            case .wall(let i)?:
                if project.room.walls.indices.contains(i) {
                    Section("Wall \(i + 1)") { MaterialPicker(selection: $project.room.walls[i], materials: project.materials) }
                    RoomInspector(project: $project, session: session)
                }
            case .floor?:
                Section("Floor") { MaterialPicker(selection: $project.room.floor, materials: project.materials) }
                RoomInspector(project: $project, session: session)
            case .ceiling?:
                Section("Ceiling") { MaterialPicker(selection: $project.room.ceiling, materials: project.materials) }
                RoomInspector(project: $project, session: session)
            default:
                RoomInspector(project: $project, session: session)
            }
            RenderInspector(project: $project, session: session)
        }
        .formStyle(.grouped)
    }
}

struct MaterialPicker: View {
    @Binding var selection: MaterialID
    var materials: [AcousticMaterial]
    var label = "Material"

    var body: some View {
        Picker(label, selection: $selection) {
            ForEach(materials) { m in
                Label { Text(m.name) } icon: { MaterialSwatch(material: m) }.tag(m.id)
            }
        }
    }
}

private struct RoomInspector: View {
    @Binding var project: SpatialProject
    var session: ProjectSession
    @State private var userPresets: [RoomPreset] = []
    @State private var measured: MeasuredRoom?

    var body: some View {
        let estimate = RoomAcoustics.reverb(room: project.room, materials: project.materials)
        Section("Room") {
            Menu("Load Room Preset") {
                ForEach(RoomPreset.builtIn + userPresets) { preset in
                    Button(preset.name) { load(preset) }
                }
                Divider()
                Button("Save This Room as Preset") { savePreset() }
            }
            MetresField(label: "Ceiling height", value: $project.room.height)
            MaterialPicker(selection: $project.room.ceiling, materials: project.materials, label: "Ceiling")
            MaterialPicker(selection: $project.room.floor, materials: project.materials, label: "Floor")
            LabeledContent("Volume") { Text("\(project.room.volume, format: .number.precision(.fractionLength(1))) m³") }
            Button("Edit Materials…") { session.isEditingMaterials = true }
        }
        Section("Reverberation estimate") {
            Chart(Array(zip(estimate.bandsHz, estimate.eyringSeconds)), id: \.0) { band, seconds in
                BarMark(x: .value("Band", band < 1000 ? "\(band.formatted())" : "\(Int(band / 1000))k"), y: .value("RT60", seconds))
            }
            .chartYAxisLabel("seconds")
            .frame(height: 110)
            .accessibilityLabel("RT60 by octave band")
            LabeledContent("Mid RT60 (Eyring)") { Text("\(estimate.midSeconds, format: .number.precision(.fractionLength(2))) s") }
            if let measured {
                LabeledContent("Mid RT60 (library)") { Text("\(measured.midSeconds, format: .number.precision(.fractionLength(2))) s") }
                if let traced = measured.tracedMidSeconds {
                    LabeledContent("Mid RT60 (traced)") { Text("\(traced, format: .number.precision(.fractionLength(2))) s") }
                }
            }
            Text("Statistical estimate from volume and absorption. The traced figure comes from the energy field the convolution reverb is built from, and settles over a few seconds after the room changes.")
                .font(.caption)
                .foregroundStyle(.secondary)
        }
        .task(id: project.room) {
            // The traced field keeps converging after a change, so keep reading while the room is on screen.
            while !Task.isCancelled {
                measured = await session.engine.measuredRoom()
                try? await Task.sleep(for: .seconds(1))
            }
        }
        .task { await reload() }
    }

    private func reload() async { userPresets = await store.load().rooms }

    private var store: PresetStore {
        PresetStore(directory: FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask)[0]
            .appendingPathComponent("Marathon Spatial Studio"))
    }

    private func load(_ preset: RoomPreset) {
        // Materials a preset names must exist: restore any the project lost.
        let used = Set(preset.room.walls + [preset.room.floor, preset.room.ceiling] + preset.room.blocks.map(\.material))
        for id in used where project.material(id) == nil {
            if let m = MaterialLibrary.acoustic.first(where: { $0.id == id }) { project.materials.append(m) }
        }
        project.room = preset.room
        session.selection = .room
    }

    private func savePreset() {
        let preset = RoomPreset(id: UUID().uuidString, name: "Room \(userPresets.count + 1)", room: project.room, isUser: true)
        Task {
            do {
                try await store.add(room: preset)
                await reload()
            } catch { session.alert = error.localizedDescription }
        }
    }
}

private struct BedInspector: View {
    @Binding var bed: AmbisonicBed
    var clip: ClipReference?

    var body: some View {
        Section("Ambisonic bed") {
            TextField("Name", text: $bed.name)
            if let clip, let order = AmbisonicBed.order(channels: clip.channels) {
                LabeledContent("Order") { Text("\(order) (\(clip.channels) channels, ACN/SN3D)") }
            }
            LabeledContent("Gain") {
                Slider(value: $bed.gainDB, in: -60...12, step: 0.5) { Text("Gain") }
                    .labelsHidden()
                Text("\(bed.gainDB, format: .number.precision(.fractionLength(1))) dB").monospacedDigit().frame(width: 62, alignment: .trailing)
            }
            LabeledContent("Rotation") {
                Slider(value: $bed.yawDegrees, in: -180...180, step: 1) { Text("Rotation") }
                    .labelsHidden()
                Text("\(Int(bed.yawDegrees))°").monospacedDigit().frame(width: 44, alignment: .trailing)
            }
            Toggle("Loop", isOn: $bed.loops)
            Toggle("Locked to the head", isOn: $bed.headLocked)
        }
    }
}

private struct ListenerInspector: View {
    @Binding var listener: Listener
    var session: ProjectSession

    var body: some View {
        @Bindable var session = session
        Section("Head tracking") {
            Toggle("Track headphones", isOn: $session.headTracking)
                .disabled(!session.headTracker.isAvailable)
            Button("Recenter") { session.headTracker.recenter() }
                .disabled(!session.headTracking)
            if !session.headTracker.isAvailable {
                Text("No headphones with motion sensors are reported.").font(.caption).foregroundStyle(.secondary)
            }
        }
        Section("Listener") {
            Vec3Fields(value: $listener.position)
            LabeledContent("Heading") {
                Slider(value: $listener.yawDegrees, in: -180...180, step: 1) { Text("Heading") }
                    .labelsHidden()
                Text("\(Int(listener.yawDegrees))°").monospacedDigit().frame(width: 44, alignment: .trailing)
            }
            LabeledContent("Pitch") {
                Slider(value: $listener.pitchDegrees, in: -90...90, step: 1) { Text("Pitch") }
                    .labelsHidden()
                Text("\(Int(listener.pitchDegrees))°").monospacedDigit().frame(width: 44, alignment: .trailing)
            }
        }
        MotionSection(motion: $listener.motion, home: listener.position)
    }
}

private struct SourceInspector: View {
    @Binding var source: SoundSource
    var project: SpatialProject
    var session: ProjectSession

    var body: some View {
        Section("Source") {
            TextField("Name", text: $source.name)
            Picker("Audio", selection: $source.content) {
                Text("None").tag(SourceContent.none)
                Section("Your audio") {
                    ForEach(project.clips) { Text($0.name).tag(SourceContent.clip($0.id)) }
                }
                Section("Test signals") {
                    ForEach(TestSignal.allCases) { Text($0.title).tag(SourceContent.signal($0)) }
                }
            }
            Button("Import Audio…") { session.isImporting = true }
            Toggle("Loop", isOn: $source.loops)
            MetresField(label: "Starts at", value: $source.startSeconds, unit: "s")
            LabeledContent("Gain") {
                Slider(value: $source.gainDB, in: -60...12, step: 0.5) { Text("Gain") }
                    .labelsHidden()
                Text("\(source.gainDB, format: .number.precision(.fractionLength(1))) dB").monospacedDigit().frame(width: 62, alignment: .trailing)
            }
            HStack {
                Toggle("Mute", isOn: $source.muted)
                Toggle("Solo", isOn: $source.solo)
            }
        }
        Section("Position") {
            Vec3Fields(value: $source.position)
        }
        MotionSection(motion: $source.motion, home: source.position)
        Section("Directivity") {
            HStack(alignment: .center, spacing: 12) {
                PolarPlot(directivity: source.directivity).frame(width: 84, height: 84)
                VStack {
                    Slider(value: $source.directivity.dipoleWeight, in: 0...1) { Text("Omni to figure-8") }
                    Slider(value: $source.directivity.power, in: 0.5...8) { Text("Focus") }
                }
            }
            Slider(value: $source.directivity.yawDegrees, in: -180...180, step: 1) { Text("Aim") }
            Slider(value: $source.directivity.pitchDegrees, in: -90...90, step: 1) { Text("Tilt") }
        }
        Section("Propagation") {
            Toggle("Doppler", isOn: $source.doppler)
            Toggle("Air absorption", isOn: $source.airAbsorption)
            Toggle("Near field", isOn: $source.nearField)
            Toggle("Reflections", isOn: $source.reflections)
        }
    }
}

private struct PolarPlot: View {
    var directivity: Directivity

    var body: some View {
        Canvas { c, size in
            let m = CGPoint(x: size.width / 2, y: size.height / 2), r = min(size.width, size.height) / 2 - 2
            c.stroke(Path(ellipseIn: CGRect(x: m.x - r, y: m.y - r, width: 2 * r, height: 2 * r)), with: .color(.secondary.opacity(0.5)))
            var lobe = Path()
            for a in stride(from: 0.0, through: 360.0, by: 4.0) {
                let g = directivity.gain(offAxisDegrees: a) * r, t = a * .pi / 180
                let q = CGPoint(x: m.x + sin(t) * g, y: m.y - cos(t) * g)
                if a == 0 { lobe.move(to: q) } else { lobe.addLine(to: q) }
            }
            lobe.closeSubpath()
            c.fill(lobe, with: .color(Palette.sourceTint.opacity(0.5)))
            c.stroke(lobe, with: .color(Palette.sourceTint))
        }
        .accessibilityHidden(true)
    }
}

private struct MotionSection: View {
    @Binding var motion: Motion
    var home: Vec3

    private enum Kind: String, CaseIterable, Identifiable {
        case still, orbit, shuttle
        var id: String { rawValue }
    }

    private var kind: Binding<Kind> {
        Binding(get: {
            switch motion {
            case .still: .still
            case .orbit: .orbit
            case .shuttle: .shuttle
            }
        }, set: {
            switch $0 {
            case .still: motion = .still
            case .orbit: motion = .orbit(center: home.plan, radius: 2, periodSeconds: 8)
            case .shuttle: motion = .shuttle(offset: Vec3(4, 0, 0), periodSeconds: 6)
            }
        })
    }

    var body: some View {
        Section("Motion") {
            Picker("Path", selection: kind) {
                ForEach(Kind.allCases) { Text($0.rawValue.capitalized).tag($0) }
            }
            .pickerStyle(.segmented)
            switch motion {
            case .still: EmptyView()
            case let .orbit(center, radius, period):
                MetresField(label: "Radius", value: Binding(get: { radius }, set: { motion = .orbit(center: center, radius: $0, periodSeconds: period) }))
                MetresField(label: "Period", value: Binding(get: { period }, set: { motion = .orbit(center: center, radius: radius, periodSeconds: max(0.1, $0)) }), unit: "s")
            case let .shuttle(offset, period):
                Vec3Fields(labels: ["Travel X", "Travel Y", "Travel Z"], value: Binding(get: { offset }, set: { motion = .shuttle(offset: $0, periodSeconds: period) }))
                MetresField(label: "Period", value: Binding(get: { period }, set: { motion = .shuttle(offset: offset, periodSeconds: max(0.1, $0)) }), unit: "s")
            }
        }
    }
}

private struct BodyInspector: View {
    @Binding var item: SoundingBody
    var materials: [PhysicalMaterial]
    var session: ProjectSession

    var body: some View {
        Section("Body") {
            TextField("Name", text: $item.name)
            Picker("Material", selection: $item.material) {
                ForEach(materials) { Text($0.name).tag($0.id) }
            }
            Vec3Fields(labels: ["Length", "Width", "Thickness"], value: $item.size)
            MetresField(label: "Wall", value: $item.wallThicknessM)
            LabeledContent("Modelled as") { Text(item.modelName) }
            LabeledContent("Mass") {
                let density = materials.first { $0.id == item.material }?.densityKgM3 ?? 0
                Text("\(item.massKg(density: density), format: .number.precision(.significantDigits(3))) kg")
            }
            LabeledContent("Gain") {
                Slider(value: $item.gainDB, in: -60...12, step: 0.5) { Text("Gain") }
                    .labelsHidden()
                Text("\(item.gainDB, format: .number.precision(.fractionLength(1))) dB").monospacedDigit().frame(width: 62, alignment: .trailing)
            }
            Button("Strike") {
                session.engine.strike(Strike(body: item.id, timeSeconds: session.playhead, u: 0.37, v: 0.41,
                                             impulseNs: session.strikeImpulseNs, contactSeconds: session.strikeContactSeconds))
            }
            .disabled(!session.engine.status.canStrike)
        }
        Section("Position") { Vec3Fields(value: $item.position) }
    }
}

private struct BlockInspector: View {
    @Binding var block: Block
    var materials: [AcousticMaterial]

    var body: some View {
        Section("Block") {
            TextField("Name", text: $block.name)
            MaterialPicker(selection: $block.material, materials: materials)
        }
        Section("Centre") { Vec3Fields(value: $block.center) }
        Section("Size") { Vec3Fields(labels: ["Width", "Height", "Depth"], value: $block.size) }
    }
}

private struct RenderInspector: View {
    @Binding var project: SpatialProject
    var session: ProjectSession
    @State private var userPresets: [RenderPreset] = []

    private var store: PresetStore {
        PresetStore(directory: FileManager.default.urls(for: .applicationSupportDirectory, in: .userDomainMask)[0]
            .appendingPathComponent("Marathon Spatial Studio"))
    }

    var body: some View {
        Section("Render") {
            Menu("Load Render Preset") {
                ForEach(RenderPreset.builtIn + userPresets) { preset in
                    Button(preset.name) {
                        var next = preset.settings
                        next.hrtfSet = project.render.hrtfSet
                        project.render = next
                    }
                }
                Divider()
                Button("Save These Settings as Preset") { save() }
            }
            Picker("Output", selection: $project.render.output) {
                ForEach(OutputMode.allCases) { Text($0.title).tag($0) }
            }
            Picker("HRTF set", selection: $project.render.hrtfSet) {
                Text("Library default").tag(String?.none)
                ForEach(session.hrtfSets) { set in
                    Text(set.name).tag(Optional(set.id))
                }
            }
            .disabled(project.render.output != .binaural)
            if let id = project.render.hrtfSet, let set = session.hrtfSets.first(where: { $0.id == id }) {
                Text("\(set.detail). Licence: \(set.licence).").font(.caption).foregroundStyle(.secondary)
            }
            Button("Import SOFA File…") { session.isImportingHRTF = true }
                .disabled(!session.engine.status.canRender)
                .fileImporter(isPresented: Binding(get: { session.isImportingHRTF }, set: { session.isImportingHRTF = $0 }),
                              allowedContentTypes: [UTType(filenameExtension: "sofa") ?? .data]) { result in
                    importHRTF(result)
                }
            Text("Turn off your headphones' own spatial audio mode: this is already binaural.")
                .font(.caption)
                .foregroundStyle(.secondary)
            Toggle("Early reflections", isOn: $project.render.earlyReflections)
            Stepper("Drawn reflection order \(project.render.reflectionOrder)", value: $project.render.reflectionOrder, in: 1...3)
                .disabled(!project.render.earlyReflections)
            Group {
                Toggle("Late reverb", isOn: $project.render.lateReverb)
                HStack {
                    Text("Reverb trim").fixedSize()
                    Slider(value: $project.render.reverbTrimDB, in: -24...12, step: 0.5) { Text("Reverb trim") }
                        .labelsHidden()
                    Text("\(project.render.reverbTrimDB, format: .number.precision(.fractionLength(1))) dB").monospacedDigit().frame(width: 62, alignment: .trailing)
                }
            }
            .disabled(!session.engine.status.hasLateReverb)
            if !session.engine.status.hasLateReverb {
                Text("The library is not linked, so nothing here renders late reverb.")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
            Toggle("Peak limiter", isOn: $project.render.limiter)
            LabeledContent("Master") {
                Slider(value: $project.render.masterGainDB, in: -60...12, step: 0.5) { Text("Master") }
                    .labelsHidden()
                Text("\(project.render.masterGainDB, format: .number.precision(.fractionLength(1))) dB").monospacedDigit().frame(width: 62, alignment: .trailing)
            }
            MetresField(label: "Timeline length", value: $project.timelineSeconds, unit: "s")
        }
        .task { userPresets = await store.load().renders }
    }

    private func importHRTF(_ result: Result<URL, Error>) {
        Task {
            do {
                let set = try await session.engine.importHRTF(result.get())
                session.hrtfSets = await session.engine.hrtfSets()
                project.render.hrtfSet = set.id
            } catch { session.alert = error.localizedDescription }
        }
    }

    private func save() {
        let preset = RenderPreset(id: UUID().uuidString, name: "Render \(userPresets.count + 1)", settings: project.render, isUser: true)
        Task {
            do {
                try await store.add(render: preset)
                userPresets = await store.load().renders
            } catch { session.alert = error.localizedDescription }
        }
    }
}

/// The project's acoustic materials, editable band by band so measured data can go in.
struct MaterialEditor: View {
    @Binding var project: SpatialProject
    @State private var selected: MaterialID?
    @Environment(\.dismiss) private var dismiss

    var body: some View {
        NavigationSplitView {
            List(project.materials, selection: $selected) { m in
                Label { Text(m.name) } icon: { MaterialSwatch(material: m) }.tag(m.id)
            }
            .safeAreaInset(edge: .bottom) {
                Button("New Material", systemImage: "plus") { add() }.padding(8)
            }
            .navigationSplitViewColumnWidth(220)
        } detail: {
            if let i = project.materials.firstIndex(where: { $0.id == selected }) {
                Form {
                    TextField("Name", text: $project.materials[i].name)
                    Section("Absorption coefficient") {
                        ForEach(AcousticMaterial.bandsHz.indices, id: \.self) { b in
                            let hz = AcousticMaterial.bandsHz[b]
                            LabeledContent(hz < 1000 ? "\(hz.formatted()) Hz" : "\(Int(hz / 1000)) kHz") {
                                Slider(value: $project.materials[i].absorption[b], in: 0...1) { Text("\(Int(hz)) hertz") }
                                    .labelsHidden()
                                TextField("Value", value: $project.materials[i].absorption[b], format: .number.precision(.fractionLength(2)))
                                    .labelsHidden()
                                    .frame(width: 54)
                            }
                        }
                    }
                    Section("Surface") {
                        LabeledContent("Scattering") {
                            Slider(value: $project.materials[i].scattering, in: 0...1) { Text("Scattering") }
                                .labelsHidden()
                            Text("\(project.materials[i].scattering, format: .number.precision(.fractionLength(2)))").monospacedDigit()
                        }
                        MetresField(label: "Transmission", value: $project.materials[i].transmissionDB, unit: "dB")
                    }
                }
                .formStyle(.grouped)
            } else {
                ContentUnavailableView("Select a Material", systemImage: "square.3.layers.3d")
            }
        }
        .frame(minWidth: 680, minHeight: 500)
        .toolbar {
            ToolbarItem(placement: .confirmationAction) { Button("Done") { dismiss() } }
        }
    }

    private func add() {
        let id = UUID().uuidString
        project.materials.append(AcousticMaterial(id: id, name: "Custom \(project.materials.count + 1)", absorption: Array(repeating: 0.1, count: AcousticMaterial.bandsHz.count),
                                                  scattering: 0.1, transmissionDB: -30, swatch: project.materials.count))
        selected = id
    }
}
