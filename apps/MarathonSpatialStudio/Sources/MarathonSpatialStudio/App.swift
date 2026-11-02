import StudioKit
import SwiftUI
import UniformTypeIdentifiers

extension UTType {
    /// Declared in Support/Info.plist; the bundle script puts it in the app.
    static let spatialProject = UTType(exportedAs: "org.marathonspatial.project", conformingTo: .package)
}

struct ProjectDocument: FileDocument {
    static var readableContentTypes: [UTType] { [.spatialProject] }

    var project: SpatialProject

    init(project: SpatialProject = .starter()) { self.project = project }

    init(configuration: ReadConfiguration) throws {
        project = try ProjectPackage.read(configuration.file, store: .shared)
    }

    func fileWrapper(configuration: WriteConfiguration) throws -> FileWrapper {
        try ProjectPackage.write(project, store: .shared)
    }
}

@main
struct MarathonSpatialStudioApp: App {
    var body: some Scene {
        DocumentGroup(newDocument: ProjectDocument()) { file in
            ProjectWindow(document: file.$document)
        }
        .defaultSize(width: 1280, height: 820)
        .commands { StudioCommands() }

        Settings { PreferencesView() }
    }
}

struct PreferencesView: View {
    @AppStorage("reflectionPathLimit") private var pathLimit = 64
    @AppStorage("meterRefreshHz") private var meterHz = 30.0

    var body: some View {
        Form {
            Stepper("Reflection paths drawn: \(pathLimit)", value: $pathLimit, in: 8...256, step: 8)
            Picker("Meter refresh", selection: $meterHz) {
                Text("15 per second").tag(15.0)
                Text("30 per second").tag(30.0)
                Text("60 per second").tag(60.0)
            }
        }
        .formStyle(.grouped)
        .frame(width: 380)
    }
}

struct StudioCommands: Commands {
    @FocusedValue(\.projectSession) private var session

    var body: some Commands {
        CommandMenu("Transport") {
            Button(session?.isPlaying == true ? "Stop" : "Play") { session?.togglePlayback() }
                .keyboardShortcut("r", modifiers: .command)
                .disabled(session?.engine.status.canPlay != true)
            Button("Return to Start") { session?.seek(to: 0) }
                .keyboardShortcut(.leftArrow, modifiers: [.command, .option])
                .disabled(session == nil)
            Toggle("Record Strikes", isOn: Binding(get: { session?.isRecordingStrikes ?? false },
                                                   set: { session?.isRecordingStrikes = $0 }))
                .keyboardShortcut("r", modifiers: [.command, .shift])
                .disabled(session == nil)
            Divider()
            Button("Export Audio…") { session?.isExporting = true }
                .keyboardShortcut("e", modifiers: .command)
                .disabled(session == nil)
        }
        CommandMenu("Scene") {
            Button("Import Audio…") { session?.isImporting = true }
                .keyboardShortcut("i", modifiers: [.command, .shift])
                .disabled(session == nil)
            Button("Edit Materials…") { session?.isEditingMaterials = true }
                .keyboardShortcut("m", modifiers: [.command, .shift])
                .disabled(session == nil)
            Divider()
            ForEach(Array(EditTool.allCases.enumerated()), id: \.element) { index, tool in
                Button(tool.title) { session?.tool = tool }
                    .keyboardShortcut(KeyEquivalent(Character("\(index + 1)")), modifiers: .command)
                    .disabled(session == nil)
            }
            Divider()
            Button("Select Next Object") { session?.cycleSelection() }
                .keyboardShortcut("]", modifiers: .command)
                .disabled(session == nil)
            Button("Delete Selection") { session?.requestDelete = true }
                .keyboardShortcut(.delete, modifiers: .command)
                .disabled(session == nil)
        }
        CommandGroup(after: .toolbar) {
            Button("Room") { session?.workspace = .room }.keyboardShortcut("1", modifiers: [.command, .option])
            Button("Playground") { session?.workspace = .playground }.keyboardShortcut("2", modifiers: [.command, .option])
            Toggle("Reflection Paths", isOn: Binding(get: { session?.showReflections ?? true },
                                                     set: { session?.showReflections = $0 }))
                .keyboardShortcut("p", modifiers: [.command, .shift])
            Button("Zoom In") { session?.zoom *= 1.25 }.keyboardShortcut("=", modifiers: .command)
            Button("Zoom Out") { session?.zoom /= 1.25 }.keyboardShortcut("-", modifiers: .command)
            Button("Fit Room") { session?.zoom = 1 }.keyboardShortcut("0", modifiers: .command)
        }
    }
}

extension FocusedValues {
    @Entry var projectSession: ProjectSession?
}
