import AppKit
import StudioKit
import SwiftUI
import UniformTypeIdentifiers

/// Offline render to WAV. It runs faster than real time and opens no output device.
struct ExportSheet: View {
    var project: SpatialProject
    var session: ProjectSession

    @State private var options = ExportOptions()
    @State private var output = OutputMode.binaural
    @State private var fraction: Double?
    @State private var task: Task<Void, Never>?
    @State private var result: ExportResult?
    @State private var failure: String?
    @Environment(\.dismiss) private var dismiss

    var body: some View {
        VStack(alignment: .leading, spacing: 12) {
            Text("Export Audio").font(.title2.bold())
            Form {
                Picker("Mix", selection: $output) {
                    ForEach(OutputMode.allCases) { Text($0.title).tag($0) }
                }
                Picker("Format", selection: $options.format) {
                    ForEach(WAVFormat.allCases) { Text($0.title).tag($0) }
                }
                Picker("Sample rate", selection: $options.sampleRate) {
                    Text("44.1 kHz").tag(44100)
                    Text("48 kHz").tag(48000)
                    Text("96 kHz").tag(96000)
                }
                MetresField(label: "Start", value: $options.startSeconds, unit: "s")
                MetresField(label: "Length", value: $options.durationSeconds, unit: "s")
                MetresField(label: "Tail", value: $options.tailSeconds, unit: "s")
            }
            .formStyle(.grouped)
            .disabled(fraction != nil)

            if let fraction {
                ProgressView("Rendering", value: fraction)
            }
            if let result {
                VStack(alignment: .leading, spacing: 2) {
                    Label("Exported \(result.url.lastPathComponent)", systemImage: "checkmark.circle")
                    Text(String(format: "%.1f s, peak %.1f dBFS, %lld clipped samples", result.seconds, decibels(result.summary.peak), result.summary.clippedSamples))
                        .font(.caption)
                        .foregroundStyle(result.summary.clippedSamples > 0 ? .red : .secondary)
                }
                .accessibilityElement(children: .combine)
            }
            if let failure {
                Label(failure, systemImage: "exclamationmark.triangle").foregroundStyle(.red)
            }
            if let limitation = session.engine.status.limitation, !session.engine.status.canRender {
                Text(limitation).font(.callout).foregroundStyle(.secondary)
            }
            HStack {
                Spacer()
                if fraction != nil {
                    Button("Cancel Export", role: .cancel) { task?.cancel() }
                } else {
                    Button("Close", role: .cancel) { dismiss() }
                        .keyboardShortcut(.cancelAction)
                    Button("Export…") { Task { await choose() } }
                        .keyboardShortcut(.defaultAction)
                        .buttonStyle(.glassProminent)
                        .disabled(!session.engine.status.canRender || options.durationSeconds <= 0)
                }
            }
        }
        .padding(20)
        .frame(width: 460)
        .onAppear {
            options.durationSeconds = project.timelineSeconds
            output = project.render.output
        }
    }

    private func choose() async {
        let panel = NSSavePanel()
        panel.allowedContentTypes = [.wav]
        panel.nameFieldStringValue = "Marathon Spatial \(output == .binaural ? "binaural" : "stereo").wav"
        guard panel.runModal() == .OK, let url = panel.url else { return }
        result = nil
        failure = nil
        fraction = 0
        var snapshot = project
        snapshot.render.output = output
        let engine = session.engine
        let options = options
        task = Task {
            do {
                let r = try await Exporter.export(project: snapshot, engine: engine, store: .shared, options: options, to: url) { f in
                    Task { @MainActor in fraction = f }
                }
                result = r
            } catch is CancellationError {
            } catch {
                failure = error.localizedDescription
            }
            fraction = nil
        }
    }
}
