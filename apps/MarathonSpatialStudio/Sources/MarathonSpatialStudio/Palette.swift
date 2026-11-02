import StudioKit
import SwiftUI

/// Material colours, each paired with a pattern so that no material depends on colour alone: a dash on the plan, a drawn surface on a swatch.
enum Palette {
    static let colors: [Color] = [
        Color(red: 0.55, green: 0.57, blue: 0.60), Color(red: 0.70, green: 0.36, blue: 0.27), Color(red: 0.88, green: 0.85, blue: 0.78),
        Color(red: 0.77, green: 0.75, blue: 0.70), Color(red: 0.45, green: 0.75, blue: 0.85), Color(red: 0.62, green: 0.42, blue: 0.22),
        Color(red: 0.76, green: 0.58, blue: 0.36), Color(red: 0.55, green: 0.30, blue: 0.55), Color(red: 0.70, green: 0.20, blue: 0.30),
        Color(red: 0.92, green: 0.92, blue: 0.95), Color(red: 0.40, green: 0.45, blue: 0.52), Color(red: 0.35, green: 0.60, blue: 0.55),
        Color(red: 0.20, green: 0.45, blue: 0.75),
    ]
    static let dashes: [[CGFloat]] = [[], [10, 3], [3, 3], [14, 3, 3, 3], [6, 6], [2, 5], [10, 3, 3, 3, 3, 3], [4, 2], [16, 4], [8, 2, 2, 2], [12, 6], [5, 3, 1, 3], [1, 3]]

    static func color(_ swatch: Int) -> Color { colors[abs(swatch) % colors.count] }
    static func dash(_ swatch: Int) -> [CGFloat] { dashes[abs(swatch) % dashes.count] }

    /// The surface of a material, in vector: the order follows `colors`.
    static func surface(_ swatch: Int, in s: CGSize) -> Path {
        let (w, h) = (s.width, s.height)
        var p = Path()
        func line(_ x0: CGFloat, _ y0: CGFloat, _ x1: CGFloat, _ y1: CGFloat) {
            p.move(to: CGPoint(x: x0 * w, y: y0 * h))
            p.addLine(to: CGPoint(x: x1 * w, y: y1 * h))
        }
        func dot(_ x: CGFloat, _ y: CGFloat, _ radius: CGFloat = 0.4) {
            p.addEllipse(in: CGRect(x: x * w - radius, y: y * h - radius, width: 2 * radius, height: 2 * radius))
        }
        func wave(_ y: CGFloat, _ depth: CGFloat, _ crests: Int) {
            p.move(to: CGPoint(x: 0, y: y * h))
            for i in 0..<crests {
                let (a, b) = (CGFloat(i) / CGFloat(crests) * w, CGFloat(i + 1) / CGFloat(crests) * w)
                p.addQuadCurve(to: CGPoint(x: b, y: y * h), control: CGPoint(x: (a + b) / 2, y: (y + (i % 2 == 0 ? -depth : depth)) * h))
            }
        }
        switch abs(swatch) % colors.count {
        case 0: // Painted concrete: aggregate.
            for (x, y) in [(0.14, 0.3), (0.3, 0.72), (0.43, 0.26), (0.56, 0.6), (0.7, 0.3), (0.84, 0.74), (0.9, 0.22), (0.2, 0.52)] { dot(x, y) }
        case 1: // Bare brick: running bond.
            for y in [1.0 / 3, 2.0 / 3] { line(0, y, 1, y) }
            for x in [0.25, 0.75] { line(x, 0, x, 1.0 / 3); line(x, 2.0 / 3, x, 1) }
            line(0.5, 1.0 / 3, 0.5, 2.0 / 3)
        case 2: // Plaster: trowel sweeps.
            wave(0.34, 0.2, 2); wave(0.72, 0.2, 3)
        case 3: // Gypsum board: a joint and its screws.
            line(0.5, 0, 0.5, 1)
            for x in [0.38, 0.62] { for y in [0.25, 0.75] { dot(x, y) } }
        case 4: // Window glass: glints.
            line(0.2, 0.85, 0.5, 0.15); line(0.36, 0.85, 0.66, 0.15); line(0.7, 0.85, 0.82, 0.57)
        case 5: // Wood panelling: upright boards.
            for x in [0.25, 0.5, 0.75] { line(x, 0, x, 1) }
            wave(0.3, 0.08, 1)
        case 6: // Wood floor: staggered planks.
            for y in [0.25, 0.5, 0.75] { line(0, y, 1, y) }
            line(0.3, 0, 0.3, 0.25); line(0.72, 0.25, 0.72, 0.5); line(0.18, 0.5, 0.18, 0.75); line(0.58, 0.75, 0.58, 1)
        case 7: // Heavy carpet: tufts.
            for row in 0..<3 {
                for i in 0..<6 {
                    let x = (CGFloat(i) + (row % 2 == 0 ? 0.3 : 0.8)) / 6, y = 0.22 + CGFloat(row) * 0.28
                    line(x, y + 0.1, x + 0.03, y - 0.1)
                }
            }
        case 8: // Heavy curtain: folds.
            for x in [0.2, 0.4, 0.6, 0.8] as [CGFloat] {
                p.move(to: CGPoint(x: x * w, y: 0))
                p.addCurve(to: CGPoint(x: x * w, y: h), control1: CGPoint(x: (x + 0.07) * w, y: 0.33 * h), control2: CGPoint(x: (x - 0.07) * w, y: 0.66 * h))
            }
        case 9: // Acoustic ceiling tile: a grid, perforated.
            for x in [1.0 / 3, 2.0 / 3] { line(x, 0, x, 1) }
            line(0, 0.5, 1, 0.5)
            for x in [1.0 / 6, 0.5, 5.0 / 6] { for y in [0.25, 0.75] { dot(x, y, 0.5) } }
        case 10: // Steel plate: rivets inside the edge.
            p.addRoundedRect(in: CGRect(x: 2.5, y: 2.5, width: w - 5, height: h - 5), cornerSize: CGSize(width: 1.5, height: 1.5))
            for x in [0.22, 0.78] { for y in [0.34, 0.66] { dot(x, y, 0.6) } }
        case 11: // Polished stone: veins.
            p.move(to: CGPoint(x: 0, y: 0.7 * h))
            p.addLine(to: CGPoint(x: 0.3 * w, y: 0.5 * h))
            p.addLine(to: CGPoint(x: 0.45 * w, y: 0.62 * h))
            p.addLine(to: CGPoint(x: 0.75 * w, y: 0.2 * h))
            p.addLine(to: CGPoint(x: w, y: 0.3 * h))
            line(0.45, 0.62, 0.62, 0.95)
        default: // Water surface.
            wave(0.3, 0.14, 4); wave(0.66, 0.14, 4)
        }
        return p
    }

    static let sourceTint = Color.orange
    static let listenerTint = Color.accentColor
    static let bodyTint = Color.mint

    /// Early paths warm, later ones cool.
    static func orderColor(_ order: Int) -> Color {
        switch order {
        case 1: .yellow
        case 2: .orange
        case 3: .pink
        default: .purple
        }
    }
}

struct MaterialSwatch: View {
    var material: AcousticMaterial
    var size: CGFloat = 16

    var body: some View {
        Canvas { context, s in
            let tile = Path(roundedRect: CGRect(origin: .zero, size: s).insetBy(dx: 0.25, dy: 0.25), cornerRadius: 4)
            context.fill(tile, with: .color(Palette.color(material.swatch)))
            var surface = context
            surface.clip(to: tile)
            // On the tile's own colour, so it reads alike in light mode, dark mode and a selected row.
            surface.stroke(Palette.surface(material.swatch, in: s), with: .color(.black.opacity(0.5)), style: StrokeStyle(lineWidth: 0.75, lineCap: .round, lineJoin: .round))
            context.stroke(tile, with: .color(.black.opacity(0.2)), lineWidth: 0.5)
        }
        .frame(width: size * 1.6, height: size)
        .accessibilityHidden(true)
    }
}

extension View {
    /// A floating Liquid Glass panel.
    func glassPanel(cornerRadius: CGFloat = 16) -> some View {
        glassEffect(.regular, in: .rect(cornerRadius: cornerRadius))
    }
}

/// A number field with a unit, for positions and sizes.
struct MetresField: View {
    var label: String
    @Binding var value: Double
    var unit = "m"

    var body: some View {
        LabeledContent(label) {
            HStack(spacing: 4) {
                TextField(label, value: $value, format: .number.precision(.fractionLength(0...3)))
                    .labelsHidden()
                    .multilineTextAlignment(.trailing)
                    .frame(maxWidth: 80)
                Text(unit).foregroundStyle(.secondary)
            }
        }
    }
}

struct Vec3Fields: View {
    var labels = ["X", "Y (height)", "Z"]
    @Binding var value: Vec3

    var body: some View {
        MetresField(label: labels[0], value: $value.x)
        MetresField(label: labels[1], value: $value.y)
        MetresField(label: labels[2], value: $value.z)
    }
}

func decibels(_ linear: Float) -> Double { linear > 1e-6 ? 20 * log10(Double(linear)) : -120 }
