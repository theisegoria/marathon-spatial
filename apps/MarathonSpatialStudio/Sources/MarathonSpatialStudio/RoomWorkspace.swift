import StudioKit
import SwiftUI

/// Maps world plan coordinates to view points. Z increases toward the bottom of the view.
struct PlanMapping {
    var center: Vec2
    var scale: Double
    var size: CGSize

    init(room: Room, size: CGSize, zoom: Double, pan: CGSize) {
        let xs = room.outline.map(\.x), zs = room.outline.map(\.z)
        let minX = xs.min() ?? -1, maxX = xs.max() ?? 1, minZ = zs.min() ?? -1, maxZ = zs.max() ?? 1
        let fit = min(size.width / max(maxX - minX, 0.1), size.height / max(maxZ - minZ, 0.1)) * 0.82
        scale = fit * zoom
        center = Vec2((minX + maxX) / 2 - pan.width / scale, (minZ + maxZ) / 2 - pan.height / scale)
        self.size = size
    }

    func point(_ p: Vec2) -> CGPoint {
        CGPoint(x: size.width / 2 + (p.x - center.x) * scale, y: size.height / 2 + (p.z - center.z) * scale)
    }
    func point(_ p: Vec3) -> CGPoint { point(p.plan) }
    func world(_ q: CGPoint) -> Vec2 {
        Vec2(center.x + (q.x - size.width / 2) / scale, center.z + (q.y - size.height / 2) / scale)
    }
}

/// Hit target under the pointer.
enum RoomHit: Equatable {
    case listener
    case source(UUID)
    case body(UUID)
    case block(UUID)
    case vertex(Int)
    case wall(Int, Vec2)
    case floor
}

struct RoomWorkspace: View {
    @Binding var project: SpatialProject
    var session: ProjectSession

    /// Gesture in progress. Commits once at the end, so one drag is one undo step.
    @State private var preview: SpatialProject?
    @State private var dragStart: DragStart?
    @State private var blockDraft: (Vec2, Vec2)?
    @State private var pan = CGSize.zero
    @State private var panStart: CGSize?
    @State private var pinchBase: Double?
    @FocusState private var focused: Bool
    @AppStorage("reflectionPathLimit") private var pathLimit = 64

    private struct DragStart {
        var hit: RoomHit?
        var world: Vec2
        var origin: Vec3
    }

    private var shown: SpatialProject { preview ?? project }

    var body: some View {
        VStack(spacing: 0) {
            GeometryReader { geo in
                let map = PlanMapping(room: shown.room, size: geo.size, zoom: session.zoom, pan: pan)
                RoomCanvas(project: shown, session: session, map: map, draft: blockDraft)
                    .contentShape(Rectangle())
                    .gesture(dragGesture(map))
                    .simultaneousGesture(pinch)
                    .focusable()
                    .focusEffectDisabled()
                    .focused($focused)
                    .onKeyPress(phases: [.down, .repeat]) { press in handle(press) }
                    .accessibilityElement(children: .contain)
                    .accessibilityLabel("Room plan")
                    .accessibilityHint("Use the arrow keys to move the selected object. Tab through objects with Command right bracket.")
                    .accessibilityChildren { RoomAccessibility(project: $project, session: session, map: map) }
                    .onTapGesture { focused = true }
            }
            .overlay(alignment: .topLeading) {
                ToolPalette(project: $project, session: session).padding(12)
            }
            .overlay(alignment: .bottomLeading) {
                ReflectionLegend(session: session).padding(12)
            }
            ElevationStrip(project: $project, session: session)
                .frame(height: 96)
                .padding(.horizontal, 12)
                .padding(.bottom, 4)
        }
        .onAppear { focused = true }
        .onChange(of: preview) { _, next in
            if let next { session.projectChanged(next, pathLimit: pathLimit) }
        }
    }

    private var pinch: some Gesture {
        MagnifyGesture()
            .onChanged { value in
                if pinchBase == nil { pinchBase = session.zoom }
                session.zoom = max(0.25, min(8, (pinchBase ?? 1) * value.magnification))
            }
            .onEnded { _ in pinchBase = nil }
    }

    // MARK: Pointer

    private func hit(_ q: CGPoint, map: PlanMapping) -> RoomHit? {
        let p = shown
        let reach = 12.0
        func near(_ w: Vec3) -> Bool { hypot(map.point(w).x - q.x, map.point(w).y - q.y) <= reach }
        if session.tool == .outline {
            for (i, v) in p.room.outline.enumerated() where hypot(map.point(v).x - q.x, map.point(v).y - q.y) <= reach { return .vertex(i) }
        }
        if session.tool != .paint && session.tool != .outline {
            for s in p.sources where near(ProjectCompiler.frame(project: p, time: session.playhead).sources.first { $0.id == s.id }?.position ?? s.position) {
                return .source(s.id)
            }
            for b in p.bodies where near(b.position) { return .body(b.id) }
            if near(p.listener.position) { return .listener }
        }
        let w = map.world(q)
        for b in p.room.blocks where abs(w.x - b.center.x) <= b.size.x / 2 && abs(w.z - b.center.z) <= b.size.z / 2 { return .block(b.id) }
        let n = p.room.outline.count
        var best: (Int, Vec2, Double)?
        for i in 0..<n {
            let c = closestOnSegment(p.room.outline[i], p.room.outline[(i + 1) % n], w).point
            let d = (c - w).length * map.scale
            if d <= 10, d < (best?.2 ?? .infinity) { best = (i, c, d) }
        }
        if let best { return .wall(best.0, best.1) }
        return Polygon2.contains(p.room.outline, w) ? .floor : nil
    }

    private func dragGesture(_ map: PlanMapping) -> some Gesture {
        DragGesture(minimumDistance: 0)
            .onChanged { value in
                focused = true
                let world = map.world(value.location)
                if dragStart == nil { begin(at: value.startLocation, map: map) }
                guard let start = dragStart else { return }
                update(start: start, world: world, moved: hypot(value.translation.width, value.translation.height) > 2,
                       translation: value.translation)
            }
            .onEnded { value in
                finish(at: value.location, map: map, moved: hypot(value.translation.width, value.translation.height) > 2)
            }
    }

    private func begin(at q: CGPoint, map: PlanMapping) {
        let h = hit(q, map: map)
        let world = map.world(q)
        var origin = Vec3(plan: world, y: 0)
        switch h {
        case .listener?:
            session.selection = .listener
            origin = project.listener.position
        case .source(let id)?:
            if session.tool != .strike { session.selection = .source(id) }
            origin = project.sources.first { $0.id == id }?.position ?? origin
        case .body(let id)?:
            if session.tool != .strike { session.selection = .body(id) }
            origin = project.bodies.first { $0.id == id }?.position ?? origin
        case .block(let id)?:
            session.selection = .block(id)
            origin = project.room.blocks.first { $0.id == id }?.center ?? origin
        case .wall(let i, _)?:
            if session.tool == .select { session.selection = .wall(i) }
        case .floor?:
            if session.tool == .select { session.selection = .floor }
        default:
            if session.tool == .select { session.selection = .room }
        }
        dragStart = DragStart(hit: h, world: world, origin: origin)
    }

    private func update(start: DragStart, world: Vec2, moved: Bool, translation: CGSize) {
        switch session.tool {
        case .select, .addSource, .strike:
            guard moved else { return }
            guard let h = start.hit, h != .floor else {
                if session.tool == .select {
                    let base = panStart ?? pan
                    panStart = base
                    pan = CGSize(width: base.width + translation.width, height: base.height + translation.height)
                }
                return
            }
            var next = shown
            let target = Vec2(start.origin.x, start.origin.z) + (world - start.world)
            let snapped = snap(target)
            switch h {
            case .listener: next.listener.position.x = snapped.x; next.listener.position.z = snapped.z
            case .source(let id):
                if let i = next.sources.firstIndex(where: { $0.id == id }) { next.sources[i].position.x = snapped.x; next.sources[i].position.z = snapped.z }
            case .body(let id):
                guard session.tool != .strike else { return }
                if let i = next.bodies.firstIndex(where: { $0.id == id }) { next.bodies[i].position.x = snapped.x; next.bodies[i].position.z = snapped.z }
            case .block(let id):
                if let i = next.room.blocks.firstIndex(where: { $0.id == id }) { next.room.blocks[i].center.x = snapped.x; next.room.blocks[i].center.z = snapped.z }
            default: return
            }
            preview = next
        case .addBlock:
            blockDraft = (start.world, world)
        case .paint:
            paint(start: start, world: world)
        case .outline:
            guard case .vertex(let i)? = start.hit, moved else { return }
            var next = shown
            next.room.outline[i] = snap(world)
            if Polygon2.isSimple(next.room.outline) { preview = next }
        }
    }

    private func paint(start: DragStart, world: Vec2) {
        var next = shown
        let n = next.room.outline.count
        var hitWall = false
        for i in 0..<n {
            let c = closestOnSegment(next.room.outline[i], next.room.outline[(i + 1) % n], world).point
            if (c - world).length < 0.25 {
                next.room.walls[i] = session.brush
                hitWall = true
            }
        }
        if !hitWall, case .block(let id)? = start.hit, let i = next.room.blocks.firstIndex(where: { $0.id == id }) {
            next.room.blocks[i].material = session.brush
        }
        preview = next
    }

    private func finish(at q: CGPoint, map: PlanMapping, moved: Bool) {
        defer {
            dragStart = nil
            blockDraft = nil
            preview = nil
            panStart = nil
        }
        let start = dragStart
        if let preview { project = preview }
        switch session.tool {
        case .addSource where !moved:
            let w = map.world(q)
            guard start?.hit == nil || start?.hit == .floor else { return }
            let s = SoundSource(name: "Source \(project.sources.count + 1)", position: Vec3(plan: snap(w), y: 1.5))
            project.sources.append(s)
            session.selection = .source(s.id)
        case .addBlock:
            guard let (a, b) = blockDraft, moved else { return }
            let size = Vec3(abs(a.x - b.x), min(1.5, project.room.height), abs(a.z - b.z))
            guard size.x > 0.05, size.z > 0.05 else { return }
            let block = Block(name: "Block \(project.room.blocks.count + 1)", center: Vec3((a.x + b.x) / 2, size.y / 2, (a.z + b.z) / 2), size: size, material: session.brush)
            project.room.blocks.append(block)
            session.selection = .block(block.id)
        case .paint where !moved:
            switch start?.hit {
            case .wall(let i, _)?: project.room.walls[i] = session.brush
            case .block(let id)?:
                if let i = project.room.blocks.firstIndex(where: { $0.id == id }) { project.room.blocks[i].material = session.brush }
            case .floor?:
                if NSEvent.modifierFlags.contains(.option) { project.room.ceiling = session.brush } else { project.room.floor = session.brush }
            default: break
            }
        case .outline where !moved:
            if case .vertex(let i)? = start?.hit, NSEvent.modifierFlags.contains(.option) {
                if !project.room.removeVertex(i) { session.alert = "That vertex cannot be removed without crossing walls." }
            } else if case .wall(let i, let c)? = start?.hit {
                project.room.splitWall(i, at: c)
            }
        case .strike where !moved:
            if case .body(let id)? = start?.hit, let body = project.bodies.first(where: { $0.id == id }) {
                let w = map.world(q)
                let u = 0.5 + (w.x - body.position.x) / max(body.sortedSize.x, 0.05)
                let v = 0.5 + (w.z - body.position.z) / max(body.sortedSize.y, 0.05)
                session.strike(body, u: min(1, max(0, u)), v: min(1, max(0, v)), project: &project)
            }
        default: break
        }
    }

    /// Hold Control to place freely; otherwise snap to 5 cm.
    private func snap(_ p: Vec2) -> Vec2 {
        if NSEvent.modifierFlags.contains(.control) { return p }
        return Vec2((p.x * 20).rounded() / 20, (p.z * 20).rounded() / 20)
    }

    // MARK: Keyboard

    private func handle(_ press: KeyPress) -> KeyPress.Result {
        let step = press.modifiers.contains(.shift) ? 1.0 : 0.1
        switch press.key {
        case .leftArrow: return nudge(Vec3(-step, 0, 0))
        case .rightArrow: return nudge(Vec3(step, 0, 0))
        case .upArrow: return nudge(press.modifiers.contains(.option) ? Vec3(0, step, 0) : Vec3(0, 0, -step))
        case .downArrow: return nudge(press.modifiers.contains(.option) ? Vec3(0, -step, 0) : Vec3(0, 0, step))
        case .delete, .deleteForward:
            session.requestDelete = true
            return .handled
        case .space:
            session.togglePlayback()
            return .handled
        default: break
        }
        if let tool = EditTool.allCases.first(where: { String($0.key) == press.characters }) {
            session.tool = tool
            return .handled
        }
        if press.characters == "q" || press.characters == "e" {
            if session.selection == .listener {
                project.listener.yawDegrees += press.characters == "q" ? -15 : 15
                return .handled
            }
        }
        return .ignored
    }

    private func nudge(_ d: Vec3) -> KeyPress.Result {
        switch session.selection {
        case .listener?: project.listener.position += d
        case .source(let id)?:
            guard let i = project.sources.firstIndex(where: { $0.id == id }) else { return .ignored }
            project.sources[i].position += d
        case .body(let id)?:
            guard let i = project.bodies.firstIndex(where: { $0.id == id }) else { return .ignored }
            project.bodies[i].position += d
        case .block(let id)?:
            guard let i = project.room.blocks.firstIndex(where: { $0.id == id }) else { return .ignored }
            project.room.blocks[i].center += d
        default: return .ignored
        }
        return .handled
    }
}

/// Everything drawn in plan.
struct RoomCanvas: View {
    var project: SpatialProject
    var session: ProjectSession
    var map: PlanMapping
    var draft: (Vec2, Vec2)?

    var body: some View {
        Canvas(rendersAsynchronously: false) { context, size in
            drawGrid(&context, size)
            drawRoom(&context)
            if session.showReflections { drawReflections(&context) }
            drawBlocks(&context)
            if let draft { drawDraft(&context, draft) }
            let frame = ProjectCompiler.frame(project: project, time: session.playhead)
            for source in project.sources {
                let p = frame.sources.first { $0.id == source.id }?.position ?? source.position
                drawSource(&context, source, position: p)
            }
            for body in project.bodies { drawBody(&context, body) }
            drawListener(&context, position: frame.listenerPosition)
            if session.tool == .outline {
                for v in project.room.outline {
                    let q = map.point(v)
                    context.fill(Path(ellipseIn: CGRect(x: q.x - 5, y: q.y - 5, width: 10, height: 10)), with: .color(.accentColor))
                }
            }
        }
        .background(Color(nsColor: .underPageBackgroundColor))
    }

    private func material(_ id: MaterialID) -> AcousticMaterial? { project.material(id) }

    private func drawGrid(_ c: inout GraphicsContext, _ size: CGSize) {
        var step = 1.0
        while step * map.scale < 24 { step *= 2 }
        let a = map.world(.zero), b = map.world(CGPoint(x: size.width, y: size.height))
        var grid = Path()
        var x = (a.x / step).rounded(.down) * step
        while x <= b.x {
            grid.move(to: map.point(Vec2(x, a.z))); grid.addLine(to: map.point(Vec2(x, b.z)))
            x += step
        }
        var z = (a.z / step).rounded(.down) * step
        while z <= b.z {
            grid.move(to: map.point(Vec2(a.x, z))); grid.addLine(to: map.point(Vec2(b.x, z)))
            z += step
        }
        c.stroke(grid, with: .color(.secondary.opacity(0.15)), lineWidth: 0.5)
    }

    private func drawRoom(_ c: inout GraphicsContext) {
        let o = project.room.outline
        var floor = Path()
        floor.addLines(o.map { map.point($0) })
        floor.closeSubpath()
        let swatch = material(project.room.floor)?.swatch ?? 0
        c.fill(floor, with: .color(Palette.color(swatch).opacity(0.22)))
        for i in o.indices {
            var wall = Path()
            wall.move(to: map.point(o[i]))
            wall.addLine(to: map.point(o[(i + 1) % o.count]))
            let m = material(project.room.walls[i])?.swatch ?? 0
            let selected = session.selection == .wall(i)
            c.stroke(wall, with: .color(Palette.color(m)), lineWidth: selected ? 9 : 6)
            c.stroke(wall, with: .color(.primary.opacity(0.8)), style: StrokeStyle(lineWidth: 1.5, dash: Palette.dash(m)))
        }
    }

    private func drawReflections(_ c: inout GraphicsContext) {
        for path in session.reflectionPaths where path.order <= session.maximumDisplayOrder {
            var p = Path()
            p.addLines(path.points.map { map.point($0) })
            let opacity = max(0.08, min(0.9, 1 + path.levelDB / 40))
            c.stroke(p, with: .color(Palette.orderColor(path.order).opacity(opacity)), lineWidth: path.order == 1 ? 1.4 : 0.9)
            for q in path.points.dropFirst().dropLast() {
                let s = map.point(q)
                c.fill(Path(ellipseIn: CGRect(x: s.x - 2, y: s.y - 2, width: 4, height: 4)), with: .color(Palette.orderColor(path.order).opacity(opacity)))
            }
        }
    }

    private func drawBlocks(_ c: inout GraphicsContext) {
        for b in project.room.blocks {
            let a = map.point(b.center - b.size / 2), z = map.point(b.center + b.size / 2)
            let r = CGRect(x: min(a.x, z.x), y: min(a.y, z.y), width: abs(a.x - z.x), height: abs(a.y - z.y))
            let m = material(b.material)?.swatch ?? 0
            c.fill(Path(r), with: .color(Palette.color(m).opacity(0.7)))
            c.stroke(Path(r), with: .color(session.selection == .block(b.id) ? .accentColor : .primary.opacity(0.6)),
                     style: StrokeStyle(lineWidth: session.selection == .block(b.id) ? 2.5 : 1, dash: Palette.dash(m)))
        }
    }

    private func drawDraft(_ c: inout GraphicsContext, _ d: (Vec2, Vec2)) {
        let a = map.point(d.0), z = map.point(d.1)
        let r = CGRect(x: min(a.x, z.x), y: min(a.y, z.y), width: abs(a.x - z.x), height: abs(a.y - z.y))
        c.stroke(Path(r), with: .color(.accentColor), style: StrokeStyle(lineWidth: 1.5, dash: [5, 3]))
    }

    private func drawSource(_ c: inout GraphicsContext, _ s: SoundSource, position: Vec3) {
        let q = map.point(position)
        let selected = session.selection == .source(s.id)
        let tint = s.muted ? Color.secondary : Palette.sourceTint
        if case let .orbit(center, radius, _) = s.motion {
            let k = radius * map.scale, m = map.point(center)
            c.stroke(Path(ellipseIn: CGRect(x: m.x - k, y: m.y - k, width: 2 * k, height: 2 * k)), with: .color(tint.opacity(0.5)), style: StrokeStyle(lineWidth: 1, dash: [4, 4]))
        } else if case let .shuttle(offset, _) = s.motion {
            var p = Path()
            p.move(to: map.point(s.position)); p.addLine(to: map.point(s.position + offset))
            c.stroke(p, with: .color(tint.opacity(0.5)), style: StrokeStyle(lineWidth: 1, dash: [4, 4]))
        }
        if s.directivity.dipoleWeight > 0 {
            var lobe = Path()
            for a in stride(from: 0.0, through: 360.0, by: 6.0) {
                let g = s.directivity.gain(offAxisDegrees: a) * 34
                let t = (a + s.directivity.yawDegrees) * .pi / 180
                let point = CGPoint(x: q.x + sin(t) * g, y: q.y - cos(t) * g)
                if a == 0 { lobe.move(to: point) } else { lobe.addLine(to: point) }
            }
            lobe.closeSubpath()
            c.fill(lobe, with: .color(tint.opacity(0.16)))
        }
        let level = decibels(session.meters.sources[s.id] ?? 0)
        let glow = 11 + 10 * max(0, min(1, (level + 60) / 60))
        c.fill(Path(ellipseIn: CGRect(x: q.x - glow, y: q.y - glow, width: 2 * glow, height: 2 * glow)), with: .color(tint.opacity(0.22)))
        c.fill(Path(ellipseIn: CGRect(x: q.x - 8, y: q.y - 8, width: 16, height: 16)), with: .color(tint))
        if selected { c.stroke(Path(ellipseIn: CGRect(x: q.x - 12, y: q.y - 12, width: 24, height: 24)), with: .color(.accentColor), lineWidth: 2) }
        c.draw(Text(s.name).font(.caption), at: CGPoint(x: q.x, y: q.y + 20))
    }

    private func drawBody(_ c: inout GraphicsContext, _ b: SoundingBody) {
        let q = map.point(b.position)
        var d = Path()
        d.addLines([CGPoint(x: q.x, y: q.y - 10), CGPoint(x: q.x + 10, y: q.y), CGPoint(x: q.x, y: q.y + 10), CGPoint(x: q.x - 10, y: q.y)])
        d.closeSubpath()
        c.fill(d, with: .color(Palette.bodyTint))
        if session.selection == .body(b.id) { c.stroke(d, with: .color(.accentColor), lineWidth: 2) }
        c.draw(Text(b.name).font(.caption), at: CGPoint(x: q.x, y: q.y + 20))
    }

    private func drawListener(_ c: inout GraphicsContext, position: Vec3) {
        let q = map.point(position)
        let f = project.listener.forward
        let len = max((f.x * f.x + f.z * f.z).squareRoot(), 1e-6)
        let dir = CGPoint(x: f.x / len, y: f.z / len)
        var arrow = Path()
        arrow.move(to: CGPoint(x: q.x + dir.x * 24, y: q.y + dir.y * 24))
        arrow.addLine(to: CGPoint(x: q.x - dir.y * 8, y: q.y + dir.x * 8))
        arrow.addLine(to: CGPoint(x: q.x + dir.y * 8, y: q.y - dir.x * 8))
        arrow.closeSubpath()
        c.fill(arrow, with: .color(Palette.listenerTint.opacity(0.8)))
        c.fill(Path(ellipseIn: CGRect(x: q.x - 9, y: q.y - 9, width: 18, height: 18)), with: .color(Palette.listenerTint))
        if session.selection == .listener {
            c.stroke(Path(ellipseIn: CGRect(x: q.x - 13, y: q.y - 13, width: 26, height: 26)), with: .color(.primary), lineWidth: 2)
        }
    }
}

/// One accessibility element per object: VoiceOver can find, read and move each without the canvas.
struct RoomAccessibility: View {
    @Binding var project: SpatialProject
    var session: ProjectSession
    var map: PlanMapping

    var body: some View {
        element("Listener", project.listener.position, kind: "listener", .listener) { project.listener.position += $0 }
        ForEach(project.sources) { s in
            element(s.name, s.position, kind: "sound source", .source(s.id)) { d in
                if let i = project.sources.firstIndex(where: { $0.id == s.id }) { project.sources[i].position += d }
            }
        }
        ForEach(project.bodies) { b in
            element(b.name, b.position, kind: "sounding body", .body(b.id)) { d in
                if let i = project.bodies.firstIndex(where: { $0.id == b.id }) { project.bodies[i].position += d }
            }
        }
    }

    private func element(_ name: String, _ p: Vec3, kind: String, _ selection: SceneSelection, move: @escaping (Vec3) -> Void) -> some View {
        let q = map.point(p)
        return Rectangle()
            .frame(width: 24, height: 24)
            .position(q)
            .accessibilityLabel("\(name), \(kind)")
            .accessibilityValue(String(format: "x %.1f, height %.1f, z %.1f metres", p.x, p.y, p.z))
            .accessibilityAddTraits(session.selection == selection ? [.isSelected, .isButton] : .isButton)
            .accessibilityAction { session.selection = selection }
            .accessibilityAction(named: "Move left") { move(Vec3(-0.5, 0, 0)) }
            .accessibilityAction(named: "Move right") { move(Vec3(0.5, 0, 0)) }
            .accessibilityAction(named: "Move forward") { move(Vec3(0, 0, -0.5)) }
            .accessibilityAction(named: "Move back") { move(Vec3(0, 0, 0.5)) }
            .accessibilityAction(named: "Move up") { move(Vec3(0, 0.25, 0)) }
            .accessibilityAction(named: "Move down") { move(Vec3(0, -0.25, 0)) }
    }
}

struct ToolPalette: View {
    @Binding var project: SpatialProject
    var session: ProjectSession

    var body: some View {
        @Bindable var session = session
        GlassEffectContainer(spacing: 8) {
            VStack(alignment: .leading, spacing: 8) {
                HStack(spacing: 2) {
                    ForEach(EditTool.allCases) { tool in
                        Button { session.tool = tool } label: {
                            Label(tool.title, systemImage: tool.symbol).labelStyle(.iconOnly).frame(width: 26, height: 22)
                        }
                            .buttonStyle(.glass)
                            .tint(session.tool == tool ? .accentColor : nil)
                            .help("\(tool.title) (\(String(tool.key).uppercased()))")
                            .accessibilityAddTraits(session.tool == tool ? .isSelected : [])
                    }
                }
                if session.tool == .paint || session.tool == .addBlock {
                    Picker("Material", selection: $session.brush) {
                        ForEach(project.materials) { m in
                            Label { Text(m.name) } icon: { MaterialSwatch(material: m) }.tag(m.id)
                        }
                    }
                    .frame(maxWidth: 240)
                    .padding(8)
                    .glassPanel(cornerRadius: 12)
                    if session.tool == .paint {
                        Text("Click or drag along walls. Click the floor to paint it; Option-click paints the ceiling.")
                            .font(.caption)
                            .frame(maxWidth: 240, alignment: .leading)
                            .padding(8)
                            .glassPanel(cornerRadius: 12)
                    }
                }
                if session.tool == .outline {
                    Text("Drag a corner to move it. Click a wall to add a corner. Option-click a corner to remove it.")
                        .font(.caption)
                        .frame(maxWidth: 240, alignment: .leading)
                        .padding(8)
                        .glassPanel(cornerRadius: 12)
                }
            }
        }
    }
}

struct ReflectionLegend: View {
    var session: ProjectSession

    var body: some View {
        @Bindable var session = session
        HStack(spacing: 10) {
            Toggle("Reflections", isOn: $session.showReflections)
            if session.showReflections {
                Stepper("Order \(session.maximumDisplayOrder)", value: $session.maximumDisplayOrder, in: 1...3)
                ForEach(1...session.maximumDisplayOrder, id: \.self) { o in
                    Rectangle().fill(Palette.orderColor(o)).frame(width: 14, height: 3).accessibilityHidden(true)
                }
                Text("\(session.reflectionPaths.count) paths shown, up to 12 per source sent to the library")
                    .font(.caption)
                    .foregroundStyle(.secondary)
            }
        }
        .controlSize(.small)
        .padding(8)
        .glassPanel(cornerRadius: 12)
    }
}

/// Side view: height of everything, in X. Drag to change height.
struct ElevationStrip: View {
    @Binding var project: SpatialProject
    var session: ProjectSession
    @State private var dragging: SceneSelection?
    /// Height under the pointer. Commits at the end of the drag, so a drag is one undo step.
    @State private var draftHeight: Double?

    private struct Item: Identifiable {
        var selection: SceneSelection
        var name: String
        var position: Vec3
        var tint: Color
        var id: SceneSelection { selection }
    }

    private var items: [Item] {
        [Item(selection: .listener, name: "Listener", position: project.listener.position, tint: Palette.listenerTint)]
            + project.sources.map { Item(selection: .source($0.id), name: $0.name, position: $0.position, tint: Palette.sourceTint) }
            + project.bodies.map { Item(selection: .body($0.id), name: $0.name, position: $0.position, tint: Palette.bodyTint) }
    }

    private func height(of item: Item) -> Double {
        dragging == item.selection ? draftHeight ?? item.position.y : item.position.y
    }

    private func bounds() -> (Double, Double) {
        let xs = project.room.outline.map(\.x)
        return (xs.min() ?? -1, xs.max() ?? 1)
    }

    var body: some View {
        GeometryReader { geo in
            let (lo, hi) = bounds()
            let w = max(hi - lo, 0.1), h = max(project.room.height, 0.1), inset = 10.0
            let x = { (v: Double) in inset + (v - lo) / w * (geo.size.width - 2 * inset) }
            let y = { (v: Double) in geo.size.height - inset - v / h * (geo.size.height - 2 * inset) }
            ZStack(alignment: .topLeading) {
                Canvas { c, _ in
                    let r = CGRect(x: x(lo), y: y(h), width: x(hi) - x(lo), height: y(0) - y(h))
                    c.stroke(Path(r), with: .color(.secondary.opacity(0.7)), lineWidth: 1)
                    for b in project.room.blocks {
                        let a = CGRect(x: x(b.center.x - b.size.x / 2), y: y(b.center.y + b.size.y / 2),
                                       width: b.size.x / w * (geo.size.width - 2 * inset), height: b.size.y / h * (geo.size.height - 2 * inset))
                        c.fill(Path(a), with: .color(.secondary.opacity(0.35)))
                    }
                    for item in items {
                        let q = CGPoint(x: x(item.position.x), y: y(height(of: item)))
                        c.fill(Path(ellipseIn: CGRect(x: q.x - 6, y: q.y - 6, width: 12, height: 12)), with: .color(item.tint))
                        if session.selection == item.selection {
                            c.stroke(Path(ellipseIn: CGRect(x: q.x - 9, y: q.y - 9, width: 18, height: 18)), with: .color(.accentColor), lineWidth: 2)
                        }
                    }
                }
                .contentShape(Rectangle())
                .gesture(DragGesture(minimumDistance: 0).onChanged { value in
                    if dragging == nil {
                        dragging = items.min {
                            hypot(x($0.position.x) - value.startLocation.x, y(height(of: $0)) - value.startLocation.y)
                                < hypot(x($1.position.x) - value.startLocation.x, y(height(of: $1)) - value.startLocation.y)
                        }?.selection
                        if let dragging { session.selection = dragging }
                    }
                    let v = (geo.size.height - inset - value.location.y) / (geo.size.height - 2 * inset) * h
                    draftHeight = max(0, min(h, (v * 20).rounded() / 20))
                }.onEnded { _ in
                    if let draftHeight { setHeight(draftHeight) }
                    draftHeight = nil
                    dragging = nil
                })
                Text("Height, side view").font(.caption2).foregroundStyle(.secondary).padding(.leading, inset + 6).padding(.top, inset + 3)
            }
            .accessibilityElement()
            .accessibilityLabel("Height view")
            .accessibilityHint("Heights can also be set in the inspector, or with Option and the arrow keys.")
        }
    }

    private func setHeight(_ v: Double) {
        switch dragging {
        case .listener?: project.listener.position.y = v
        case .source(let id)?:
            if let i = project.sources.firstIndex(where: { $0.id == id }) { project.sources[i].position.y = v }
        case .body(let id)?:
            if let i = project.bodies.firstIndex(where: { $0.id == id }) { project.bodies[i].position.y = v }
        default: break
        }
    }
}
