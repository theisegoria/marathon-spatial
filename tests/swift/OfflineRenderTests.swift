import Foundation
import Testing

import MarathonSpatial
import MarathonSpatialHRTFData

// These render into memory. SpatialOutput is never created here.
@Suite struct OfflineRender {
    @Test func aSourceToTheRightIsLouderInTheRightEar() throws {
        let url = try #require(HRTFData.url(forSet: "sadie2_d1_ku100"))
        let engine = try SpatialEngine(.init(limiter: false))
        let set = try HRTF(contentsOf: url, sampleRate: engine.sampleRate)
        #expect(set.info.directions == 8802)
        try engine.setHRTF(set)

        let tone = (0..<4800).map { 0.5 * sin(Float($0) * 2 * .pi * 3000 / 48000) }
        let buffer = try engine.makeBuffer(interleaved: tone, channels: 1, sampleRate: 48000)
        let source = try engine.makeSource()
        var parameters = SpatialEngine.defaultSourceParameters
        parameters.position = .init(x: 1, y: 0, z: 0)
        try engine.set(source, parameters: parameters)
        try engine.attach(buffer, to: source)
        try engine.play(source)
        try engine.commit()

        let stereo = engine.renderInterleaved(frames: 4800)
        var left: Float = 0, right: Float = 0
        for i in stride(from: 0, to: stereo.count, by: 2) {
            left += stereo[i] * stereo[i]
            right += stereo[i + 1] * stereo[i + 1]
        }
        #expect(right > 4 * left)
        #expect(engine.stats.faults == 0)
    }
}
