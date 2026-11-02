import Foundation

/// The HRTF sets that ship with the library. See hrtf/PROVENANCE.md for each set's origin and licence.
public enum HRTFData {
    public static let setNames = ["sadie2_d1_ku100", "sadie2_d2_kemar", "sadie2_h3", "mit_kemar"]

    public static func url(forSet name: String) -> URL? {
        Bundle.module.url(forResource: name, withExtension: "mshrtf", subdirectory: "hrtf")
    }
}
