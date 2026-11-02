#!/bin/sh
# Assemble "build/Marathon Spatial Studio.app" with an ad hoc signature. It builds the app and never launches it.
# The document type lives in Info.plist, so documents save and open only from this bundle, not from `swift run`.
set -eu
cd "$(dirname "$0")/.."
export MARATHON_SPATIAL_REQUIRE_LINKED=1
swift package --manifest-cache none describe >/dev/null
swift build -c "${1:-release}" --product MarathonSpatialStudio
BIN="$(swift build -c "${1:-release}" --show-bin-path)/MarathonSpatialStudio"
APP="build/Marathon Spatial Studio.app"
rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS"
cp "$BIN" "$APP/Contents/MacOS/MarathonSpatialStudio"
cp Support/Info.plist "$APP/Contents/Info.plist"
mkdir -p "$APP/Contents/Resources/HRTF"
# Icon Composer's layered source: Assets.car for Liquid Glass, AppIcon.icns for everything else.
xcrun actool --compile "$APP/Contents/Resources" --platform macosx --minimum-deployment-target 27.0 --app-icon AppIcon \
    --output-partial-info-plist build/AppIcon-partial.plist Support/AppIcon.icon >/dev/null
cp ../../data/hrtf/*.mshrtf "$APP/Contents/Resources/HRTF/" 2>/dev/null || echo "no HRTF packs in data/hrtf; the app will offer the analytic sphere"
# What the licences require any copy of the app to carry.
cp ../../LICENSE ../../NOTICE ../../THIRD_PARTY_NOTICES.md "$APP/Contents/Resources/"
cp ../../third_party/libmysofa/LICENSE "$APP/Contents/Resources/LICENSE.libmysofa"
cp -R ../../LICENSES "$APP/Contents/Resources/LICENSES"
METALLIB="$(find "${MARATHON_SPATIAL_LIB_DIR:-../../build}" -name acoustic_trace.metallib -print -quit 2>/dev/null || true)"
if [ -n "$METALLIB" ]; then
    cp "$METALLIB" "$APP/Contents/Resources/acoustic_trace.metallib"
else
    echo "no acoustic_trace.metallib found; the app will trace reflections on the CPU"
fi
plutil -lint "$APP/Contents/Info.plist"
codesign --force --sign - "$APP"
echo "$APP"
