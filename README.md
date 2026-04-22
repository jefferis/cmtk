# cmtk

[![CMake CI](https://github.com/jefferis/cmtk/actions/workflows/cmake-multi-platform.yml/badge.svg?branch=natdev)](https://github.com/jefferis/cmtk/actions/workflows/cmake-multi-platform.yml?query=branch%3Anatdev)
[![Release Builds](https://github.com/jefferis/cmtk/actions/workflows/release-builds.yml/badge.svg?branch=natdev)](https://github.com/jefferis/cmtk/actions/workflows/release-builds.yml?query=branch%3Anatdev)

## README

This README is designed to be viewed on GitHub. See the following project files
for information:

  * [core/README.txt](core/README.txt) Main README
  * [core/CHANGELOG](core/CHANGELOG) Changes for each CMTK version
  * [core/LICENSE](core/LICENSE) See also [core/COPYING.txt](core/COPYING.txt)

## Getting CMTK

Prebuilt binaries are available from several places:

  * [NITRC downloads](https://www.nitrc.org/projects/cmtk/) remain the main upstream download location.
  * [Latest tagged release assets](https://github.com/jefferis/cmtk/releases/latest) are published on GitHub Releases.
  * [Latest `natdev` release assets](https://github.com/jefferis/cmtk/releases/tag/natdev-latest) are published as a rolling prerelease on GitHub Releases.
    Rolling `natdev` package filenames use the form `cmtk-<version>-dev-...` to distinguish them from tagged releases.
  * The [Release Builds workflow](https://github.com/jefferis/cmtk/actions/workflows/release-builds.yml?query=branch%3Anatdev) still keeps per-run workflow artifacts for debugging and reproducibility, but GitHub release assets are the preferred download format.

The current GitHub release-artifact builds target:

  * Linux `x86_64` with OpenMP enabled, currently publishing `.deb`, `.rpm`, and `.tar.gz` packages
  * macOS `arm64` with Grand Central Dispatch enabled, built on `macos-15` and targeting `macOS 12.0+`, currently publishing `.pkg` and `.tar.gz` packages. The `.pkg` is the preferred macOS download.

Current macOS packages are unsigned and not notarized. Gatekeeper will warn that the `.pkg` is from an unidentified developer.

For the current unsigned macOS `.pkg`, use Finder to open it via the context menu:

  1. Control-click or right-click the downloaded `.pkg`
  2. Choose `Open`
  3. Confirm the additional prompt to proceed

If macOS still blocks the installer, open `System Settings` > `Privacy & Security` and use `Open Anyway` for the blocked package.

Source builds can be configured from the repository root with:

```bash
cmake -B build -S core -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Getting Help

If you have any trouble compiling or running CMTK, please write to the 
[CMTK User Forum](https://www.nitrc.org/forum/forum.php?forum_id=857).
