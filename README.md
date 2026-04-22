# cmtk

[![CMake CI](https://github.com/jefferis/cmtk/actions/workflows/cmake-multi-platform.yml/badge.svg?branch=natdev)](https://github.com/jefferis/cmtk/actions/workflows/cmake-multi-platform.yml?query=branch%3Anatdev)

## README

This README is designed to be viewed on GitHub. See the following project files
for information:

  * [core/README.txt](core/README.txt) Main README
  * [core/CHANGELOG](core/CHANGELOG) Changes for each CMTK version
  * [core/LICENSE](core/LICENSE) See also [core/COPYING.txt](core/COPYING.txt)

## Getting CMTK

Prebuilt binaries are available from several places:

  * [NITRC downloads](https://www.nitrc.org/projects/cmtk/) remain the main upstream download location.
  * [GitHub Releases](https://github.com/jefferis/cmtk/releases) are the intended home for tagged release artifacts in this mirror.
  * The latest `natdev` build artifacts are published by the [Release Builds workflow](https://github.com/jefferis/cmtk/actions/workflows/release-builds.yml?query=branch%3Anatdev). Open the most recent successful `natdev` run and download the workflow artifacts:
    * `cmtk-linux-x86_64-openmp`
    * `cmtk-macos-arm64-gcd`

The current GitHub release-artifact builds target:

  * Linux `x86_64` with OpenMP enabled
  * macOS `arm64` with Grand Central Dispatch enabled, built on `macos-15` and targeting `macOS 12.0+`

Source builds can be configured from the repository root with:

```bash
cmake -B build -S core -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build --output-on-failure
```

## Getting Help

If you have any trouble compiling or running CMTK, please write to the 
[CMTK User Forum](https://www.nitrc.org/forum/forum.php?forum_id=857).
