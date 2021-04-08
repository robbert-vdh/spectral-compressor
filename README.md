# Spectral Compressor

[![Automated builds](https://github.com/robbert-vdh/spectral-compressor/workflows/Automated%20builds/badge.svg?branch=master&event=push)](https://github.com/robbert-vdh/spectral-compressor/actions?query=workflow%3A%22Automated+builds%22+branch%3Amaster)

Just some experimentation with spectral operations.

## Download

You can download the latest development binaries for Linux, Windows and macOS
from the [automated
builds](https://github.com/robbert-vdh/spectral-compressor/actions?query=workflow%3A%22Automated+builds%22+branch%3Amaster)
page. GitHub requires you to be signed in to be able to download these files. I
sadly cannot provide any support for the Windows and macOS versions at the
moment, but the plugins should work!

## Building

To build the VST3 plugin, you'll need [CMake 3.15 or
higher](https://cliutils.gitlab.io/modern-cmake/chapters/intro/installing.html)
and a recent C++ compiler.

```shell
cmake -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

You'll find the compiled plugin in `build/SpectralCompressor_artefacts/Release/VST3`.

### Static linking dependencies

By default this project will use the system's copy of FFTW (`fftw3f`) if it's
available, and it will fall back to building a static library and linking to
that. Adding `-DFORCE_STATIC_LINKING=ON` to the command line forces static
linking for distribution. This will also statically linking to the MSVC++
runtime on Windows.
