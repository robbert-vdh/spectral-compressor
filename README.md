# Spectral Compressor

Just some experimentation with spectral operations.

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
