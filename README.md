# Spectral Compressor

Just some experimentation with spectral operations.

## Building

This project builds a VST3 plugin using JUCE. You'll need [CMake 3.15 or higher](https://cliutils.gitlab.io/modern-cmake/chapters/intro/installing.html).

```shell
cmake -Bbuild -DCMAKE_BUILD_TYPE=Release
cmake --build build --parallel
```

TODO: Instructions on obtaining the VST3 SDK, not sure how JUCE's CMake API handles this
TODO: Instructions on what to do with the built plugin.
TODO: Static linking on Linux
