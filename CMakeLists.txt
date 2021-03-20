cmake_minimum_required(VERSION 3.15)

project(spectral_compressor VERSION 0.0.1 LANGUAGES CXX)

# During development you'll probably want to add
# `-DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DJUCE_ENABLE_MODULE_SOURCE_GROUPS=1`
# to the CMake command line
set_property(GLOBAL PROPERTY USE_FOLDERS YES)

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-fdiagnostics-color=always)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_compile_options(-fcolor-diagnostics)
endif()

# Fetch JUCE and other dependencies. For non-JUCE dependencies we'll check
# whether the package is available locally with pkgconfig before trying to
# download and build it ourselves.

find_package(PkgConfig)
if(PkgConfig_FOUND)
  pkg_search_module(FFTW fftw3f IMPORTED_TARGET)
endif()

include(cmake/CPM.cmake)
CPMAddPackage("gh:juce-framework/JUCE#6.0.7")

# TODO: Add a way to force statically linking for distribution
if(FFTW_FOUND)
  set(fftw_target PkgConfig::FFTW)
  set(use_shared_fftw TRUE)
else()
  message(STATUS "fftw3f not found using pkgconfig, falling back to a statically linked local build")

  # TODO: Is there no way to add a checksum or something?
  CPMAddPackage(
    NAME FFTW
    VERSIOIN 3.3.9
    URL http://fftw.org/fftw-3.3.9.tar.gz
    OPTIONS
    "BUILD_SHARED_LIBS OFF"
    "BUILD_TESTS OFF"
    "ENABLE_FLOAT ON"
    "ENABLE_SSE2 ON"
    "ENABLE_AVX2 ON"
  )

  # CMake builds static libraries without -fPIC by default
  set_target_properties(fftw3f PROPERTIES POSITION_INDEPENDENT_CODE ON)
  set(fftw_target fftw3f)
  set(use_static_fftw TRUE)
endif()

# We only expose a single VST3 plugin

juce_add_plugin(SpectralCompressor
  PRODUCT_NAME "Spectral Compressor"
  COMPANY_NAME "Robbert van der Helm"
  FORMATS VST3 AU

  PLUGIN_MANUFACTURER_CODE RvdH
  PLUGIN_CODE Spcc

  IS_SYNTH FALSE
  NEEDS_MIDI_INPUT FALSE
  NEEDS_MIDI_OUTPUT FALSE
  IS_MIDI_EFFECT FALSE
  EDITOR_WANTS_KEYBOARD_FOCUS FALSE

  VST3_CATEGORIES Fx Dynamics)

target_sources(SpectralCompressor PRIVATE
  src/editor.cpp
  src/processor.cpp)

target_compile_features(SpectralCompressor PUBLIC cxx_std_20)
set_target_properties(SpectralCompressor PROPERTIES CXX_EXTENSIONS OFF)

target_compile_definitions(SpectralCompressor PUBLIC
  JUCE_WEB_BROWSER=0
  JUCE_USE_CURL=0
  JUCE_VST3_CAN_REPLACE_VST2=0
  # We're licensed under the GPL
  JUCE_DISPLAY_SPLASH_SCREEN=0
  $<$<BOOL:${use_shared_fftw}>:JUCE_DSP_USE_SHARED_FFTW=1>
  $<$<BOOL:${use_static_fftw}>:JUCE_DSP_USE_STATIC_FFTW=1>)

# TODO: Fix static linking. With JUCE's default settings we dynamically link to
#       a ton of stuff (and stuff that we don't even use).
target_link_libraries(SpectralCompressor
  PUBLIC
    juce::juce_recommended_warning_flags
    juce::juce_recommended_lto_flags
    juce::juce_recommended_config_flags

  PRIVATE
    juce::juce_audio_utils
    juce::juce_dsp
    ${fftw_target})
