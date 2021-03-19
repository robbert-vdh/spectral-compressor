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

# Fetch JUCE and other dependencies

include(cmake/CPM.cmake)

CPMAddPackage("gh:juce-framework/JUCE#6.0.7")
# TODO: How can we fall back to linking to system FFTW? This seems to be quite a
#       mess, and the CMake config for FFTW that ships with most distros is
#       invalid
# TODO: Is there no way to add a checksum or something?
CPMAddPackage(
  NAME fftw
  VERSIOIN 3.3.9
  URL http://fftw.org/fftw-3.3.9.tar.gz
  OPTIONS
    "BUILD_TESTS OFF"
    "ENABLE_FLOAT ON"
    "ENABLE_SSE2 ON"
    "ENABLE_AVX ON"
    "ENABLE_AVX2 ON"
)

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
  JUCE_DSP_USE_STATIC_FFTW=1)

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
    fftw3f)
