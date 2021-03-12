cmake_minimum_required(VERSION 3.15)

project(spectral_compressor VERSION 0.0.1 LANGUAGES CXX)

include(FetchContent)
FetchContent_Declare(
  JUCE
  GIT_REPOSITORY https://github.com/juce-framework/JUCE.git
  GIT_TAG        6.0.7
)
FetchContent_MakeAvailable(juce)

# We only expose a single VST3 plugin

juce_add_plugin(SpectralCompressor
  PRODUCT_NAME "Spectral Compressor"
  COMPANY_NAME "Robbert van der Helm"
  FORMATS VST3

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
  # We're using the GPL
  JUCE_DISPLAY_SPLASH_SCREEN=0)

# TODO: Fix static linking. With JUCE's default settings we dynamically link to
#       a ton of stuff (and stuff that we don't even use).

target_link_libraries(SpectralCompressor PUBLIC
  juce::juce_recommended_warning_flags
  juce::juce_recommended_lto_flags
  juce::juce_recommended_config_flags

  PRIVATE
  juce::juce_audio_utils
  juce::juce_dsp)
