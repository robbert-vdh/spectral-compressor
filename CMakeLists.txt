cmake_minimum_required(VERSION 3.15)

project(spectral_compressor VERSION 0.0.1 LANGUAGES CXX)

# TODO: Figure out a clean way to allow maintainers to disable all static
#       linking and downloading
option(FORCE_STATIC_LINKING "Statically link all dependencies, for distribution" OFF)
# This option should not be necessary as FFTW should at runtime choose the
# correct version based on the current CPU, but apparently on MacOS it doesn't
# do that. This option is there specifically to support a single mid 2012
# MacBook Pro that does not yet support AVX2. It also only does anything unless
# `FORCE_STATIC_LINKING` is also enabled.
option(WITH_FFTW_AVX2 "Enable AVX2 support. By default both AVX and AVX2 are enabled." ON)

# CMake for some reason doesn't enable diagnostic colors by default
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  add_compile_options(-fdiagnostics-color=always)
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
  add_compile_options(-fcolor-diagnostics)
endif()

# Statically link the STL on Windows for the CI builds, and target a lower macOS version
if(FORCE_STATIC_LINKING)
  set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
  set(CMAKE_OSX_DEPLOYMENT_TARGET "10.13")
endif()

#
# Dependencies
#

# Fetch JUCE and other dependencies. For non-JUCE dependencies we'll check
# whether the package is available locally with pkgconfig before trying to
# download and build it ourselves.
include(cmake/CPM.cmake)
# FIXME: JUCE 6.1.3 and 6.1.4 will trigger an X11 error upon loading the plugin
# FIXME: Oh and JUCE 6.1.2 doesn't compile on macOS
if(APPLE)
  CPMAddPackage("gh:juce-framework/JUCE#6.1.3")
else()
  CPMAddPackage("gh:juce-framework/JUCE#6.1.2")
endif()

if(NOT FORCE_STATIC_LINKING)
  find_package(PkgConfig)
  if(PkgConfig_FOUND)
    pkg_search_module(FFTW IMPORTED_TARGET fftw3f)
    pkg_search_module(function2 IMPORTED_TARGET function2)
  endif()
endif()

# Either dynamically link to the system FFTW, or statically link for
# distribution (FFTW also seems to rarely be shipped with static archives)
if(FFTW_FOUND AND NOT FORCE_STATIC_LINKING)
  set(fftw_target PkgConfig::FFTW)
  set(use_shared_fftw TRUE)
else()
  message(STATUS "fftw3f not found using pkgconfig, falling back to a statically linked local build")

  CPMAddPackage(
    NAME FFTW
    VERSIOIN 3.3.9
    URL http://fftw.org/fftw-3.3.9.tar.gz
    URL_HASH SHA256=bf2c7ce40b04ae811af714deb512510cc2c17b9ab9d6ddcf49fe4487eea7af3d
    OPTIONS
    "BUILD_SHARED_LIBS OFF"
    "BUILD_TESTS OFF"
    "ENABLE_FLOAT ON"
    "ENABLE_AVX ON"
    "ENABLE_AVX2 ${WITH_FFTW_AVX2}"
  )

  # CMake builds static libraries without -fPIC by default
  set_target_properties(fftw3f PROPERTIES POSITION_INDEPENDENT_CODE ON)
  set(fftw_target fftw3f)
  set(use_static_fftw TRUE)
endif()

if(NOT function2_FOUND)
  message(STATUS "function2 not found using pkgconfig, downloading from GitHub")
  CPMAddPackage("gh:Naios/function2#4.1.0")
endif()

#
# Plugins
#

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
  src/processor.cpp
  src/utils.cpp)

target_compile_definitions(SpectralCompressor PUBLIC
  JUCE_WEB_BROWSER=0
  JUCE_USE_CURL=0
  JUCE_VST3_CAN_REPLACE_VST2=0
  # We're licensed under the GPL
  JUCE_DISPLAY_SPLASH_SCREEN=0
  $<$<BOOL:${use_shared_fftw}>:JUCE_DSP_USE_SHARED_FFTW=1>
  $<$<BOOL:${use_static_fftw}>:JUCE_DSP_USE_STATIC_FFTW=1>)

target_compile_features(SpectralCompressor PUBLIC cxx_std_20)
set_target_properties(SpectralCompressor PROPERTIES CXX_EXTENSIONS OFF)

# GCC 7+ no longer emits instructions for 128-bit compare-and-swaps and instead
# uses libatomic for this
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(SpectralCompressor PRIVATE -latomic)
endif()

# Statically link the STL on Linux for the CI builds
if(FORCE_STATIC_LINKING AND CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
  target_link_libraries(SpectralCompressor PRIVATE -static-libstdc++)
endif()

target_link_libraries(SpectralCompressor
  PUBLIC
    juce::juce_recommended_warning_flags
    juce::juce_recommended_lto_flags
    juce::juce_recommended_config_flags

  PRIVATE
    juce::juce_audio_utils
    juce::juce_dsp
    ${fftw_target}
    function2)
