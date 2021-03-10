// Spectral Compressor: an FFT based compressor
// Copyright (C) 2021 Robbert van der Helm
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.

#include "processor.h"

#include "editor.h"

// TODO: Rewrite, this is from the example

using juce::uint32;

//==============================================================================
SpectralCompressorProcessor::SpectralCompressorProcessor()
    : AudioProcessor(
          BusesProperties()
#if !JucePlugin_IsMidiEffect
#if !JucePlugin_IsSynth
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
#endif
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)
#endif
              ),
      fft(fft_order) {
}

SpectralCompressorProcessor::~SpectralCompressorProcessor() {}

//==============================================================================
const juce::String SpectralCompressorProcessor::getName() const {
    return JucePlugin_Name;
}

bool SpectralCompressorProcessor::acceptsMidi() const {
#if JucePlugin_WantsMidiInput
    return true;
#else
    return false;
#endif
}

bool SpectralCompressorProcessor::producesMidi() const {
#if JucePlugin_ProducesMidiOutput
    return true;
#else
    return false;
#endif
}

bool SpectralCompressorProcessor::isMidiEffect() const {
#if JucePlugin_IsMidiEffect
    return true;
#else
    return false;
#endif
}

double SpectralCompressorProcessor::getTailLengthSeconds() const {
    return 0.0;
}

int SpectralCompressorProcessor::getNumPrograms() {
    return 1;
}

int SpectralCompressorProcessor::getCurrentProgram() {
    return 0;
}

void SpectralCompressorProcessor::setCurrentProgram(int index) {
    juce::ignoreUnused(index);
}

const juce::String SpectralCompressorProcessor::getProgramName(int /*index*/) {
    return "default";
}

void SpectralCompressorProcessor::changeProgramName(
    int /*index*/,
    const juce::String& /*newName*/) {}

void SpectralCompressorProcessor::prepareToPlay(
    double sampleRate,
    int maximumExpectedSamplesPerBlock) {
    // TODO: The FFT settings are now fixed, we'll want to make this
    //       configurable later

    // JUCE's FFT class interleaves the real and imaginary numbers, so this
    // buffer should be twice the window size in size
    fft_scratch_buffer.resize(fft_window_size * 2);

    // Every FFT bin on both channels gets its own compressor, hooray!
    // TODO: Make the compressor settings configurable
    juce::dsp::Compressor<float> compressor{};
    compressor.setThreshold(-10);
    compressor.setRatio(3.0);
    compressor.setAttack(10);
    compressor.setRelease(30);
    compressor.prepare(juce::dsp::ProcessSpec{
        .sampleRate = sampleRate,
        .maximumBlockSize = static_cast<uint32>(maximumExpectedSamplesPerBlock),
        .numChannels = static_cast<uint32>(getMainBusNumInputChannels())});

    spectral_compressors.resize(static_cast<size_t>(getTotalNumInputChannels()),
                                std::vector(fft_window_size, compressor));
}

void SpectralCompressorProcessor::releaseResources() {
    fft_scratch_buffer.clear();
    spectral_compressors.clear();
}

bool SpectralCompressorProcessor::isBusesLayoutSupported(
    const BusesLayout& layouts) const {
#if JucePlugin_IsMidiEffect
    juce::ignoreUnused(layouts);
    return true;
#else
    // TODO: Why does the example not make whether the input layout is mono or
    //       stereo?
    if (layouts.getMainOutputChannelSet() != juce::AudioChannelSet::mono() &&
        layouts.getMainOutputChannelSet() != juce::AudioChannelSet::stereo()) {
        return false;
    }

#if !JucePlugin_IsSynth
    if (layouts.getMainOutputChannelSet() != layouts.getMainInputChannelSet()) {
        return false;
    }
#endif

    return true;
#endif
}

void SpectralCompressorProcessor::processBlock(
    juce::AudioBuffer<float>& buffer,
    juce::MidiBuffer& /*midiMessages*/) {
    juce::ScopedNoDenormals noDenormals;
    const int input_channels = getTotalNumInputChannels();
    const int output_channels = getTotalNumOutputChannels();

    // Zero out all unused channels
    for (auto i = input_channels; i < output_channels; i++) {
        buffer.clear(i, 0, buffer.getNumSamples());
    }

    for (int channel = 0; channel < input_channels; channel++) {
        float* channelData = buffer.getWritePointer(channel);
        juce::ignoreUnused(channelData);
    }
}

//==============================================================================
bool SpectralCompressorProcessor::hasEditor() const {
    return true;
}

juce::AudioProcessorEditor* SpectralCompressorProcessor::createEditor() {
    return new SpectralCompressorEditor(*this);
}

//==============================================================================
void SpectralCompressorProcessor::getStateInformation(
    juce::MemoryBlock& /*destData*/) {
    // You should use this method to store your parameters in the memory block.
    // You could do that either as raw data, or use the XML or ValueTree classes
    // as intermediaries to make it easy to save and load complex data.
    // TODO: See above
}

void SpectralCompressorProcessor::setStateInformation(const void* /*data*/,
                                                      int /*sizeInBytes*/) {
    // You should use this method to restore your parameters from this memory
    // block, whose contents will have been created by the getStateInformation()
    // call.
    // TODO: Same
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
    return new SpectralCompressorProcessor();
}
