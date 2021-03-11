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
    setLatencySamples(4096);
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
    fft_scratch_buffer.resize(static_cast<size_t>(getTotalNumInputChannels()),
                              std::vector(fft_window_size * 2, 0.0f));

    // Every FFT bin on both channels gets its own compressor, hooray!
    // TODO: Make the compressor settings configurable
    // TODO: These settings are also very extreme
    juce::dsp::Compressor<float> compressor{};
    compressor.setThreshold(-20);
    compressor.setRatio(50.0);
    compressor.setAttack(10);
    compressor.setRelease(100);
    compressor.prepare(juce::dsp::ProcessSpec{
        .sampleRate = sampleRate,
        .maximumBlockSize = static_cast<uint32>(maximumExpectedSamplesPerBlock),
        .numChannels = static_cast<uint32>(getMainBusNumInputChannels())});

    spectral_compressors.resize(fft_window_size, compressor);

    first_iteration = true;

    // We use ring buffers to fill our FFT buffers
    ring_buffers.resize(static_cast<size_t>(getTotalNumInputChannels()),
                        std::vector(fft_window_size, 0.0f));
    ring_buffer_pos.resize(static_cast<size_t>(getTotalNumInputChannels()));
}

void SpectralCompressorProcessor::releaseResources() {
    fft_scratch_buffer.clear();
    spectral_compressors.clear();
    ring_buffers.clear();
    ring_buffer_pos.clear();
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

    const size_t input_channels =
        static_cast<size_t>(getTotalNumInputChannels());
    const size_t output_channels =
        static_cast<size_t>(getTotalNumOutputChannels());
    const size_t num_samples = static_cast<size_t>(buffer.getNumSamples());

    // Zero out all unused channels
    for (auto i = input_channels; i < output_channels; i++) {
        buffer.clear(i, 0, buffer.getNumSamples());
    }

    // TODO: Add oversampling, potentially reduce latency
    // TODO: Handle buffers that are smaller than our FFT window
    // TODO: Handle buffers that are not a clean multiple of our FFT window
    // TODO: First iteration should be skipped
    for (size_t channel = 0; channel < input_channels; channel++) {
        const size_t current_ring_buffer_pos = ring_buffer_pos[channel];

        // If our ring buffer got filled, then we can do another round of FFT
        // processing
        // HACK: This only works because we can divide `fft_window_size` by our
        //       buffer size
        if (current_ring_buffer_pos == 0 && !first_iteration) {
            std::copy(ring_buffers[channel].begin(),
                      ring_buffers[channel].end(),
                      fft_scratch_buffer[channel].begin());

            fft.performRealOnlyForwardTransform(
                fft_scratch_buffer[channel].data(), true);

            // We'll compress every FTT bin individually
            for (size_t i = 0; i < fft_scratch_buffer[channel].size(); i += 2) {
                // TODO: Or are these reversed? Doesn't really matter
                const float real = fft_scratch_buffer[channel][i];
                const float imag = fft_scratch_buffer[channel][i + 1];
                const float magnitude =
                    std::sqrt((real * real) + (imag * imag));

                // The real and imaginary parts are interleaved, so ever bin
                // spans two values in the scratch buffer
                const size_t compressor_idx = i / 2;
                const float compressed_magnitude =
                    spectral_compressors[compressor_idx].processSample(
                        channel, magnitude);

                // We need to scale both components by the same value
                const float compression_multiplier =
                    magnitude != 0.0f ? compressed_magnitude / magnitude : 1.0f;

                fft_scratch_buffer[channel][i] *= compression_multiplier;
                fft_scratch_buffer[channel][i + 1] *= compression_multiplier;
            }

            fft.performRealOnlyInverseTransform(
                fft_scratch_buffer[channel].data());
        }

        // Copy the input audio into our ring buffer and increment the position
        // (we can do this because we already store this iteration's position)
        float* input_samples = buffer.getWritePointer(channel);
        std::copy(input_samples, input_samples + num_samples,
                  &ring_buffers[channel][current_ring_buffer_pos]);
        ring_buffer_pos[channel] += num_samples;
        if (ring_buffer_pos[channel] >= ring_buffers[channel].size()) {
            // TODO: This too of course only works when things divide cleanly!
            ring_buffer_pos[channel] = 0;
        }

        // Copy over the audio from the previous FFT processing cycle.
        // Everything uses the same ring buffer position.
        // TODO: This too of course only works when the buffer size is cleanly
        //       divisible by the FFT window size
        // TODO: The gain multiplier is for makeup gain. This of course needs to
        //       be dependent on the compressor settings
        buffer.copyFrom(channel, 0,
                        &fft_scratch_buffer[channel][current_ring_buffer_pos],
                        num_samples, 10.0f);

        first_iteration = false;
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
