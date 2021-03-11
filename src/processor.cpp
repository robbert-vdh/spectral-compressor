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
    // The `(fft_window_size / 2) - 1` is because the first bin is the DC offset
    // and shouldn't be compressed, and the bins after the Nyquist frequency are
    // the same as the first half but in reverse order.
    // TODO: Make the compressor settings configurable
    // TODO: These settings are also very extreme
    juce::dsp::Compressor<float> compressor{};
    compressor.setRatio(50.0);
    // TODO: I don't think the attack and release times are correct with the way
    //       we only sporadically use our compressors
    compressor.setAttack(10.0);
    compressor.setRelease(50.0);
    compressor.prepare(juce::dsp::ProcessSpec{
        .sampleRate = sampleRate,
        .maximumBlockSize = static_cast<uint32>(maximumExpectedSamplesPerBlock),
        .numChannels = static_cast<uint32>(getMainBusNumInputChannels())});

    spectral_compressors.resize((fft_window_size / 2) - 1, compressor);

    // The thresholds are set to match pink noise.
    constexpr float base_threshold_dbfs = 0.0f;
    const float frequency_increment = sampleRate / fft_window_size;
    for (size_t compressor_idx = 0;
         compressor_idx < spectral_compressors.size(); compressor_idx++) {
        // The first bin doesn't get a compressor
        const size_t bin_idx = compressor_idx + 1;
        const float frequency = frequency_increment * bin_idx;

        // This starts at 1 for 0 Hz (DC)
        const float octave = std::log2(frequency + 2);

        // The 3 dB is to compensate for bin 0
        const float threshold = (base_threshold_dbfs + 3.0f) - (3.0f * octave);
        spectral_compressors[compressor_idx].setThreshold(threshold);

        std::cerr << "Compressor " << compressor_idx << ", frequency "
                  << frequency << " Hz, octave " << octave << ", threshold "
                  << threshold << " dBFS" << std::endl;
    }

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
        if (current_ring_buffer_pos == 0) {
            std::copy(ring_buffers[channel].begin(),
                      ring_buffers[channel].end(),
                      fft_scratch_buffer[channel].begin());

            fft.performRealOnlyForwardTransform(
                fft_scratch_buffer[channel].data());

            // We'll compress every FTT bin individually. Bin 0 is the DC offset
            // and should be skipped, and the latter half of the FFT bins should
            // be processed in the same way as the first half but in reverse
            // order.
            const size_t num_bins = fft_scratch_buffer[channel].size();
            for (size_t compressor_idx = 0;
                 compressor_idx < spectral_compressors.size();
                 compressor_idx++) {
                const size_t bin_idx = (compressor_idx * 2) + 2;

                // The real and imaginary parts are interleaved, so ever bin
                // spans two values in the scratch buffer
                // TODO: Or are these reversed? Doesn't really matter
                // TODO: Are these _really_ exactly the same in the second half
                //       ergo this single magnitude is sufficient?
                const float real = fft_scratch_buffer[channel][bin_idx];
                const float imag = fft_scratch_buffer[channel][bin_idx + 1];
                const float magnitude =
                    std::sqrt((real * real) + (imag * imag));

                const float compressed_magnitude =
                    spectral_compressors[compressor_idx].processSample(
                        channel, magnitude);

                // We need to scale both the imaginary and real components of
                // the bins at the start and end of the spectrum by the same
                // value
                const float compression_multiplier =
                    magnitude != 0.0f ? compressed_magnitude / magnitude : 1.0f;

                // FIXME: Remove this after everything's A-Ok
                if (!std::isnormal(compression_multiplier)) {
                    std::cerr << "Skipping multiplier "
                              << compression_multiplier << " @ " << channel
                              << ":" << compressor_idx << std::endl;
                    continue;
                }

                // The same operation should be applied to the mirrored bins at
                // the end of the FFT window, except for if this is the last bin
                fft_scratch_buffer[channel][bin_idx] *= compression_multiplier;
                fft_scratch_buffer[channel][bin_idx + 1] *=
                    compression_multiplier;
                if (compressor_idx != spectral_compressors.size() - 1) {
                    fft_scratch_buffer[channel][num_bins - bin_idx] *=
                        compression_multiplier;
                    fft_scratch_buffer[channel][num_bins - bin_idx + 1] *=
                        compression_multiplier;
                }
            }

            // FIXME: Remove this after everything's A-Ok
            for (size_t i = 0; i < fft_scratch_buffer[channel].size(); i++) {
                if (fft_scratch_buffer[channel][i] != 0 &&
                    !std::isnormal(fft_scratch_buffer[channel][i])) {
                    std::cerr << "Post-FFT non-normal "
                              << fft_scratch_buffer[channel][i] << " @ "
                              << channel << ":" << i << std::endl;
                    switch (std::fpclassify(fft_scratch_buffer[channel][i])) {
                        case FP_INFINITE:
                            std::cerr << "Inf" << std::endl;
                            break;
                        case FP_NAN:
                            std::cerr << "NaN" << std::endl;
                            break;
                        case FP_NORMAL:
                            std::cerr << "normal" << std::endl;
                            break;
                        case FP_SUBNORMAL:
                            std::cerr << "subnormal" << std::endl;
                            break;
                        case FP_ZERO:
                            std::cerr << "zero" << std::endl;
                            break;
                        default:
                            std::cerr << "unknown" << std::endl;
                            break;
                    }
                    fft_scratch_buffer[channel][i] = 0.0f;
                }
            }

            fft.performRealOnlyInverseTransform(
                fft_scratch_buffer[channel].data());

            // TODO: Makeup gain
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
        buffer.copyFrom(channel, 0,
                        &fft_scratch_buffer[channel][current_ring_buffer_pos],
                        num_samples);
    }
}

//==============================================================================
bool SpectralCompressorProcessor::hasEditor() const {
    // TODO: Add an editor at some point
    return false;
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
