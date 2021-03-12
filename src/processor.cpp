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

#include <span>

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
      windowing_function(
          fft_window_size,
          juce::dsp::WindowingFunction<float>::WindowingMethod::hamming,
          // TODO: Or should we leave normalization enabled?
          false),
      fft(fft_order) {
    setLatencySamples(fft_window_size);
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

    // We use ring buffers to fill our FFT buffers and to write the results to
    // after we process windows of the input. They need to be able to contain as
    // many samples as we'll process with FFT plus an additional window so we
    // can clear out our buffers.
    // TODO: Do we need this extra capacity?
    input_ring_buffers.resize(
        static_cast<size_t>(getTotalNumInputChannels()),
        std::vector(fft_window_size + windowing_interval, 0.0f));
    output_ring_buffers.resize(
        static_cast<size_t>(getTotalNumInputChannels()),
        std::vector(fft_window_size + windowing_interval, 0.0f));
    ring_buffer_pos.resize(static_cast<size_t>(getTotalNumInputChannels()));
}

void SpectralCompressorProcessor::releaseResources() {
    fft_scratch_buffer.clear();
    spectral_compressors.clear();
    input_ring_buffers.clear();
    output_ring_buffers.clear();
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

    // FIXME: Since we're cutting some corners still this might only work at 512
    //        or 1024 samples/buffer
    jassert((windowing_interval % num_samples) == 0);

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
        if ((current_ring_buffer_pos % windowing_interval) == 0) {
            const size_t copy_input_samples = std::min(
                static_cast<size_t>(fft_window_size),
                input_ring_buffers[channel].size() - current_ring_buffer_pos);
            const size_t remaining_input_samples =
                fft_window_size - copy_input_samples;
            std::copy_n(input_ring_buffers[channel].begin() +
                            static_cast<int>(current_ring_buffer_pos),
                        copy_input_samples,
                        fft_scratch_buffer[channel].begin());
            if (remaining_input_samples > 0) {
                std::copy_n(input_ring_buffers[channel].begin(),
                            remaining_input_samples,
                            fft_scratch_buffer[channel].begin() +
                                static_cast<int>(copy_input_samples));
            }

            windowing_function.multiplyWithWindowingTable(
                fft_scratch_buffer[channel].data(),
                fft_scratch_buffer[channel].size());
            fft.performRealOnlyForwardTransform(
                fft_scratch_buffer[channel].data());

            // We'll compress every FTT bin individually. Bin 0 is the DC offset
            // and should be skipped, and the latter half of the FFT bins should
            // be processed in the same way as the first half but in reverse
            // order.
            // The real and imaginary parts are interleaved, so ever bin spans
            // two values in the scratch buffer. We can 'safely' do this cast so
            // we can use the STL's complex value functions.
            std::span<std::complex<float>> fft_buffer(
                reinterpret_cast<std::complex<float>*>(
                    fft_scratch_buffer[channel].data()),
                fft_window_size);
            for (size_t compressor_idx = 0;
                 compressor_idx < spectral_compressors.size();
                 compressor_idx++) {
                // We don't have a compressor for the first bin
                const size_t bin_idx = compressor_idx + 1;

                // TODO: Are these _really_ exactly the same in the second half
                //       ergo this single magnitude is sufficient?
                const float magnitude = std::abs(fft_buffer[bin_idx]);
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
                fft_buffer[bin_idx] *= compression_multiplier;
                // TODO: Is this mirrored part necessary?
                if (compressor_idx != spectral_compressors.size() - 1) {
                    const size_t mirrored_bin_idx = fft_window_size - bin_idx;
                    fft_buffer[mirrored_bin_idx] *= compression_multiplier;
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
            windowing_function.multiplyWithWindowingTable(
                fft_scratch_buffer[channel].data(),
                fft_scratch_buffer[channel].size());

            // TODO: Makeup gain

            // After processing the windowed data, we'll add it to our output
            // ring buffer.
            // TODO: We probably need to increase this buffer for more overlap?
            juce::FloatVectorOperations::add(
                &output_ring_buffers[channel][current_ring_buffer_pos],
                fft_scratch_buffer[channel].data(), copy_input_samples);
            if (remaining_input_samples > 0) {
                juce::FloatVectorOperations::add(
                    &output_ring_buffers[channel][0],
                    &fft_scratch_buffer[channel][copy_input_samples],
                    remaining_input_samples);
            }
        }

        // Copy the input audio into our ring buffer and increment the position
        // (we can do this because we already store this iteration's position)
        // FIXME: This is of course only safe because of the assertion at the
        //        start of the function! Otherwise we might have to split this
        //        up.
        float* input_samples = buffer.getWritePointer(channel);
        std::copy_n(input_samples, num_samples,
                    &input_ring_buffers[channel][current_ring_buffer_pos]);
        ring_buffer_pos[channel] += num_samples;
        if (ring_buffer_pos[channel] >= input_ring_buffers[channel].size()) {
            // TODO: This too of course only works when things divide cleanly!
            ring_buffer_pos[channel] = 0;
        }

        // Copy over the audio from the previous the added processed windows
        // TODO: This too of course only works when the buffer size is cleanly
        //       divisible by the FFT window size
        buffer.copyFrom(channel, 0,
                        &output_ring_buffers[channel][current_ring_buffer_pos],
                        num_samples);

        // Clear out the part of the output buffer we just output, so we can
        // write new processed windowed samples there
        // TODO: Where should we clear this out? It kind of makes sense to clear
        //       out this window of the output buffer just after outputting it
        std::fill_n(&output_ring_buffers[channel][current_ring_buffer_pos],
                    num_samples, 0.0f);
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
