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
          juce::dsp::WindowingFunction<float>::WindowingMethod::hann,
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
    // TODO: The user should be able to configure their own slope (or free
    //       drawn)
    // TODO: Setting the thresholds based on a sidechain signal would be super
    //       cool
    // TODO: And we should be doing both upwards and downwards compression,
    //       OTT-style
    juce::dsp::Compressor<float> compressor{};
    compressor.setRatio(50.0);
    compressor.setAttack(50.0);
    compressor.setRelease(5000.0);
    compressor.prepare(juce::dsp::ProcessSpec{
        // We only process everything once every `windowing_interval`, otherwise
        // our attack and release times will be all messed up
        .sampleRate = sampleRate / windowing_interval,
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

    // We use ring buffers to store the samples we'll process using FFT and also
    // to store the samples that should be played back to.
    input_ring_buffers.resize(static_cast<size_t>(getTotalNumInputChannels()),
                              std::vector(fft_window_size, 0.0f));
    output_ring_buffers.resize(static_cast<size_t>(getTotalNumInputChannels()),
                               std::vector(fft_window_size, 0.0f));
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

    // Zero out all unused channels
    for (auto i = input_channels; i < output_channels; i++) {
        buffer.clear(i, 0, buffer.getNumSamples());
    }

    // TODO: Add oversampling, potentially reduce latency
    // FIXME: Handling arbitrary buffer sizes almost works! But the latency
    //        compensation is incorrect with buffer sizes that aren't a power of
    //        2 >= 512, and odd sizes between 1024 and 2048 samples/buffer cause
    //        segfaults or artifacts
    for (size_t channel = 0; channel < input_channels; channel++) {
        // Since with large audio buffers we could be processing multiple
        // windows in a single audio processing call, we'll increase this
        // position as we process windows and write the resulting value back at
        // the end of this loop
        size_t current_ring_buffer_pos = ring_buffer_pos[channel];

        float* sample_buffer = buffer.getWritePointer(channel);

        // We process incoming audio in windows of `windowing_interval`, and
        // when using non-power of 2 buffer sizes of buffers that are smaller
        // than `windowing_interval` it can happen that we have to copy over
        // already processed audio before processing a new window
        const size_t already_processed_samples =
            current_ring_buffer_pos % windowing_interval;
        const size_t samples_to_be_processed =
            num_samples - already_processed_samples;
        const int windows_to_process = std::ceil(
            static_cast<float>(num_samples - already_processed_samples) /
            windowing_interval);

        // Copying from the input buffer to our input ring buffer, copying from
        // our output ring buffer to the output buffer, and clearing the output
        // buffer to prevent feedback is always done in sync
        if (already_processed_samples > 0) {
            // TODO: Improve naming for these things
            const size_t copy_samples = std::min(
                already_processed_samples,
                output_ring_buffers[channel].size() - current_ring_buffer_pos);
            const size_t remaining_samples =
                already_processed_samples - copy_samples;
            std::copy_n(sample_buffer, copy_samples,
                        &input_ring_buffers[channel][current_ring_buffer_pos]);
            std::copy_n(&output_ring_buffers[channel][current_ring_buffer_pos],
                        copy_samples, sample_buffer);
            // And clear out the part of the output buffer we just wrote, so we
            // can write new processed windowed samples there
            std::fill_n(&output_ring_buffers[channel][current_ring_buffer_pos],
                        copy_samples, 0.0f);
            if (remaining_samples > 0) {
                std::copy_n(sample_buffer + copy_samples, remaining_samples,
                            &input_ring_buffers[channel][0]);
                std::copy_n(&output_ring_buffers[channel][0], remaining_samples,
                            sample_buffer + copy_samples);
                std::fill_n(&output_ring_buffers[channel][0], remaining_samples,
                            0.0f);
            }
        }

        // Now if `windows_to_process > 0`, the current ring buffer position
        // will align with a window and we can start doing our FFT magic
        current_ring_buffer_pos += already_processed_samples;
        for (int window_idx = 0; window_idx < windows_to_process;
             window_idx++) {
            // This is actual processing
            {
                const size_t copy_input_samples =
                    std::min(static_cast<size_t>(fft_window_size),
                             input_ring_buffers[channel].size() -
                                 current_ring_buffer_pos);
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

                // We'll compress every FTT bin individually. Bin 0 is the DC
                // offset and should be skipped, and the latter half of the FFT
                // bins should be processed in the same way as the first half
                // but in reverse order. The real and imaginary parts are
                // interleaved, so ever bin spans two values in the scratch
                // buffer. We can 'safely' do this cast so we can use the STL's
                // complex value functions.
                std::span<std::complex<float>> fft_buffer(
                    reinterpret_cast<std::complex<float>*>(
                        fft_scratch_buffer[channel].data()),
                    fft_window_size);
                for (size_t compressor_idx = 0;
                     compressor_idx < spectral_compressors.size();
                     compressor_idx++) {
                    // We don't have a compressor for the first bin
                    const size_t bin_idx = compressor_idx + 1;

                    // TODO: Are these _really_ exactly the same in the second
                    //       half ergo this single magnitude is sufficient?
                    const float magnitude = std::abs(fft_buffer[bin_idx]);
                    const float compressed_magnitude =
                        spectral_compressors[compressor_idx].processSample(
                            channel, magnitude);

                    // We need to scale both the imaginary and real components
                    // of the bins at the start and end of the spectrum by the
                    // same value
                    const float compression_multiplier =
                        magnitude != 0.0f ? compressed_magnitude / magnitude
                                          : 1.0f;

                    // The same operation should be applied to the mirrored bins
                    // at the end of the FFT window, except for if this is the
                    // last bin
                    fft_buffer[bin_idx] *= compression_multiplier;
                    // TODO: Is this mirrored part necessary?
                    if (compressor_idx != spectral_compressors.size() - 1) {
                        const size_t mirrored_bin_idx =
                            fft_window_size - bin_idx;
                        fft_buffer[mirrored_bin_idx] *= compression_multiplier;
                    }
                }

                // TODO: Should we also use this window function after
                //       processing?
                fft.performRealOnlyInverseTransform(
                    fft_scratch_buffer[channel].data());
                windowing_function.multiplyWithWindowingTable(
                    fft_scratch_buffer[channel].data(),
                    fft_scratch_buffer[channel].size());

                // TODO: Makeup gain, and when implementing this, take into
                //       account that the 4x overlap also multiplies the volume
                //       by 4 so we should subtrace 6 dB from the makeup gain

                // After processing the windowed data, we'll add it to our
                // output ring buffer
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

            // Copy the input audio into our ring buffer and copy the processed
            // audio into the output buffer
            const size_t samples_to_process_this_iteration =
                window_idx < (windows_to_process - 1)
                    ? windowing_interval
                    : static_cast<size_t>(
                          static_cast<int>(samples_to_be_processed) -
                          (windowing_interval * window_idx));
            const size_t copy_samples = std::min(
                samples_to_process_this_iteration,
                output_ring_buffers[channel].size() - current_ring_buffer_pos);
            const size_t remaining_samples =
                samples_to_process_this_iteration - copy_samples;
            // TODO: Maybe we should maintain the current offset in some
            //       variable
            std::copy_n(sample_buffer + already_processed_samples +
                            (windowing_interval * window_idx),
                        copy_samples,
                        &input_ring_buffers[channel][current_ring_buffer_pos]);
            std::copy_n(&output_ring_buffers[channel][current_ring_buffer_pos],
                        copy_samples,
                        sample_buffer + already_processed_samples +
                            (windowing_interval * window_idx));
            // And clear out the part of the output buffer we just wrote, so we
            // can write new processed windowed samples there
            std::fill_n(&output_ring_buffers[channel][current_ring_buffer_pos],
                        copy_samples, 0.0f);
            if (remaining_samples > 0) {
                std::copy_n(sample_buffer + already_processed_samples +
                                (windowing_interval * window_idx) +
                                copy_samples,
                            remaining_samples, &input_ring_buffers[channel][0]);
                std::copy_n(&output_ring_buffers[channel][0], remaining_samples,
                            sample_buffer + already_processed_samples +
                                (windowing_interval * window_idx) +
                                copy_samples);
                std::fill_n(&output_ring_buffers[channel][0], remaining_samples,
                            0.0f);
            }

            // Depending on whether this is the last iteration and the audio
            // buffer settings, this can either leave us at the next windowing
            // interval or somewhere inbetween windows
            current_ring_buffer_pos += samples_to_process_this_iteration;
            if (current_ring_buffer_pos >= input_ring_buffers[channel].size()) {
                current_ring_buffer_pos -= input_ring_buffers[channel].size();
            }
        }

        // As mentioned at the start of the loop
        ring_buffer_pos[channel] = current_ring_buffer_pos;
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
