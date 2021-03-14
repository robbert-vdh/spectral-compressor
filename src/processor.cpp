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
    fft_scratch_buffer.resize(fft_window_size * 2, 0.0f);

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
    // TODO: Change the calculations so that the base threshold parameter is
    //       centered around some frequency
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
    }

    // We use ring buffers to store the samples we'll process using FFT and also
    // to store the samples that should be played back to.
    input_ring_buffers.resize(static_cast<size_t>(getTotalNumInputChannels()),
                              RingBuffer<float>(fft_window_size));
    output_ring_buffers.resize(static_cast<size_t>(getTotalNumInputChannels()),
                               RingBuffer<float>(fft_window_size));
}

void SpectralCompressorProcessor::releaseResources() {
    // TODO: Clearing a vector actually doesn't really do anything. So either
    //       don't do anything here or actually swap these vectors with new
    //       vectors/shrink to fit.
    fft_scratch_buffer.clear();
    spectral_compressors.clear();
    input_ring_buffers.clear();
    output_ring_buffers.clear();
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

void SpectralCompressorProcessor::processBlockBypassed(
    juce::AudioBuffer<float>& buffer,
    juce::MidiBuffer& /*midiMessages*/) {
    // We need to maintain the same latency when bypassed, so we'll reuse most
    // of the processing logic
    process(buffer, true);
}

void SpectralCompressorProcessor::processBlock(
    juce::AudioBuffer<float>& buffer,
    juce::MidiBuffer& /*midiMessages*/) {
    process(buffer, false);
}

bool SpectralCompressorProcessor::hasEditor() const {
    // TODO: Add an editor at some point
    return false;
}

juce::AudioProcessorEditor* SpectralCompressorProcessor::createEditor() {
    return new SpectralCompressorEditor(*this);
}

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

void SpectralCompressorProcessor::process(juce::AudioBuffer<float>& buffer,
                                          bool bypassed) {
    juce::ScopedNoDenormals noDenormals;

    const size_t input_channels =
        static_cast<size_t>(getMainBusNumInputChannels());
    const size_t output_channels =
        static_cast<size_t>(getTotalNumOutputChannels());
    const size_t num_samples = static_cast<size_t>(buffer.getNumSamples());

    // Zero out all unused channels
    for (auto channel = input_channels; channel < output_channels; channel++) {
        buffer.clear(channel, 0.0f, num_samples);
    }

    // TODO: Add oversampling, potentially reduce latency
    for (size_t channel = 0; channel < input_channels; channel++) {
        // The sample buffer contains input audio and anything written to it
        // will be passed to the host. We'll keep track of the current offset in
        // the sample buffer since we'll be reading and writing in chunks.
        float* sample_buffer = buffer.getWritePointer(channel);
        size_t sample_buffer_offset = 0;

        // We process incoming audio in windows of `windowing_interval`, and
        // when using non-power of 2 buffer sizes of buffers that are smaller
        // than `windowing_interval` it can happen that we have to copy over
        // already processed audio before processing a new window
        const size_t already_processed_samples = std::min(
            num_samples,
            (windowing_interval -
             (input_ring_buffers[channel].pos() % windowing_interval)) %
                windowing_interval);
        const size_t samples_to_be_processed =
            num_samples - already_processed_samples;
        const int windows_to_process = std::ceil(
            static_cast<float>(samples_to_be_processed) / windowing_interval);

        // Copying from the input buffer to our input ring buffer, copying from
        // our output ring buffer to the output buffer, and clearing the output
        // buffer to prevent feedback is always done in sync
        if (already_processed_samples > 0) {
            input_ring_buffers[channel].read_n_from(sample_buffer,
                                                    already_processed_samples);
            output_ring_buffers[channel].copy_n_to(
                sample_buffer, already_processed_samples, true);
            sample_buffer_offset += already_processed_samples;
        }

        // Now if `windows_to_process > 0`, the current ring buffer position
        // will align with a window and we can start doing our FFT magic
        for (int window_idx = 0; window_idx < windows_to_process;
             window_idx++) {
            // This is actual processing
            if (!bypassed) {
                input_ring_buffers[channel].copy_last_n_to(
                    fft_scratch_buffer.data(), fft_window_size);

                windowing_function.multiplyWithWindowingTable(
                    fft_scratch_buffer.data(), fft_scratch_buffer.size());
                fft.performRealOnlyForwardTransform(fft_scratch_buffer.data());

                // We'll compress every FTT bin individually. Bin 0 is the DC
                // offset and should be skipped, and the latter half of the FFT
                // bins should be processed in the same way as the first half
                // but in reverse order. The real and imaginary parts are
                // interleaved, so ever bin spans two values in the scratch
                // buffer. We can 'safely' do this cast so we can use the STL's
                // complex value functions.
                std::span<std::complex<float>> fft_buffer(
                    reinterpret_cast<std::complex<float>*>(
                        fft_scratch_buffer.data()),
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
                fft.performRealOnlyInverseTransform(fft_scratch_buffer.data());
                windowing_function.multiplyWithWindowingTable(
                    fft_scratch_buffer.data(), fft_scratch_buffer.size());

                // TODO: Makeup gain, and when implementing this, take into
                //       account that the 4x overlap also multiplies the volume
                //       by 4 so we should subtrace 6 dB from the makeup gain.
                //       We should probably apply this when copying data from
                //       the output ring buffer to the sample buffer.

                // After processing the windowed data, we'll add it to our
                // output ring buffer
                output_ring_buffers[channel].add_n_from_in_place(
                    fft_scratch_buffer.data(), fft_window_size);
            } else {
                // We don't have a way to directly copy between buffers, but
                // most hosts should not actually hit this bypassed state
                // anyways
                // TODO: At some point, do implement this without using the
                //       scratch buffer
                input_ring_buffers[channel].copy_last_n_to(
                    fft_scratch_buffer.data(), fft_window_size);
                output_ring_buffers[channel].read_n_from_in_place(
                    fft_scratch_buffer.data(), fft_window_size);
            }

            // Copy the input audio into our ring buffer and copy the processed
            // audio into the output buffer
            const size_t samples_to_process_this_iteration = std::min(
                windowing_interval, num_samples - sample_buffer_offset);
            input_ring_buffers[channel].read_n_from(
                sample_buffer + sample_buffer_offset,
                samples_to_process_this_iteration);
            output_ring_buffers[channel].copy_n_to(
                sample_buffer + sample_buffer_offset,
                samples_to_process_this_iteration, true);
            sample_buffer_offset += samples_to_process_this_iteration;
        }

        jassert(sample_buffer_offset == num_samples);
    }
}

//==============================================================================
// This creates new instances of the plugin..
juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
    return new SpectralCompressorProcessor();
}
