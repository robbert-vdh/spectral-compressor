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

#pragma once

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "ring.h"
#include "utils.h"

/**
 * All of the buffers, compressors and other miscellaneous object we'll need to
 * do our FFT audio processing. This will be used together with
 * `AtomicResizable<T>` so it can be resized depending on the current FFT window
 * settings.
 */
struct ProcessData {
    /**
     * The size of the FFT window used in this `ProcessData` object. We store
     * this here so that we always refer to the correct window size when doing
     * swaps.
     */
    size_t fft_window_size;

    /**
     * The numbers of windows already processed. We use this to reduce clicks by
     * not copying over audio to the output during the first
     * `windowing_overap_times` windows.
     */
    int num_windows_processed;

    /**
     * We'll process the signal with overlapping windows that are added to each
     * other to form the output signal. See `input_ring_buffers` for more
     * information on how we'll do this.
     */
    std::optional<juce::dsp::WindowingFunction<float>> windowing_function;

    /**
     * The FFT processor.
     */
    std::optional<juce::dsp::FFT> fft;

    /**
     * We need a scratch buffer that can contain `fft_window_size * 2` samples.
     */
    std::vector<float> fft_scratch_buffer;

    /**
     * This will contain `(fft_window_size / 2) - 1` compressors. The
     * compressors are already multichannel so we don't need a nested vector
     * here. We'll compress the magnitude of every FFT bin (`sqrt(i^2 + r^2)`)
     * individually, and then scale both the real and imaginary components by
     * the ratio of their magnitude and the compressed value. Bin 0 is the DC
     * offset and the bins in the second half should be processed the same was
     * as the bins in the first half but mirrored.
     */
    std::vector<juce::dsp::Compressor<float>> spectral_compressors;

    /**
     * When setting compressor thresholds based on a sidechain signal we should
     * be taking the average bin magnitudes of all channels. This buffer
     * accumulates `spectral_compressors.size()` threshold values while
     * iterating over the channels of the sidechain signal so we can then
     * average them and configure the compressors based on that.
     */
    std::vector<float> spectral_compressor_sidechain_thresholds;

    /**
     * A ring buffer of size `fft_window_size` for every channel. Every
     * `windowing_interval` we'll copy the last `fft_window_size` samples to
     * `fft_scratch_buffers` using a window function, process it, and then add
     * the results to `output_ring_buffers`.
     */
    std::vector<RingBuffer<float>> input_ring_buffers;
    /**
     * The processed results as described in the docstring of
     * `input_ring_buffers`. Samples from this buffer will be written to the
     * output.
     */
    std::vector<RingBuffer<float>> output_ring_buffers;
    /**
     * These ring buffers are identical to `input_ring_buffers`, but with data
     * from the sidechain input. When sidechaining is enabled, we set the
     * compressor thresholds based on the magnitudes from the same FFT analysis
     * applied to the sidechain input.
     */
    std::vector<RingBuffer<float>> sidechain_ring_buffers;
};

/**
 * Run some function on the message thread. This function will be executed
 * synchronously and should thus run in constant time.
 */
class LambdaAsyncUpdater : public juce::AsyncUpdater {
   public:
    LambdaAsyncUpdater(fu2::unique_function<void()> callback);

    void handleAsyncUpdate() override;

   private:
    fu2::unique_function<void()> callback;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LambdaAsyncUpdater)
};

/**
 * Run some function whenever a parameter changes. This function will be
 * executed synchronously and should thus run in constant time.
 */
class LambdaParameterListener
    : public juce::AudioProcessorValueTreeState::Listener {
   public:
    LambdaParameterListener(
        fu2::unique_function<void(const juce::String&, float)> callback);

    void parameterChanged(const juce::String& parameterID,
                          float newValue) override;

   private:
    fu2::unique_function<void(const juce::String&, float)> callback;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LambdaParameterListener)
};

class SpectralCompressorProcessor : public juce::AudioProcessor {
   public:
    SpectralCompressorProcessor();
    ~SpectralCompressorProcessor() override;

    void prepareToPlay(double sampleRate,
                       int maximumExpectedSamplesPerBlock) override;
    void releaseResources() override;

    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;

    void processBlockBypassed(juce::AudioBuffer<float>& buffer,
                              juce::MidiBuffer& midiMessages) override;
    using AudioProcessor::processBlockBypassed;
    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    using AudioProcessor::processBlock;

    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram(int index) override;
    const juce::String getProgramName(int index) override;
    void changeProgramName(int index, const juce::String& newName) override;

    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;

   private:
    /**
     * (Re)initialize a process data object for the given FFT order. If the new
     * FFT order is 0, then the object will be cleared instead.
     */
    void initialize_process_data(ProcessData& inactive, size_t new_fft_order);

    /**
     * Calculate new compressor thresholds and other settings based on the
     * current parameters.
     */
    void update_compressors(ProcessData& data);

    /**
     * Process audio using a short term Fourier transform. This involves using
     * the input ring buffers of `ProcessingData` to buffer audio, processing
     * that audio in windows, adding up those windows in the output ring
     * buffers, and then writing those outputs to `buffer`'s outputs. This
     * function handles all of the boilerplate outside of the actual FFT
     * operations.
     *
     * @param buffer The current processing cycle's buffers. This should contain
     *   input, output, and sidechain busses with an equal number of channels
     *   for each bus.
     * @param data The current processing data.
     * @param process_fn A function that performs a forward FFT on the input
     *   ring buffer, processes the results, performs an IFFT, and then adds the
     *   results to the output ring buffer using a windowing function.
     *
     * @tparam F A function that takes the current processing data, and the
     *   number of input channels as its arguments.
     */
    template <std::invocable<ProcessData&, size_t> F>
    void do_stft(juce::AudioBuffer<float>& buffer,
                 ProcessData& data,
                 F process_fn) {
        juce::ScopedNoDenormals noDenormals;

        juce::AudioBuffer<float> main_io = getBusBuffer(buffer, true, 0);
        juce::AudioBuffer<float> sidechain_io = getBusBuffer(buffer, true, 1);

        const size_t input_channels =
            static_cast<size_t>(getMainBusNumInputChannels());
        const size_t output_channels =
            static_cast<size_t>(getMainBusNumOutputChannels());
        const size_t num_samples = static_cast<size_t>(buffer.getNumSamples());

        // Zero out all unused channels
        for (auto channel = input_channels; channel < output_channels;
             channel++) {
            buffer.clear(channel, 0.0f, num_samples);
        }

        // We'll process audio in lockstep to make it easier to use processors
        // that require lookahead and thus induce latency. Every this many
        // samples we'll process a new window of input samples. The results will
        // be added to the output ring buffers.
        const size_t windowing_interval =
            data.fft_window_size / static_cast<size_t>(windowing_overlap_times);

        // We process incoming audio in windows of `windowing_interval`, and
        // when using non-power of 2 buffer sizes of buffers that are smaller
        // than `windowing_interval` it can happen that we have to copy over
        // already processed audio before processing a new window
        const size_t already_processed_samples =
            std::min(num_samples,
                     (windowing_interval -
                      (data.input_ring_buffers[0].pos() % windowing_interval)) %
                         windowing_interval);
        const size_t samples_to_be_processed =
            num_samples - already_processed_samples;
        const int windows_to_process = std::ceil(
            static_cast<float>(samples_to_be_processed) / windowing_interval);

        // Since we're processing audio in small chunks, we need to keep track
        // of the current sample offset in `buffers` we should use for our
        // actual audio input and output
        size_t sample_buffer_offset = 0;

        // Copying from the input buffer to our input ring buffer, copying from
        // our output ring buffer to the output buffer, and clearing the output
        // buffer to prevent feedback is always done in sync
        if (already_processed_samples > 0) {
            for (size_t channel = 0; channel < input_channels; channel++) {
                data.input_ring_buffers[channel].read_n_from(
                    main_io.getReadPointer(channel), already_processed_samples);
                if (data.num_windows_processed >= windowing_overlap_times) {
                    data.output_ring_buffers[channel].copy_n_to(
                        main_io.getWritePointer(channel),
                        already_processed_samples, true);
                } else {
                    main_io.clear(channel, 0, already_processed_samples);
                }
                if (sidechain_active) {
                    data.sidechain_ring_buffers[channel].read_n_from(
                        sidechain_io.getReadPointer(channel),
                        already_processed_samples);
                }
            }

            sample_buffer_offset += already_processed_samples;
        }

        // We'll update the compressor settings just before processing if the
        // settings have changed or if the sidechaining has been disabled.
        bool expected = true;
        if (windows_to_process > 0 &&
            compressor_settings_changed.compare_exchange_strong(expected,
                                                                false)) {
            update_compressors(data);
        }

        // Now if `windows_to_process > 0`, the current ring buffer position
        // will align with a window and we can start doing our FFT magic
        for (int window_idx = 0; window_idx < windows_to_process;
             window_idx++) {
            // This is where the actual processing happens. This function should
            // take
            process_fn(data, input_channels);

            // We don't copy over anything to the outputs until we processed a
            // full buffer
            data.num_windows_processed += 1;

            // Copy the input audio into our ring buffer and copy the processed
            // audio into the output buffer
            const size_t samples_to_process_this_iteration = std::min(
                windowing_interval, num_samples - sample_buffer_offset);
            for (size_t channel = 0; channel < input_channels; channel++) {
                data.input_ring_buffers[channel].read_n_from(
                    main_io.getReadPointer(channel) + sample_buffer_offset,
                    samples_to_process_this_iteration);
                if (data.num_windows_processed >= windowing_overlap_times) {
                    data.output_ring_buffers[channel].copy_n_to(
                        main_io.getWritePointer(channel) + sample_buffer_offset,
                        samples_to_process_this_iteration, true);
                } else {
                    main_io.clear(channel, sample_buffer_offset,
                                  samples_to_process_this_iteration);
                }
                if (sidechain_active) {
                    data.sidechain_ring_buffers[channel].read_n_from(
                        sidechain_io.getReadPointer(channel) +
                            sample_buffer_offset,
                        samples_to_process_this_iteration);
                }
            }

            sample_buffer_offset += samples_to_process_this_iteration;
        }

        jassert(sample_buffer_offset == num_samples);
    }

    /**
     * This contains all of our scratch buffers, ring buffers, compressors, and
     * everything else that depends on the FFT window size.
     */
    AtomicResizable<ProcessData> process_data;

    /**
     * Will be set during `prepareToPlay()`, needed to initialize compressors
     * when resizing our buffers.
     */
    juce::uint32 max_samples_per_block;

    juce::AudioProcessorValueTreeState parameters;

    juce::AudioParameterBool& sidechain_active;
    std::atomic<float>& compressor_ratio;
    /**
     * Try to automatically compensate for low thresholds. Doesn't do anything
     * when sidechaining is active.
     */
    juce::AudioParameterBool& auto_makeup_gain;

    /**
     * Will be set in `CompressorSettingsListener` when any of the compressor
     * related settings change so we can update our compressors. We'll
     * initialize this to true so the compressors will be initialized during the
     * first processing cycle.
     */
    std::atomic_bool compressor_settings_changed = true;
    /**
     * Makeup gain to be applied after compression, where 1.0 mean no gain
     * applied. Depends on the current active modes and whether the makeup gain
     * parameters.
     *
     * The computed value also takes the 4x overlap into account.
     */
    float makeup_gain;

    /**
     * The order (where `fft_window_size = 1 << fft_order`) for our spectral
     * operations. When this gets changed, we'll resize all of our buffers and
     * atomically swap the current and the resized buffers.
     */
    juce::AudioParameterInt& fft_order;
    /**
     * The amount of overlap in the windowing. We end up processing the signal
     * in `fft_window_size` windows every `fft_window_size /
     * windowing_overlap_times` samples. When this setting gets changed, we'll
     * also have to update our compressors since the effective sample rate also
     * changes.
     */
    juce::AudioParameterInt& windowing_overlap_times;

    /**
     * Will cause the compressor settings to be updated on the next processing
     * cycle whenever a compressor parameter changes.
     */
    LambdaParameterListener compressor_settings_listener;

    /**
     * When the FFT order parameter changes, we'll have to create a new
     * `ProcessData` object for the new FFT window size (or rather, resize an
     * inactive one to match the new size).
     */
    LambdaParameterListener fft_order_listener;

    /**
     * Atomically resizes the object `ProcessData` from a background thread.
     */
    LambdaAsyncUpdater process_data_resizer;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpectralCompressorProcessor)
};
