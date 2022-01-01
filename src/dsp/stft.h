// Spectral Compressor: an FFT based compressor
// Copyright (C) 2021-2022 Robbert van der Helm
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

#include <span>

#include <juce_dsp/juce_dsp.h>

#include "../ring.h"

/**
 * Process an audio source in the frequency domain using the overlap-add method.
 *
 * @tparam with_sidechain Whether to also do a parallel analysis on a sidechain
 *   input source. If this is enabled, then you will be able to analyze an FFT
 *   buffer from a sidechain source before processing the main signal.
 */
template <bool with_sidechain = false>
class STFT {
   public:
    /**
     * Initialize this processor for the given FFT order.
     *
     * @param num_channels The number of channels. This should be equal for the
     *   input, sidechain, and output busses.
     * @param fft_order The order of the FFT window. The actual size of the
     *   window will be `1 << fft_order`.
     */
    STFT(size_t num_channels, size_t fft_order)
        : fft_window_size(1 << fft_order),
          fft_(fft_order),
          windowing_function_(
              fft_window_size,
              juce::dsp::WindowingFunction<float>::WindowingMethod::hann,
              // TODO: Or should we leave normalization enabled?
              false),
          // JUCE's FFT class interleaves the real and imaginary numbers, so
          // this buffer should be twice the window size in size
          fft_scratch_buffer_(fft_window_size * 2),
          input_ring_buffers_(num_channels, RingBuffer<float>(fft_window_size)),
          sidechain_ring_buffers_(with_sidechain ? num_channels : 0,
                                  with_sidechain
                                      ? RingBuffer<float>(fft_window_size)
                                      : RingBuffer<float>()),
          output_ring_buffers_(num_channels,
                               RingBuffer<float>(fft_window_size)) {}

    /**
     * The latency introduced by this processor, in samples.
     */
    inline int latency_samples() const { return fft_window_size; }

    /**
     * Process audio using a short term Fourier transform. This involves using
     * the input ring buffers to buffer audio, processing that audio in windows,
     * adding up those windows in the output ring buffers, and then finally
     * writing those outputs to `buffer`'s outputs. The supplied function can be
     * used to actually process the data.
     *
     * @param main_io The current processing cycle's buffers for the main input
     *   and output busses. This should contain an input and an output bus with
     *   an equal number of channels for each bus.
     * @param windowing_overlap_times How much overlap we should be using in the
     *   overlap-add process. This should be a power of two.
     * @param gain Gain to apply to every processed window before adding it to
     *   the output. If set to 1.0, no gain will be added.
     * @param preprocess_fn A function that receives a window of raw samples
     *   just before the FFT processing. The windowing function will have
     *   already applied at this point.
     * @param process_fn A function that receives and modifies an FFT buffer.
     *   The results will be written back to `buffer`'s outputs using the
     *   overlap-add method at an `fft_window_size` sample delay.
     * @param postprocess_fn A function that receives raw samples just after the
     *   FFT processing but before they are added to the output ring buffers.
     *   Windowing will have already been applied at this point.
     *
     * @tparam FPreProcess A function of type `void(std::span<float>& fft,
     *   size_t channel)`.
     * @tparam FProcess A function of type `void(std::span<std::complex<float>>&
     *   fft, size_t channel)`.
     * @tparam FPostProcess A function of type `void(std::span<float>& fft,
     *   size_t channel)`.
     */
    template <typename FPreProcess, typename FProcess, typename FPostProcess>
    void process(juce::AudioBuffer<float>& main_io,
                 int windowing_overlap_times,
                 float gain,
                 FPreProcess preprocess_fn,
                 FProcess process_fn,
                 FPostProcess postprocess_fn) {
        do_process<false, false>(
            main_io, main_io, windowing_overlap_times, gain, [](auto&, auto) {},
            []() {}, std::move(preprocess_fn), std::move(process_fn),
            std::move(postprocess_fn));
    }

    /**
     * Process audio using a short term Fourier transform. This involves using
     * the input ring buffers to buffer audio, processing that audio in windows,
     * adding up those windows in the output ring buffers, and then finally
     * writing those outputs to `buffer`'s outputs. The supplied function can be
     * used to actually process the data.
     *
     * This version lets you analyze a sidechain signal before processing the
     * main signal.
     *
     * @param main_io The current processing cycle's buffers for the main input
     *   and output busses. This should contain an input and an output bus with
     *   an equal number of channels for each bus.
     * @param sidechain_io The current processing cycle's buffers for the
     *   sidechain input busses. This should have the same number of channels as
     *   `main_io`.
     * @param windowing_overlap_times How much overlap we should be using in the
     *   overlap-add process. This should be a power of two.
     * @param gain Gain to apply to every processed window before adding it to
     *   the output. If set to 1.0, no gain will be added.
     * @param sidechain_fn A function that receives an FFT buffer obtained from
     *   the sidechain signal that can be used for analysis.
     * @param post_sidechain_fn A function called after `sidechain_fn` has been
     *   called for every channel. Can be used to aggregate per-channel data.
     * @param preprocess_fn A function that receives a window of raw samples
     *   just before the FFT processing. The windowing function will have
     *   already applied at this point.
     * @param process_fn A function that receives and modifies an FFT buffer.
     *   The results will be written back to `buffer`'s outputs using the
     *   overlap-add method at an `fft_window_size` sample delay.
     * @param postprocess_fn A function that receives raw samples just after the
     *   FFT processing but before they are added to the output ring buffers.
     *   Windowing will have already been applied at this point.
     *
     * @tparam FSidechain A function of type `void(const
     *   std::span<std::complex<float>>& fft, size_t channel)`.
     * @tparam FPostSidechain A `void()` function.
     * @tparam FPreProcess A function of type `void(std::span<float>& fft,
     *   size_t channel)`.
     * @tparam FProcess A function of type `void(std::span<std::complex<float>>&
     *   fft, size_t channel)`.
     * @tparam FPostProcess A function of type `void(std::span<float>& fft,
     *   size_t channel)`.
     */
    template <typename FSidechain,
              typename FPostSidechain,
              typename FPreProcess,
              typename FProcess,
              typename FPostProcess,
              typename = std::enable_if_t<with_sidechain>>
    void process(juce::AudioBuffer<float>& main_io,
                 const juce::AudioBuffer<float>& sidechain_io,
                 int windowing_overlap_times,
                 float gain,
                 FSidechain sidechain_fn,
                 FPostSidechain post_sidechain_fn,
                 FPreProcess preprocess_fn,
                 FProcess process_fn,
                 FPostProcess postprocess_fn) {
        do_process<false, true>(main_io, sidechain_io, windowing_overlap_times,
                                gain, std::move(sidechain_fn),
                                std::move(post_sidechain_fn),
                                std::move(preprocess_fn), std::move(process_fn),
                                std::move(postprocess_fn));
    }

    /**
     * Don't do any processing, but still keep the same amount of latency as if
     * we were calling `process()`.
     *
     * @param main_io The current processing cycle's buffers for the main input
     *   and output busses. This should contain an input and an output bus with
     *   an equal number of channels for each bus.
     */
    void process_bypassed(juce::AudioBuffer<float>& main_io) {
        do_process<true, false>(
            main_io, main_io, 1, 1.0f, [](auto&, auto) {}, []() {},
            [](auto&, auto) {}, [](auto&, auto) {}, [](auto&, auto) {});
    }

    /**
     * The size of the FFT window used.
     */
    const size_t fft_window_size;

   private:
    /**
     * Depending on `with_sidechain`, there are a few different ways to process
     * a buffer. To avoid duplication, this function has two `bypassed` and
     * `sidechain_active` template constants that control the order through this
     * function. These booleans control whether we do any FFT operations at all,
     * and whether we touch the read from the sidechain input and call the
     * sidechain analysis functions.
     */
    template <bool bypassed,
              bool sidechain_active,
              typename FSidechain,
              typename FPostSidechain,
              typename FPreProcess,
              typename FProcess,
              typename FPostProcess>
    void do_process(
        juce::AudioBuffer<float>& main_io,
        [[maybe_unused]] const juce::AudioBuffer<float>& sidechain_io,
        int windowing_overlap_times,
        float gain,
        [[maybe_unused]] FSidechain sidechain_fn,
        [[maybe_unused]] FPostSidechain post_sidechain_fn,
        FPreProcess preprocess_fn,
        FProcess process_fn,
        FPostProcess postprocess_fn) {
        juce::ScopedNoDenormals noDenormals;

        const size_t num_channels =
            static_cast<size_t>(main_io.getNumChannels());
        const size_t num_samples = static_cast<size_t>(main_io.getNumSamples());
        if constexpr (sidechain_active) {
            jassert(sidechain_io.getNumChannels() ==
                    static_cast<int>(num_channels));
            jassert(sidechain_io.getNumSamples() ==
                    static_cast<int>(num_samples));
        }

        // We'll process audio in lockstep to make it easier to use processors
        // that require lookahead and thus induce latency. Every this many
        // samples we'll process a new window of input samples. The results will
        // be added to the output ring buffers.
        const size_t windowing_interval =
            fft_window_size / static_cast<size_t>(windowing_overlap_times);

        // We process incoming audio in windows of `windowing_interval`, and
        // when using non-power of 2 buffer sizes of buffers that are smaller
        // than `windowing_interval` it can happen that we have to copy over
        // already processed audio before processing a new window
        const size_t already_processed_samples = std::min(
            num_samples, (windowing_interval -
                          (input_ring_buffers_[0].pos() % windowing_interval)) %
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
            for (size_t channel = 0; channel < num_channels; channel++) {
                input_ring_buffers_[channel].read_n_from(
                    main_io.getReadPointer(channel), already_processed_samples);
                if (num_windows_processed_ >= windowing_overlap_times) {
                    output_ring_buffers_[channel].copy_n_to(
                        main_io.getWritePointer(channel),
                        already_processed_samples, true);
                } else {
                    main_io.clear(channel, 0, already_processed_samples);
                }
                if constexpr (sidechain_active) {
                    sidechain_ring_buffers_[channel].read_n_from(
                        sidechain_io.getReadPointer(channel),
                        already_processed_samples);
                }
            }

            sample_buffer_offset += already_processed_samples;
        }

        // Now if `windows_to_process > 0`, the current ring buffer position
        // will align with a window and we can start doing our FFT magic
        for (int window_idx = 0; window_idx < windows_to_process;
             window_idx++) {
            if constexpr (sidechain_active && !bypassed) {
                // The sidechain input is only used for analysis
                for (size_t channel = 0; channel < num_channels; channel++) {
                    sidechain_ring_buffers_[channel].copy_last_n_to(
                        fft_scratch_buffer_.data(), fft_window_size);
                    windowing_function_.multiplyWithWindowingTable(
                        fft_scratch_buffer_.data(), fft_window_size);
                    // TODO: We can skip negative frequencies here, right?
                    fft_.performRealOnlyForwardTransform(
                        fft_scratch_buffer_.data(), true);

                    const std::span<std::complex<float>> fft_buffer(
                        reinterpret_cast<std::complex<float>*>(
                            fft_scratch_buffer_.data()),
                        fft_window_size);
                    sidechain_fn(fft_buffer, channel);
                }

                // The user might want to do some aggregation after processing
                // every channel
                post_sidechain_fn();
            }

            // This is where the magic happens!
            for (size_t channel = 0; channel < num_channels; channel++) {
                if constexpr (!bypassed) {
                    // Depending on what stage of the transformation process
                    // we're in, our scratch buffer will contain either samples
                    // or complex frequency bins. The caller should get a chance
                    // to preprocess the (windowed) samples, process the
                    // transformed data, and the postprocess the results after
                    // the windowing function has been applied after the inverse
                    // transformation.
                    std::span<float> sample_buffer(fft_scratch_buffer_.data(),
                                                   fft_window_size);
                    std::span<std::complex<float>> fft_buffer(
                        reinterpret_cast<std::complex<float>*>(
                            fft_scratch_buffer_.data()),
                        fft_window_size);

                    input_ring_buffers_[channel].copy_last_n_to(
                        fft_scratch_buffer_.data(), fft_window_size);
                    windowing_function_.multiplyWithWindowingTable(
                        fft_scratch_buffer_.data(), fft_window_size);
                    preprocess_fn(sample_buffer, channel);

                    fft_.performRealOnlyForwardTransform(
                        fft_scratch_buffer_.data());
                    process_fn(fft_buffer, channel);

                    fft_.performRealOnlyInverseTransform(
                        fft_scratch_buffer_.data());
                    windowing_function_.multiplyWithWindowingTable(
                        fft_scratch_buffer_.data(), fft_window_size);
                    postprocess_fn(sample_buffer, channel);

                    // After processing the windowed data, we'll add it to our
                    // output ring buffer with any (automatic) makeup gain
                    // applied
                    output_ring_buffers_[channel].add_n_from_in_place(
                        fft_scratch_buffer_.data(), fft_window_size, gain);
                } else {
                    // TODO: Implement the bypass to copy directly between the
                    //       ring buffers instead of going through the scratch
                    //       buffer
                    input_ring_buffers_[channel].copy_last_n_to(
                        fft_scratch_buffer_.data(), windowing_interval);
                    output_ring_buffers_[channel].read_n_from_in_place(
                        fft_scratch_buffer_.data(), windowing_interval);
                }
            }

            // We don't copy over anything to the outputs until we processed a
            // full buffer
            num_windows_processed_ += 1;

            // Copy the input audio into our ring buffer and copy the processed
            // audio into the output buffer
            const size_t samples_to_process_this_iteration = std::min(
                windowing_interval, num_samples - sample_buffer_offset);
            for (size_t channel = 0; channel < num_channels; channel++) {
                input_ring_buffers_[channel].read_n_from(
                    main_io.getReadPointer(channel) + sample_buffer_offset,
                    samples_to_process_this_iteration);
                if (num_windows_processed_ >= windowing_overlap_times) {
                    output_ring_buffers_[channel].copy_n_to(
                        main_io.getWritePointer(channel) + sample_buffer_offset,
                        samples_to_process_this_iteration, true);
                } else {
                    main_io.clear(channel, sample_buffer_offset,
                                  samples_to_process_this_iteration);
                }
                if constexpr (sidechain_active) {
                    sidechain_ring_buffers_[channel].read_n_from(
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
     * The numbers of windows already processed. We use this to reduce clicks by
     * not copying over audio to the output during the first
     * `windowing_overlap_times` windows.
     */
    int num_windows_processed_ = 0;

    /**
     * The FFT processor.
     */
    juce::dsp::FFT fft_;

    /**
     * We'll process the signal with overlapping windows that are added to each
     * other to form the output signal. See `input_ring_buffers` for more
     * information on how we'll do this.
     */
    juce::dsp::WindowingFunction<float> windowing_function_;

    /**
     * We need a scratch buffer that can contain `fft_window_size * 2` samples
     * for `fft` to work in.
     */
    std::vector<float> fft_scratch_buffer_;

    /**
     * A ring buffer of size `fft_window_size` for every channel. Every
     * `windowing_interval` we'll copy the last `fft_window_size` samples to
     * `fft_scratch_buffers` using a window function, process it, and then add
     * the results to `output_ring_buffers`.
     */
    std::vector<RingBuffer<float>> input_ring_buffers_;
    /**
     * These ring buffers are identical to `input_ring_buffers`, but with data
     * from the sidechain input. When sidechaining is enabled, we set the
     * compressor thresholds based on the magnitudes from the same FFT analysis
     * applied to the sidechain input.
     */
    std::vector<RingBuffer<float>> sidechain_ring_buffers_;
    /**
     * The processed results as described in the docstring of
     * `input_ring_buffers`. Samples from this buffer will be written to the
     * output.
     */
    std::vector<RingBuffer<float>> output_ring_buffers_;
};
