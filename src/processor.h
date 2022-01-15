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

#include <optional>

#include <juce_audio_processors/juce_audio_processors.h>
#include <juce_dsp/juce_dsp.h>

#include "dsp/compressor.h"
#include "dsp/stft.h"
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
     * This is where the magic happens. Performs the entire STFT and overlap-add
     * process for us. See the `STFT` class for more information.
     */
    std::optional<STFT<true>> stft;

    /**
     * This will contain `fft_window_size / 2` compressors. The compressors are
     * already multichannel so we don't need a nested vector here. We'll
     * compress the magnitude of every FFT bin (`sqrt(i^2 + r^2)`) individually,
     * and then scale both the real and imaginary components by the ratio of
     * their magnitude and the compressed value. Bin 0 is the DC offset and the
     * bins in the second half should be processed the same was as the bins in
     * the first half but mirrored.
     */
    std::vector<MultiwayCompressor<float>> spectral_compressors;

    /**
     * When setting compressor thresholds based on a sidechain signal we should
     * be taking the average bin magnitudes of all channels. This buffer
     * accumulates `spectral_compressors.size()` threshold values while
     * iterating over the channels of the sidechain signal so we can then
     * average them and configure the compressors based on that.
     */
    std::vector<float> spectral_compressor_sidechain_thresholds;
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
     * (Re)initialize a process data object and all compressors within it for
     * the current FFT order on the next audio processing cycle. The inactive
     * object we're modifying will be swapped with the active object on the next
     * call to `process_data.get()`. This should not be called from the audio
     * thread.
     */
    void update_and_swap_process_data();

    /**
     * This contains all of our scratch buffers, ring buffers, compressors, and
     * everything else that depends on the FFT window size.
     */
    AtomicallySwappable<ProcessData> process_data_;

    /**
     * A dry-wet mixer we'll use to be able to blend the processed and the
     * unprocessed signals.
     */
    juce::dsp::DryWetMixer<float> mixer_;

    /**
     * Will be set during `prepareToPlay()`, needed to initialize compressors
     * when resizing our buffers.
     */
    juce::uint32 max_samples_per_block_ = 0;
    /**
     * The 'effective sample rate' (sample rate divided by the windowing
     * interval) for the last processing cycle. If this changes, then we'll need
     * to adjust our compressors accordingly.
     */
    double last_effective_sample_rate_ = 0.0;

    juce::AudioProcessorValueTreeState parameters_;

    /**
     * This is applied after windowing, just before the forward FFT
     * transformation.
     */
    std::atomic<float>& input_gain_db_;
    /**
     * This is essentially the makeup gain, in dB. When automatic makeup gain is
     * enabled this is added on top of that.
     */
    std::atomic<float>& output_gain_db_;
    /**
     * Try to automatically compensate for low thresholds. Doesn't do anything
     * when sidechaining is active.
     */
    juce::AudioParameterBool& auto_makeup_gain_;
    /**
     * Set the DC bin to 0 when enabled. This tends to otherwise just make the
     * signal a lot louder.
     */
    juce::AudioParameterBool& dc_filter_;
    /**
     * How much of the dry signal to mix in with the processed signal. This
     * mixing is done after applying the output gain.
     */
    std::atomic<float>& dry_wet_ratio_;

    /**
     * If true, set the compressor thresholds based on the sidechain signal.
     */
    juce::AudioParameterBool& sidechain_active_;
    /**
     * If set to true, the compressor will be configured based on raw
     * magnitudes, even though it expects decibels. This is incorrect, but it
     * does result in some interesting sound design possibilities.
     *
     * TODO: Decide on if we want to keep this
     */
    juce::AudioParameterBool& sidechain_exponential_;
    // TODO: Add threshold offsets, right now you need to do the oldschool
    //       compressor thing of simultaneously modifying the input and output
    //       gains
    /**
     * The compressor's mode. Either downwards, upwards, or simultaneous upwards
     * and downwards compression.
     */
    juce::AudioParameterChoice& compressor_mode_;
    /**
     * The range in dB around the threshold in decibel where the compressor
     * should not be active when the compressor is set to multiway mode.
     */
    std::atomic<float>& compressor_multiway_deadzone_;
    /**
     * Compressor ratio, where everything above 1.0 means that the signal will
     * be compressed above the threshold.
     */
    std::atomic<float>& compressor_ratio_;
    /**
     * Compressor release time in milliseconds.
     */
    std::atomic<float>& compressor_attack_ms_;
    /**
     * Compressor attack time in milliseconds.
     */
    std::atomic<float>& compressor_release_ms_;
    /**
     * Will cause the compressor settings to be updated on the next processing
     * cycle whenever a compressor parameter changes.
     */
    LambdaParameterListener compressor_settings_listener_;
    /**
     * Will be set in `CompressorSettingsListener` when any of the compressor
     * related settings change so we can update our compressors. We'll
     * initialize this to true so the compressors will be initialized during the
     * first processing cycle.
     */
    std::atomic_bool compressor_settings_changed_ = true;

    /**
     * The order (where `fft_window_size = 1 << fft_order`) for our spectral
     * operations. When this gets changed, we'll resize all of our buffers and
     * atomically swap the current and the resized buffers.
     */
    juce::AudioParameterInt& fft_order_;
    /**
     * The order of the overlap for the windowing (where
     * `windowing_overlap_times = 1 1 << windowing_overlap_order`). We end up
     * processing the signal in `fft_window_size` windows every `fft_window_size
     * / windowing_overlap_times` samples. When this setting gets changed, we'll
     * also have to update our compressors since the effective sample rate also
     * changes.
     */
    juce::AudioParameterInt& windowing_overlap_order_;
    /**
     * Atomically resizes the object `ProcessData` from a background thread.
     */
    LambdaAsyncUpdater process_data_updater_;
    /**
     * When the FFT order parameter changes, we'll have to create a new
     * `ProcessData` object for the new FFT window size (or rather, resize an
     * inactive one to match the new size).
     */
    LambdaParameterListener fft_order_listener_;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpectralCompressorProcessor)
};
