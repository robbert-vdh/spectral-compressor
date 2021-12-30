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

#include <juce_dsp/juce_dsp.h>

/**
 * A compressor that simultaneously performs both upwards and downwards
 * compression. Based on `juce::dsp::Compressor`.
 */
template <typename T>
class MultiwayCompressor {
   public:
    MultiwayCompressor() { update(); }

    /**
     * Modes for downwards, upwards, or simulteneous upwards and downwards
     * compression. In the last mode the multiway deadzone parameter acts as an
     * bidirectional offset to the threshold where the compressor doesn't do
     * anything.
     */
    enum class Mode { downwards, upwards, multiway };

    /**
     * Set the compressor's mode.
     */
    void set_mode(Mode mode) { mode_ = mode; }

    /**
     * Set the compressor's deadzone when using the multiway mode. Must not be
     * negative.
     */
    void set_multiway_deadzone(T deadzone_db) {
        jassert(deadzone_db >= 0);

        multiway_deadzone_db_ = deadzone_db;
        update();
    }

    /**
     * Set the compressor's threshold in dB.
     */
    void set_threshold(T threshold_db) {
        threshold_db_ = threshold_db;
        update();
    }

    /**
     * Set the compressor's ratio (must be higher or equal to 1).
     */
    void set_ratio(T ratio) {
        jassert(ratio >= 1);

        ratio_ = ratio;
        update();
    }

    /**
     * Set the compressor's attack time in milliseconds.
     */
    void set_attack(T attack) {
        jassert(attack >= 0);

        attack_time_ = attack;
        update();
    }

    /**
     * Set the compressor's release time in milliseconds.
     */
    void set_release(T release) {
        jassert(release >= 0);

        release_time_ = release;
        update();
    }

    /**
     * Initialize the processor.
     */
    void prepare(const juce::dsp::ProcessSpec& spec) {
        jassert(spec.sampleRate > 0);
        jassert(spec.numChannels > 0);

        sample_rate_ = spec.sampleRate;
        envelope_filter_.prepare(spec);

        update();
        reset();
    }

    /**
     * Reset the internal state variables of the processor.
     */
    void reset() { envelope_filter_.reset(); }

    /**
     * Process the input and output samples supplied in the processing context.
     */
    template <typename ProcessContext>
    void process(const ProcessContext& context) noexcept {
        const auto& input_block = context.getInputBlock();
        auto& output_block = context.getOutputBlock();
        const auto num_channels = output_block.getNumChannels();
        const auto num_samples = output_block.getNumSamples();

        jassert(input_block.getNumChannels() == num_channels);
        jassert(input_block.getNumSamples() == num_samples);

        if (context.isBypassed) {
            output_block.copyFrom(input_block);
            return;
        }

        for (size_t channel = 0; channel < num_channels; ++channel) {
            auto* input_samples = input_block.getChannelPointer(channel);
            auto* output_samples = output_block.getChannelPointer(channel);

            for (size_t i = 0; i < num_samples; ++i)
                output_samples[i] = process_sample(channel, input_samples[i]);
        }
    }

    /**
     * Process a single sample.
     */
    T process_sample(int channel, T input) {
        // Ballistics filter with peak rectifier
        auto env = envelope_filter_.processSample(channel, input);

        // VCA
        auto gain = (env < threshold_)
                        ? static_cast<T>(1.0)
                        : std::pow(env * threshold_inverse_,
                                   ratio_inverse_ - static_cast<T>(1.0));

        return input * gain;
    }

   private:
    void update() {
        // The deadzone acts in both directions, so it needs to be divided by
        // two
        multiway_deadzone_ =
            std::abs(static_cast<T>(1.0) -
                     juce::Decibels::decibelsToGain(multiway_deadzone_db_)) /
            static_cast<T>(2.0);
        threshold_ = juce::Decibels::decibelsToGain(threshold_db_,
                                                    static_cast<T>(-200.0));
        threshold_inverse_ = static_cast<T>(1.0) / threshold_;
        ratio_inverse_ = static_cast<T>(1.0) / ratio_;

        envelope_filter_.setAttackTime(attack_time_);
        envelope_filter_.setReleaseTime(release_time_);
    }

    Mode mode_ = Mode::downwards;
    double sample_rate_ = 44100.0;
    T threshold_db_ = 0.0;
    T multiway_deadzone_db_ = 0.0;
    T ratio_ = 1.0;
    T attack_time_ = 1.0;
    T release_time_ = 100.0;

    T threshold_ = 1.0;
    T threshold_inverse_ = 1.0;
    T multiway_deadzone_ = 0.0;
    T ratio_inverse_ = 1.0;
    juce::dsp::BallisticsFilter<T> envelope_filter_;
};
