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
     * The lowest sample value we will try to upwards compress, otherwise we
     * could get infinite gain ratios and it would be waste to try to make
     * silence louder.
     */
    static constexpr T epsilon =
        // GCC supports constexpr `std::pow`, MSVC and Clang don't
        // std::pow(static_cast<T>(10.0), static_cast<T>(-100 * 0.05));
        1e-05;

    /**
     * The maximum gain value, to minimize ear damage when doing upwards
     * compression.
     */
    static constexpr T gain_limit = 200;

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
        T env = envelope_filter_.processSample(channel, input);

        // The VCA can push the gain either upwards, downwards, or do nothing
        // depending on the settings and the compressor's mode
        // TODO: This can be optimized a bit, but it's fine for now
        T gain = 1.0;
        if (mode_ != Mode::upwards && env > (threshold_ + multiway_deadzone_)) {
            // Downwards compression
            gain = std::pow((env - multiway_deadzone_) * threshold_inverse_,
                            ratio_inverse_ - static_cast<T>(1.0));
        } else if (mode_ != Mode::downwards && env > epsilon &&
                   env < (threshold_ - multiway_deadzone_)) {
            // Upwards compression
            gain = std::pow((env + multiway_deadzone_) * threshold_inverse_,
                            ratio_inverse_ - static_cast<T>(1.0));

            // When levels drop very low crazy things start happening. At that
            // point it's best to just cap the gain ratio.
            if (gain > gain_limit) {
                gain = gain_limit;
            }
        }

        return input * gain;
    }

   private:
    void update() {
        // The deadzone acts in both directions, so it needs to be divided by
        // two. If multiway mode is not enabled then the deadzone is always set
        // to 0 to simplify the calculations.
        multiway_deadzone_ =
            mode_ == Mode::multiway
                ? std::abs(static_cast<T>(1.0) - juce::Decibels::decibelsToGain(
                                                     multiway_deadzone_db_)) /
                      static_cast<T>(2.0)
                : 0.0;
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
    /**
     * Will always be set to 0 when the mode is not set to multiway regardless
     * of the value of `multiway_deadzone_db_`.
     */
    T multiway_deadzone_ = 0.0;
    T ratio_inverse_ = 1.0;
    juce::dsp::BallisticsFilter<T> envelope_filter_;
};
