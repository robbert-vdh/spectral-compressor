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

using juce::uint32;

constexpr char input_gain_db_param_name[] = "input_gain";
constexpr char output_gain_db_param_name[] = "output_gain";
constexpr char dry_wet_ratio_param_name[] = "mix";
constexpr char auto_makeup_gain_param_name[] = "auto_makeup_gain";

constexpr char compressor_settings_group_name[] = "compressors";
constexpr char sidechain_active_param_name[] = "sidechain_active";
constexpr char compressor_ratio_param_name[] = "compressor_ratio";
constexpr char compressor_attack_ms_param_name[] = "compressor_attack";
constexpr char compressor_release_ms_param_name[] = "compressor_release";

constexpr char spectral_settings_group_name[] = "spectral";
constexpr char fft_order_param_name[] = "fft_size";
constexpr char windowing_overlap_order_param_name[] = "windowing_order";

constexpr int fft_order_minimum = 12;
constexpr int fft_order_maximum = 15;

SpectralCompressorProcessor::SpectralCompressorProcessor()
    : AudioProcessor(
          BusesProperties()
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)
              .withInput("Sidechain", juce::AudioChannelSet::stereo(), true)),
      mixer(1 << fft_order_maximum),
      parameters(
          *this,
          nullptr,
          "parameters",
          {
              std::make_unique<juce::AudioProcessorParameterGroup>(
                  compressor_settings_group_name,
                  "Master",
                  " | ",
                  std::make_unique<juce::AudioParameterFloat>(
                      input_gain_db_param_name,
                      "Input Gain",
                      juce::NormalisableRange<float>(-50, 50, 0.1),
                      0,
                      " dB"),
                  std::make_unique<juce::AudioParameterFloat>(
                      output_gain_db_param_name,
                      "Output Gain",
                      juce::NormalisableRange<float>(-50, 50, 0.1),
                      0,
                      " dB"),
                  std::make_unique<juce::AudioParameterBool>(
                      auto_makeup_gain_param_name,
                      "Auto Makeup Gain",
                      true),
                  std::make_unique<juce::AudioParameterFloat>(
                      dry_wet_ratio_param_name,
                      "Mix",
                      juce::NormalisableRange<float>(0.0, 1.0, 0.01),
                      1.0,
                      "%",
                      juce::AudioProcessorParameter::genericParameter,
                      [](float value, int /*max_length*/) -> juce::String {
                          return juce::String(value * 100.0f, 0);
                      },
                      [](const juce::String& text) -> float {
                          return text.getFloatValue() / 100.0f;
                      })),
              std::make_unique<juce::AudioProcessorParameterGroup>(
                  compressor_settings_group_name,
                  "Compressors",
                  " | ",
                  std::make_unique<juce::AudioParameterBool>(
                      sidechain_active_param_name,
                      "Sidechain Active",
                      false),
                  std::make_unique<juce::AudioParameterFloat>(
                      compressor_ratio_param_name,
                      "Ratio",
                      juce::NormalisableRange<float>(1.0, 300.0, 0.1, 0.25),
                      50.0),
                  std::make_unique<juce::AudioParameterFloat>(
                      compressor_attack_ms_param_name,
                      "Attack",
                      juce::NormalisableRange<float>(0.0, 10000.0, 1.0, 0.2),
                      50.0,
                      " ms",
                      juce::AudioProcessorParameter::genericParameter),
                  std::make_unique<juce::AudioParameterFloat>(
                      compressor_release_ms_param_name,
                      "Release",
                      juce::NormalisableRange<float>(0.0, 10000.0, 1.0, 0.2),
                      5000.0,
                      " ms",
                      juce::AudioProcessorParameter::genericParameter)),
              std::make_unique<juce::AudioProcessorParameterGroup>(
                  spectral_settings_group_name,
                  "Spectral Settings",
                  " | ",
                  std::make_unique<juce::AudioParameterInt>(
                      fft_order_param_name,
                      "Resolution",
                      9,
                      fft_order_maximum,
                      fft_order_minimum,
                      "",
                      [](int value, int /*max_length*/) -> juce::String {
                          return juce::String(1 << value);
                      },
                      [](const juce::String& text) -> int {
                          return std::log2(text.getIntValue());
                      }),
                  std::make_unique<juce::AudioParameterInt>(
                      windowing_overlap_order_param_name,
                      "Overlap",
                      2,
                      6,
                      2,
                      "x",
                      [&](int value, int /*max_length*/) -> juce::String {
                          return juce::String(1 << value);
                      },
                      [&](const juce::String& text) -> int {
                          return std::log2(text.getIntValue());
                      })),
          }),
      // TODO: Is this how you're supposed to retrieve non-float parameters?
      //       Seems a bit excessive
      input_gain_db(*parameters.getRawParameterValue(input_gain_db_param_name)),
      output_gain_db(
          *parameters.getRawParameterValue(output_gain_db_param_name)),
      auto_makeup_gain(*dynamic_cast<juce::AudioParameterBool*>(
          parameters.getParameter(auto_makeup_gain_param_name))),
      dry_wet_ratio(*parameters.getRawParameterValue(dry_wet_ratio_param_name)),
      sidechain_active(*dynamic_cast<juce::AudioParameterBool*>(
          parameters.getParameter(sidechain_active_param_name))),
      compressor_ratio(
          *parameters.getRawParameterValue(compressor_ratio_param_name)),
      compressor_attack_ms(
          *parameters.getRawParameterValue(compressor_attack_ms_param_name)),
      compressor_release_ms(
          *parameters.getRawParameterValue(compressor_release_ms_param_name)),
      compressor_settings_listener(
          [&](const juce::String& /*parameterID*/, float /*newValue*/) {
              compressor_settings_changed = true;
          }),
      fft_order(*dynamic_cast<juce::AudioParameterInt*>(
          parameters.getParameter(fft_order_param_name))),
      windowing_overlap_order(*dynamic_cast<juce::AudioParameterInt*>(
          parameters.getParameter(windowing_overlap_order_param_name))),
      process_data_updater([&]() {
          update_and_swap_process_data();

          const size_t new_window_size = 1 << fft_order;
          setLatencySamples(new_window_size);
      }),
      fft_order_listener(
          [&](const juce::String& /*parameterID*/, float /*newValue*/) {
              process_data_updater.triggerAsyncUpdate();
          }) {
    // TODO: Move the latency computation elsewhere
    const size_t new_window_size = 1 << fft_order;
    setLatencySamples(new_window_size);

    // XXX: There doesn't seem to be a fool proof way to just iterate over all
    //      parameters in a group, right?
    for (const auto& compressor_param_name :
         {sidechain_active_param_name, compressor_ratio_param_name,
          compressor_attack_ms_param_name, compressor_release_ms_param_name}) {
        parameters.addParameterListener(compressor_param_name,
                                        &compressor_settings_listener);
    }

    parameters.addParameterListener(fft_order_param_name, &fft_order_listener);
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
    max_samples_per_block = static_cast<uint32>(maximumExpectedSamplesPerBlock);

    // This is used to set the correct 'effective' sample rate on our
    // compressors during the processing loop
    last_effective_sample_rate = 0.0;

    // TODO: We may be doing double work here when `process_data_updater`
    //       changes the latency and the host restarts playback
    // After initializing the process data we make an explicit call to
    // `process_data.get()` to swap the two filters in case we get a
    // parameter change before the first processing cycle
    update_and_swap_process_data();
    process_data.get();

    mixer.prepare(juce::dsp::ProcessSpec{
        .sampleRate = sampleRate,
        .maximumBlockSize = static_cast<uint32>(maximumExpectedSamplesPerBlock),
        .numChannels = static_cast<uint32>(getMainBusNumInputChannels())});
}

void SpectralCompressorProcessor::releaseResources() {
    process_data.clear([](ProcessData& process_data) {
        process_data.stft.reset();

        process_data.spectral_compressors.clear();
        process_data.spectral_compressors.shrink_to_fit();
        process_data.spectral_compressor_sidechain_thresholds.clear();
        process_data.spectral_compressor_sidechain_thresholds.shrink_to_fit();
    });
    mixer.reset();
}

bool SpectralCompressorProcessor::isBusesLayoutSupported(
    const BusesLayout& layouts) const {
    // We can support any number of channels, as long as the main input, main
    // output, and sidechain input have the same number of channels
    const juce::AudioChannelSet sidechain_channel_set =
        layouts.getChannelSet(true, 1);
    return (layouts.getMainInputChannelSet() ==
            layouts.getMainOutputChannelSet()) &&
           (sidechain_channel_set == layouts.getMainInputChannelSet()) &&
           !layouts.getMainInputChannelSet().isDisabled();
}

void SpectralCompressorProcessor::processBlockBypassed(
    juce::AudioBuffer<float>& buffer,
    juce::MidiBuffer& /*midiMessages*/) {
    juce::AudioBuffer<float> main_io = getBusBuffer(buffer, true, 0);

    // We need to maintain the same latency when bypassed, so we'll reuse most
    // of the processing logic
    ProcessData& process_data = this->process_data.get();
    process_data.stft->process_bypassed(main_io);
}

void SpectralCompressorProcessor::processBlock(
    juce::AudioBuffer<float>& buffer,
    juce::MidiBuffer& /*midiMessages*/) {
    juce::AudioBuffer<float> main_io = getBusBuffer(buffer, true, 0);
    juce::AudioBuffer<float> sidechain_io = getBusBuffer(buffer, true, 1);

    juce::dsp::AudioBlock<float> main_block(main_io);
    mixer.setWetMixProportion(dry_wet_ratio);
    mixer.pushDrySamples(main_block);

    ProcessData& process_data = this->process_data.get();
    const double effective_sample_rate =
        getSampleRate() /
        (static_cast<double>(process_data.stft->fft_window_size) /
         (1 << windowing_overlap_order));
    const float fft_frequency_increment =
        getSampleRate() / process_data.stft->fft_window_size;

    // We have two different gain stages: just before the FFT transformations,
    // after the FFT transformations (the makeup gain). As part of the makeup
    // gain we also compensate for the overlap caused by our windowing. We don't
    // need any manual ramps or fades here because that's already included in
    // our Hanning windows.
    // TODO: We should probably also compensate for different FFT window sizes
    const float input_gain =
        juce::Decibels::decibelsToGain(static_cast<float>(input_gain_db));
    float makeup_gain =
        (1.0f / (1 << windowing_overlap_order)) +
        juce::Decibels::decibelsToGain(static_cast<float>(output_gain_db));
    if (auto_makeup_gain) {
        if (sidechain_active) {
            // Not really sure what makes sense here
            // TODO: Take base threshold into account
            makeup_gain *= (compressor_ratio + 24.0f) / 25.0f;
        } else {
            // TODO: Make this smarter, make it take all of the compressor
            //       parameters into account. It will probably start making
            //       sense once we add parameters for the threshold and ratio.
            // FIXME: This makes zero sense! But it works for our current
            //        parameters.
            makeup_gain *=
                (std::log10(compressor_ratio * 100.00f) * 200.0f) - 399.0f;
        }
    }

    auto preprocess_fn = [input_gain](std::span<float>& samples,
                                      size_t /*channel*/) {
        // We apply the input gain after the windowing, just before the forward
        // FFT transformation
        juce::FloatVectorOperations::multiply(samples.data(), input_gain,
                                              samples.size());
    };

    auto process_fn = [this, effective_sample_rate, fft_frequency_increment,
                       &process_data](std::span<std::complex<float>>& fft,
                                      size_t channel) {
        // We'll update the compressor settings just before processing if the
        // settings have changed or if the sidechaining has been disabled
        bool expected = true;
        const bool update_compressors_now =
            compressor_settings_changed.compare_exchange_weak(expected, false);

        // If any timing related settings change (so the FFT window size or the
        // amount of overlap), we'll need to adjust our compressors accordingly.
        // Since this process can cause pops and clicks, we only do it when
        // necessary.
        const bool update_sample_rate_now =
            last_effective_sample_rate != effective_sample_rate;
        last_effective_sample_rate = effective_sample_rate;

        // We'll compress every FTT bin individually. Bin 0 is the DC offset and
        // should be skipped, and the latter half of the FFT bins should be
        // processed in the same way as the first half but in reverse order. The
        // real and imaginary parts are interleaved, so ever bin spans two
        // values in the scratch buffer. We can 'safely' do this cast so we can
        // use the STL's complex value functions.
        // TODO: It might be nice to add a DC filter, which would be
        //       very cheap since we're already doing FFT anyways
        for (size_t compressor_idx = 0;
             compressor_idx < process_data.spectral_compressors.size();
             compressor_idx++) {
            auto& compressor =
                process_data.spectral_compressors[compressor_idx];
            // We don't have a compressor for the first bin
            const size_t bin_idx = compressor_idx + 1;

            if (update_compressors_now) {
                compressor.setRatio(compressor_ratio);
                compressor.setAttack(compressor_attack_ms);
                compressor.setRelease(compressor_release_ms);
                // TODO: The user should be able to configure their own slope
                //       (or free drawn)
                // TODO: Change the calculations so that the base threshold
                //       parameter is centered around some frequency
                // TODO: And we should be doing both upwards and downwards
                //       compression, OTT-style
                if (!sidechain_active) {
                    constexpr float base_threshold_dbfs = 0.0f;
                    const float frequency = fft_frequency_increment * bin_idx;

                    // This starts at 1 for 0 Hz (DC)
                    const float octave = std::log2(frequency + 2);

                    // The 3 dB is to compensate for bin 0
                    compressor.setThreshold((base_threshold_dbfs + 3.0f) -
                                            (3.0f * octave));
                }
            }

            if (update_sample_rate_now) {
                // TODO: This prepare resets the envelope follower, which is not
                //       what we want. In our own compressor we should have a
                //       way to just change the sample rate.
                // TODO: Now that the timings are compensated for changing
                //       window intervals, we might not need this to be
                //       configurable anymore can just leave this fixed at 4x.
                compressor.prepare(juce::dsp::ProcessSpec{
                    // We only process everything once every
                    // `windowing_interval`, otherwise our attack and release
                    // times will be all messed up
                    .sampleRate = effective_sample_rate,
                    .maximumBlockSize = max_samples_per_block,
                    .numChannels =
                        static_cast<uint32>(getMainBusNumInputChannels())});
            }

            const float magnitude = std::abs(fft[bin_idx]);
            const float compressed_magnitude =
                compressor.processSample(channel, magnitude);

            // We need to scale both the imaginary and real components of the
            // bins at the start and end of the spectrum by the same value
            // TODO: Add stereo linking
            const float compression_multiplier =
                magnitude != 0.0f ? compressed_magnitude / magnitude : 1.0f;

            // Since we're usign the real-only FFT operations we don't need to
            // touch the second, mirrored half of the FFT bins
            fft[bin_idx] *= compression_multiplier;
        }

        // TODO: We might need some kind of optional limiting stage to
        //       be safe
        // TODO: We should definitely add a way to recover transients
        //       from the original input audio, that sounds really good
    };

    auto postprocess_fn = [](std::span<float>& /*samples*/,
                             size_t /*channel*/) {};

    // We'll process the input signal in windows, using overlap-add
    if (sidechain_active) {
        process_data.stft->process(
            main_io, sidechain_io, 1 << windowing_overlap_order, makeup_gain,
            [&process_data](const std::span<std::complex<float>>& fft,
                            size_t /*channel*/) {
                // If sidechaining is active, we set the compressor thresholds
                // based on a sidechain signal. Since compression is already
                // ballistics based we don't need any additional smoothing when
                // updating those thresholds.
                for (size_t compressor_idx = 0;
                     compressor_idx < process_data.spectral_compressors.size();
                     compressor_idx++) {
                    const size_t bin_idx = compressor_idx + 1;
                    const float magnitude = std::abs(fft[bin_idx]);

                    // We'll set the compressor threshold based on the
                    // arithmetic mean of the magnitudes of all channels. As
                    // a slight premature optimization (sorry) we'll reset
                    // these magnitudes after using them to avoid the
                    // conditional here.
                    process_data.spectral_compressor_sidechain_thresholds
                        [compressor_idx] += magnitude;
                }
            },
            [&process_data, num_channels = sidechain_io.getNumChannels()]() {
                // After adding up the magnitudes for each bin in
                // `process_data.spectral_compressor_sidechain_thresholds` we
                // want to actually configure the compressor thresholds based on
                // the mean across the different channels
                for (size_t compressor_idx = 0;
                     compressor_idx < process_data.spectral_compressors.size();
                     compressor_idx++) {
                    process_data.spectral_compressors[compressor_idx]
                        .setThreshold(
                            process_data
                                .spectral_compressor_sidechain_thresholds
                                    [compressor_idx] /
                            num_channels);
                    process_data.spectral_compressor_sidechain_thresholds
                        [compressor_idx] = 0;
                }
            },
            preprocess_fn, process_fn, postprocess_fn);
    } else {
        process_data.stft->process(main_io, 1 << windowing_overlap_order,
                                   makeup_gain, preprocess_fn, process_fn,
                                   postprocess_fn);
    }

    mixer.setWetLatency(process_data.stft->latency_samples());
    mixer.mixWetSamples(main_block);
}

bool SpectralCompressorProcessor::hasEditor() const {
    return true;
}

juce::AudioProcessorEditor* SpectralCompressorProcessor::createEditor() {
    // TODO: Add an editor at some point
    // return new SpectralCompressorEditor(*this);
    return new juce::GenericAudioProcessorEditor(*this);
}

void SpectralCompressorProcessor::getStateInformation(
    juce::MemoryBlock& destData) {
    const std::unique_ptr<juce::XmlElement> xml =
        parameters.copyState().createXml();
    copyXmlToBinary(*xml, destData);
}

void SpectralCompressorProcessor::setStateInformation(const void* data,
                                                      int sizeInBytes) {
    const auto xml = getXmlFromBinary(data, sizeInBytes);
    if (xml && xml->hasTagName(parameters.state.getType())) {
        parameters.replaceState(juce::ValueTree::fromXml(*xml));
    }

    // TODO: Should we do this here, is will `prepareToPlay()` always be called
    //       between loading presets and audio processing starting?
    update_and_swap_process_data();

    // TODO: Do parameter listeners get triggered? Or alternatively, can this be
    //       called during playback (without `prepareToPlay()` being called
    //       first)?
    // TODO: Move the latency computation elsewhere
    const size_t new_window_size = 1 << fft_order;
    setLatencySamples(new_window_size);
}

void SpectralCompressorProcessor::update_and_swap_process_data() {
    process_data.modify_and_swap([this](ProcessData& process_data) {
        process_data.stft.emplace(getMainBusNumInputChannels(), fft_order);

        // Every FFT bin on both channels gets its own compressor, hooray! The
        // `fft_window_size / 2` is because the first bin is the DC offset and
        // shouldn't be compressed, and the bins after the Nyquist frequency are
        // the same as the first half but in reverse order. The compressor
        // settings will be set in `update_compressors()`, which is triggered on
        // the next processing cycle by setting `compressor_settings_changed`
        // below.
        process_data.spectral_compressors.resize(
            process_data.stft->fft_window_size / 2);
        process_data.spectral_compressor_sidechain_thresholds.resize(
            process_data.spectral_compressors.size());

        // After resizing the compressors are uninitialized and should be
        // reinitialized
        compressor_settings_changed = true;
        last_effective_sample_rate = 0.0;
    });
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
    return new SpectralCompressorProcessor();
}
