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

constexpr char compressor_settings_group_name[] = "compressors";
constexpr char sidechain_active_param_name[] = "sidechain_active";
constexpr char compressor_ratio_param_name[] = "compressor_ratio";
constexpr char auto_makeup_gain_param_name[] = "auto_makeup_gain";

constexpr char spectral_settings_group_name[] = "spectral";
constexpr char fft_order_param_name[] = "fft_size";
constexpr char windowing_overlap_times_param_name[] = "windowing_times";

LambdaAsyncUpdater::LambdaAsyncUpdater(fu2::unique_function<void()> callback)
    : callback(std::move(callback)) {}

void LambdaAsyncUpdater::handleAsyncUpdate() {
    callback();
}

LambdaParameterListener::LambdaParameterListener(
    fu2::unique_function<void(const juce::String&, float)> callback)
    : callback(std::move(callback)) {}

void LambdaParameterListener::parameterChanged(const juce::String& parameterID,
                                               float newValue) {
    callback(parameterID, newValue);
}

SpectralCompressorProcessor::SpectralCompressorProcessor()
    : AudioProcessor(
          BusesProperties()
              .withInput("Input", juce::AudioChannelSet::stereo(), true)
              .withOutput("Output", juce::AudioChannelSet::stereo(), true)
              .withInput("Sidechain", juce::AudioChannelSet::stereo(), true)),
      process_data([&](ProcessData& inactive, size_t new_fft_order) {
          initialize_process_data(inactive, new_fft_order);
      }),
      parameters(
          *this,
          nullptr,
          "parameters",
          {
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
                      "Compressor Ratio",
                      juce::NormalisableRange<float>(1.0, 300.0, 0.1, 0.25),
                      50.0),
                  std::make_unique<juce::AudioParameterBool>(
                      auto_makeup_gain_param_name,
                      "Auto Makeup Gain",
                      true)),
              std::make_unique<juce::AudioProcessorParameterGroup>(
                  spectral_settings_group_name,
                  "Spectral Settings",
                  " | ",
                  std::make_unique<juce::AudioParameterInt>(
                      fft_order_param_name,
                      "Frequency Resolution",
                      9,
                      15,
                      12,
                      "",
                      [](int value, int /*max_length*/) -> juce::String {
                          return juce::String(1 << value);
                      },
                      [](const juce::String& text) -> int {
                          return std::log2(text.getIntValue());
                      }),
                  // TODO: Change this to disallow non-power of 2 values
                  std::make_unique<juce::AudioParameterInt>(
                      windowing_overlap_times_param_name,
                      "Time Resolution",
                      2,
                      64,
                      4,
                      "x"
                      // TODO: We should show this in the GUI
                      // [&](int value, int /*max_length*/) -> juce::String {
                      //     return juce::String((1 << fft_order) / value);
                      // },
                      // [&](const juce::String& text) -> int {
                      //     return (1 << fft_order) / text.getIntValue();
                      // }
                      )),
          }),
      // TODO: Is this how you're supposed to retrieve non-float parameters?
      //       Seems a bit excessive
      sidechain_active(*dynamic_cast<juce::AudioParameterBool*>(
          parameters.getParameter(sidechain_active_param_name))),
      compressor_ratio(
          *parameters.getRawParameterValue(compressor_ratio_param_name)),
      auto_makeup_gain(*dynamic_cast<juce::AudioParameterBool*>(
          parameters.getParameter(auto_makeup_gain_param_name))),
      fft_order(*dynamic_cast<juce::AudioParameterInt*>(
          parameters.getParameter(fft_order_param_name))),
      windowing_overlap_times(*dynamic_cast<juce::AudioParameterInt*>(
          parameters.getParameter(windowing_overlap_times_param_name))),
      compressor_settings_listener(
          [&](const juce::String& /*parameterID*/, float /*newValue*/) {
              compressor_settings_changed = true;
          }),
      fft_order_listener(
          [&](const juce::String& /*parameterID*/, float /*newValue*/) {
              process_data_resizer.triggerAsyncUpdate();
              // FIXME: This should be done asynchronously, but this only works
              //        while the editor is open
              //        https://forum.juce.com/t/messages-are-only-being-handled-while-the-editor-is-open/45165/7
              if (!(getActiveEditor() && getActiveEditor()->isShowing())) {
                  process_data_resizer.handleUpdateNowIfNeeded();
              }
          }),
      process_data_resizer([&]() {
          process_data.resize_and_clear(static_cast<size_t>(fft_order));

          // The window size will also have changed, and we should thus update
          // our compressors to match the new effective sample rate.
          compressor_settings_changed = true;

          // TODO: The window size has to be changed after we change the FFT
          //       window size
          const size_t new_window_size = 1 << fft_order;
          setLatencySamples(new_window_size);
      }) {
    // TODO: Move the latency computation elsewhere
    const size_t new_window_size = 1 << fft_order;
    setLatencySamples(new_window_size);

    // XXX: There doesn't seem to be a fool proof way to just iterate over all
    //      parameters in a group, right?
    for (const auto& compressor_param_name :
         {sidechain_active_param_name, compressor_ratio_param_name,
          auto_makeup_gain_param_name, windowing_overlap_times_param_name}) {
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
    double /*sampleRate*/,
    int maximumExpectedSamplesPerBlock) {
    max_samples_per_block = static_cast<uint32>(maximumExpectedSamplesPerBlock);

    // We're explicitly initialzing the process data instead of using the resize
    // and clear here since resize_and_clear will temporarily set the resize
    // flag to false, and if `process_data_resizer` procs just before the first
    // block of audio gets processed then the swap never gets performed
    initialize_process_data(process_data.get(), static_cast<size_t>(fft_order));
}

void SpectralCompressorProcessor::releaseResources() {
    process_data.clear();
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
    ProcessData& data = process_data.get();

    // We need to maintain the same latency when bypassed, so we'll reuse most
    // of the processing logic
    do_stft(buffer, data, [this](ProcessData& data, size_t input_channels) {
        const size_t windowing_interval =
            data.fft_window_size / static_cast<size_t>(windowing_overlap_times);

        for (size_t channel = 0; channel < input_channels; channel++) {
            // We don't have a way to directly copy between buffers, but
            // most hosts should not actually hit this bypassed state
            // anyways
            // TODO: At some point, do implement this without using the scratch
            //       buffer
            data.input_ring_buffers[channel].copy_last_n_to(
                data.fft_scratch_buffer.data(), windowing_interval);
            data.output_ring_buffers[channel].read_n_from_in_place(
                data.fft_scratch_buffer.data(), windowing_interval);
        }
    });
}

void SpectralCompressorProcessor::processBlock(
    juce::AudioBuffer<float>& buffer,
    juce::MidiBuffer& /*midiMessages*/) {
    ProcessData& data = process_data.get();

    // We'll update the compressor settings just before processing if the
    // settings have changed or if the sidechaining has been disabled
    bool expected = true;
    if (compressor_settings_changed.compare_exchange_strong(expected, false)) {
        update_compressors(data);
    }

    // This function will let us process the input signal in windows, using
    // overlap-add
    do_stft(buffer, data, [this](ProcessData& data, size_t input_channels) {
        // If sidechaining is active, we set the compressor thresholds based on
        // a sidechain signal. Since compression is already ballistics based we
        // don't need any additional smoothing here.
        if (sidechain_active) {
            for (size_t channel = 0; channel < input_channels; channel++) {
                data.sidechain_ring_buffers[channel].copy_last_n_to(
                    data.fft_scratch_buffer.data(), data.fft_window_size);
                // TODO: We can skip negative frequencies here, right?
                data.fft->performRealOnlyForwardTransform(
                    data.fft_scratch_buffer.data(), true);

                // The version below is better annotated
                std::span<std::complex<float>> fft_buffer(
                    reinterpret_cast<std::complex<float>*>(
                        data.fft_scratch_buffer.data()),
                    data.fft_window_size);
                for (size_t compressor_idx = 0;
                     compressor_idx < data.spectral_compressors.size();
                     compressor_idx++) {
                    const size_t bin_idx = compressor_idx + 1;
                    const float magnitude = std::abs(fft_buffer[bin_idx]);

                    // We'll set the compressor threshold based on the
                    // arithmetic mean of the magnitudes of all channels. As a
                    // slight premature optimization (sorry) we'll reset these
                    // magnitudes after using them to avoid the conditional
                    // here.
                    data.spectral_compressor_sidechain_thresholds
                        [compressor_idx] += magnitude;
                }
            }

            for (size_t compressor_idx = 0;
                 compressor_idx < data.spectral_compressors.size();
                 compressor_idx++) {
                data.spectral_compressors[compressor_idx].setThreshold(
                    data.spectral_compressor_sidechain_thresholds
                        [compressor_idx] /
                    input_channels);
                data.spectral_compressor_sidechain_thresholds[compressor_idx] =
                    0;
            }
        }

        for (size_t channel = 0; channel < input_channels; channel++) {
            data.input_ring_buffers[channel].copy_last_n_to(
                data.fft_scratch_buffer.data(), data.fft_window_size);
            data.fft->performRealOnlyForwardTransform(
                data.fft_scratch_buffer.data());

            // We'll compress every FTT bin individually. Bin 0 is the DC offset
            // and should be skipped, and the latter half of the FFT bins should
            // be processed in the same way as the first half but in reverse
            // order. The real and imaginary parts are interleaved, so ever bin
            // spans two values in the scratch buffer. We can 'safely' do this
            // cast so we can use the STL's complex value functions.
            std::span<std::complex<float>> fft_buffer(
                reinterpret_cast<std::complex<float>*>(
                    data.fft_scratch_buffer.data()),
                data.fft_window_size);
            for (size_t compressor_idx = 0;
                 compressor_idx < data.spectral_compressors.size();
                 compressor_idx++) {
                // We don't have a compressor for the first bin
                const size_t bin_idx = compressor_idx + 1;

                // TODO: Are these _really_ exactly the same in the second
                //       half ergo this single magnitude is sufficient?
                const float magnitude = std::abs(fft_buffer[bin_idx]);
                const float compressed_magnitude =
                    data.spectral_compressors[compressor_idx].processSample(
                        channel, magnitude);

                // We need to scale both the imaginary and real components of
                // the bins at the start and end of the spectrum by the same
                // value
                // TODO: Add stereo linking
                const float compression_multiplier =
                    magnitude != 0.0f ? compressed_magnitude / magnitude : 1.0f;

                // The same operation should be applied to the mirrored bins at
                // the end of the FFT window, except for if this is the last bin
                fft_buffer[bin_idx] *= compression_multiplier;
                // TODO: Is this mirrored part necessary?
                if (compressor_idx != data.spectral_compressors.size() - 1) {
                    const size_t mirrored_bin_idx =
                        data.fft_window_size - bin_idx;
                    fft_buffer[mirrored_bin_idx] *= compression_multiplier;
                }
            }

            data.fft->performRealOnlyInverseTransform(
                data.fft_scratch_buffer.data());
            data.windowing_function->multiplyWithWindowingTable(
                data.fft_scratch_buffer.data(), data.fft_window_size);

            // After processing the windowed data, we'll add it to our output
            // ring buffer with any (automatic) makeup gain applied
            // TODO: We might need some kind of optional limiting stage to be
            //       safe
            data.output_ring_buffers[channel].add_n_from_in_place(
                data.fft_scratch_buffer.data(), data.fft_window_size,
                makeup_gain);
        }
    });
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

    compressor_settings_changed = true;

    // TODO: Do parameter listeners get triggered? Or alternatively, can this be
    //       called during playback (without `prepareToPlay()` being called
    //       first)?
    // TODO: Move the latency computation elsewhere
    const size_t new_window_size = 1 << fft_order;
    setLatencySamples(new_window_size);
}

void SpectralCompressorProcessor::initialize_process_data(
    ProcessData& inactive,
    size_t new_fft_order) {
    inactive.fft_window_size = 1 << new_fft_order;
    inactive.num_windows_processed = 0;

    if (fft_order > 0) {
        inactive.windowing_function.emplace(
            inactive.fft_window_size,
            juce::dsp::WindowingFunction<float>::WindowingMethod::hann,
            // TODO: Or should we leave normalization enabled?
            false);
        inactive.fft.emplace(new_fft_order);
    } else {
        inactive.windowing_function.reset();
        inactive.fft.reset();
    }

    // To prevent pops and artifacts caused by old samples, we'll clear all of
    // our buffers on a resize
    inactive.fft_scratch_buffer.clear();
    inactive.spectral_compressors.clear();
    inactive.spectral_compressor_sidechain_thresholds.clear();
    inactive.input_ring_buffers.clear();
    inactive.output_ring_buffers.clear();
    inactive.sidechain_ring_buffers.clear();
    if (new_fft_order <= 0) {
        inactive.fft_scratch_buffer.shrink_to_fit();
        inactive.spectral_compressors.shrink_to_fit();
        inactive.spectral_compressor_sidechain_thresholds.shrink_to_fit();
        inactive.input_ring_buffers.shrink_to_fit();
        inactive.output_ring_buffers.shrink_to_fit();
        inactive.sidechain_ring_buffers.shrink_to_fit();
        return;
    }

    // JUCE's FFT class interleaves the real and imaginary numbers, so this
    // buffer should be twice the window size in size
    inactive.fft_scratch_buffer.resize(inactive.fft_window_size * 2);

    // Every FFT bin on both channels gets its own compressor, hooray!  The
    // `(fft_window_size / 2) - 1` is because the first bin is the DC offset and
    // shouldn't be compressed, and the bins after the Nyquist frequency are the
    // same as the first half but in reverse order. The compressor settings will
    // be set in `update_compressors()`, which is triggered on the next
    // processing cycle by setting `compressor_settings_changed` below.
    juce::dsp::Compressor<float> compressor{};
    inactive.spectral_compressors.resize((inactive.fft_window_size / 2) - 1,
                                         compressor);
    inactive.spectral_compressor_sidechain_thresholds.resize(
        inactive.spectral_compressors.size());

    // We use ring buffers to store the samples we'll process using FFT and also
    // to store the samples that should be played back to
    inactive.input_ring_buffers.resize(
        static_cast<size_t>(getMainBusNumInputChannels()),
        RingBuffer<float>(inactive.fft_window_size));
    inactive.output_ring_buffers.resize(
        static_cast<size_t>(getMainBusNumOutputChannels()),
        RingBuffer<float>(inactive.fft_window_size));
    inactive.sidechain_ring_buffers.resize(
        static_cast<size_t>(getChannelCountOfBus(true, 1)),
        RingBuffer<float>(inactive.fft_window_size));

    // After resizing the compressors are uninitialized and should be
    // reinitialized
    compressor_settings_changed = true;
}

void SpectralCompressorProcessor::update_compressors(ProcessData& data) {
    const double effective_sample_rate =
        getSampleRate() /
        (static_cast<double>(data.fft_window_size) / windowing_overlap_times);
    for (size_t compressor_idx = 0;
         compressor_idx < data.spectral_compressors.size(); compressor_idx++) {
        auto& compressor = data.spectral_compressors[compressor_idx];

        // TODO: Make the timings configurable
        compressor.setRatio(compressor_ratio);
        compressor.setAttack(50.0);
        compressor.setRelease(5000.0);
        // TODO: This prepare resets the envelope follower, which is not what we
        //       want. In our own compressor we should have a way to just change
        //       the sample rate.
        // TODO: Now that the timings are compensated for changing window
        //       intervals, we might not need this to be configurable anymore
        //       can just leave this fixed at 4x.
        compressor.prepare(juce::dsp::ProcessSpec{
            // We only process everything once every `windowing_interval`,
            // otherwise our attack and release times will be all messed up
            .sampleRate = effective_sample_rate,
            .maximumBlockSize = max_samples_per_block,
            .numChannels = static_cast<uint32>(getMainBusNumInputChannels())});
    }

    // TODO: The user should be able to configure their own slope (or free
    //       drawn)
    // TODO: And we should be doing both upwards and downwards compression,
    //       OTT-style
    constexpr float base_threshold_dbfs = 0.0f;
    if (!sidechain_active) {
        // The thresholds are set to match pink noise.
        // TODO: Change the calculations so that the base threshold parameter is
        //       centered around some frequency
        const float frequency_increment =
            getSampleRate() / data.fft_window_size;
        for (size_t compressor_idx = 0;
             compressor_idx < data.spectral_compressors.size();
             compressor_idx++) {
            // The first bin doesn't get a compressor
            const size_t bin_idx = compressor_idx + 1;
            const float frequency = frequency_increment * bin_idx;

            // This starts at 1 for 0 Hz (DC)
            const float octave = std::log2(frequency + 2);

            // The 3 dB is to compensate for bin 0
            const float threshold =
                (base_threshold_dbfs + 3.0f) - (3.0f * octave);
            data.spectral_compressors[compressor_idx].setThreshold(threshold);
        }
    }

    // We need to compensate for the extra gain added by 4x overlap
    makeup_gain = 1.0f / windowing_overlap_times;
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
}

juce::AudioProcessor* JUCE_CALLTYPE createPluginFilter() {
    return new SpectralCompressorProcessor();
}
