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

// TODO: Move this back to 4096 and introduce latency
constexpr int fft_window_size = 512;
constexpr int fft_order = 9;

static_assert(1 << fft_order == fft_window_size,
              "The FFT order and FFT window sizes don't match up");

// TODO: Rewrite, this is from the example

//==============================================================================
class SpectralCompressorProcessor : public juce::AudioProcessor {
   public:
    //==============================================================================
    SpectralCompressorProcessor();
    ~SpectralCompressorProcessor() override;

    //==============================================================================
    void prepareToPlay(double sampleRate,
                       int maximumExpectedSamplesPerBlock) override;
    void releaseResources() override;

    bool isBusesLayoutSupported(const BusesLayout& layouts) const override;

    void processBlock(juce::AudioBuffer<float>&, juce::MidiBuffer&) override;
    using AudioProcessor::processBlock;

    //==============================================================================
    juce::AudioProcessorEditor* createEditor() override;
    bool hasEditor() const override;

    //==============================================================================
    const juce::String getName() const override;

    bool acceptsMidi() const override;
    bool producesMidi() const override;
    bool isMidiEffect() const override;
    double getTailLengthSeconds() const override;

    //==============================================================================
    int getNumPrograms() override;
    int getCurrentProgram() override;
    void setCurrentProgram(int index) override;
    const juce::String getProgramName(int index) override;
    void changeProgramName(int index, const juce::String& newName) override;

    //==============================================================================
    void getStateInformation(juce::MemoryBlock& destData) override;
    void setStateInformation(const void* data, int sizeInBytes) override;

   private:
    juce::dsp::FFT fft;
    /**
     * We need a scratch buffer that can contain `fft_window_size * 2` samples.
     */
    std::vector<float> fft_scratch_buffer;

    /**
     * This will contain `fft_window_size` compressors. We'll compress the
     * magnitude of every FFT bin (`sqrt(i^2 + r^2)`) individually, and then
     * scale both the real and imaginary components by the ratio of their
     * magnitude and the compressed value.
     */
    std::vector<juce::dsp::Compressor<float>> spectral_compressors;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(SpectralCompressorProcessor)
};
