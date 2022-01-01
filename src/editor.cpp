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

#include "editor.h"

#include "processor.h"

SpectralCompressorEditor::SpectralCompressorEditor(
    SpectralCompressorProcessor& p)
    : AudioProcessorEditor(&p), processor_(p) {
    setSize(400, 300);
}

SpectralCompressorEditor::~SpectralCompressorEditor() {}

//==============================================================================
void SpectralCompressorEditor::paint(juce::Graphics& g) {
    // TODO: Replace with something else. Or drop the GUI.
    g.fillAll(
        getLookAndFeel().findColour(juce::ResizableWindow::backgroundColourId));

    g.setColour(juce::Colours::white);
    g.setFont(15.0f);
    g.drawFittedText("Hello World!", getLocalBounds(),
                     juce::Justification::centred, 1);
}

void SpectralCompressorEditor::resized() {
    // TODO
}
