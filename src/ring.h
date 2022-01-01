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

#include <juce_audio_basics/juce_audio_basics.h>
#include <juce_core/juce_core.h>

/**
 * A simple resizeable ring buffer that allows copying up to its size number of
 * samples to and from the buffer at a time.
 *
 * @tparam T The element type of this ring buffer. Because of the operations
 *   used, this can only be `float` or `double`.
 */
template <typename T>
class RingBuffer {
   public:
    /**
     * The default constructor doesn't initialize the ring buffer.
     * `RingBuffer::resize()` should be called before actually using this.
     *
     * @see RingBuffer::resize()
     */
    RingBuffer() {}

    /**
     * Initialize the ring buffer to contain `size` `T`s.
     */
    RingBuffer(size_t size) : buffer_(size, 0.0) {}

    /**
     * Resize the ring buffer to be able to contain `new_size` elements. This
     * will reset the current position to 0. Existing data will not be cleared.
     */
    void resize(size_t new_size) {
        buffer_.resize(new_size);
        current_pos_ = 0;
    }

    /**
     * Returns the ring buffer's current size.
     */
    inline size_t size() const { return buffer_.size(); }

    /**
     * Returns the current head position in the ring buffer.
     */
    inline size_t pos() const { return current_pos_; }

    /**
     * Copy `num` samples from `src` into the ring buffer, starting at `pos()`.
     *
     * This advances the current position by `num`.
     *
     * @param src The buffer to copy from.
     * @param num How many elements to read, should not exceed `size()`.
     *
     * @return The number of elements read.
     *
     * @throw std::invalid_argument When `num > size()`.
     */
    size_t read_n_from(const T* src, size_t num) {
        if (num > buffer_.size()) {
            throw std::invalid_argument(
                "num > size() in RingBuffer::read_n_from()");
        }

        const auto& [num_to_end, num_from_start] =
            split_range_from(current_pos_, num);
        std::copy_n(src, num_to_end, &buffer_[current_pos_]);
        std::copy_n(src + num_to_end, num_from_start, &buffer_[0]);

        current_pos_ += num;
        if (current_pos_ >= buffer_.size()) {
            current_pos_ -= buffer_.size();
        }

        return num;
    }

    /**
     * Copy `num` samples (starting at `pos()`) to `dst`.
     *
     * This advances the current position by `num`.
     *
     * @param dst Where to write the data to.
     * @param num How many elements to write, should not exceed `size()`.
     * @param clear Whether we should overwrite the in our buffer we just copied
     *   to `dst` with `0.0`. We'll use this when writing processed output back
     *   since we're adding up the results of overlapping processed regions and
     *   we definitely don't want any feedback.
     *
     * @return The number of elements copied.
     *
     * @throw std::invalid_argument When `num > size()`.
     */
    size_t copy_n_to(T* dst, size_t num, bool clear) {
        if (num > buffer_.size()) {
            throw std::invalid_argument(
                "num > size() in RingBuffer::copy_n_to()");
        }

        const auto& [num_to_end, num_from_start] =
            split_range_from(current_pos_, num);
        std::copy_n(&buffer_[current_pos_], num_to_end, dst);
        std::copy_n(&buffer_[0], num_from_start, dst + num_to_end);
        if (clear) {
            std::fill_n(&buffer_[current_pos_], num_to_end, 0.0);
            std::fill_n(&buffer_[0], num_from_start, 0.0);
        }

        current_pos_ += num;
        if (current_pos_ >= buffer_.size()) {
            current_pos_ -= buffer_.size();
        }

        return num;
    }

    // The following operations are similar to the reading and writing functions
    // above, but they do not move the current position

    /**
     * Add `num` samples from `src` to the existing values in the ring buffer,
     * starting at `pos()`. We use this to process overlapping windows. The
     * `clear` flag in `copy_n_to()` is used to prevent these additions from
     * causing feedback.
     *
     * This does not advance the current position.
     *
     * @param src The buffer to copy from.
     * @param num How many elements to read, should not exceed `size()`.
     * @param gain A gain multiplier before adding the values. If set to `1.0`
     *   then no gain will be added.
     *
     * @return The number of elements copied.
     *
     * @throw std::invalid_argument When `num > size()`.
     */
    size_t add_n_from_in_place(const T* src, size_t num, float gain = 1.0) {
        if (num > buffer_.size()) {
            throw std::invalid_argument(
                "num > size() in RingBuffer::copy_n_to()");
        }

        const auto& [num_to_end, num_from_start] =
            split_range_from(current_pos_, num);
        if (gain == 1.0) {
            juce::FloatVectorOperations::add(&buffer_[current_pos_], src,
                                             num_to_end);
            juce::FloatVectorOperations::add(&buffer_[0], src + num_to_end,
                                             num_from_start);
        } else {
            juce::FloatVectorOperations::addWithMultiply(&buffer_[current_pos_],
                                                         src, gain, num_to_end);
            juce::FloatVectorOperations::addWithMultiply(
                &buffer_[0], src + num_to_end, gain, num_from_start);
        }

        return num;
    }

    /**
     * Copy `num` samples from `src` to the ring buffer, starting at `pos()`.
     * This is only used when the plugin is bypassed to maintain the proper
     * latency.
     *
     * This does not advance the current position.
     *
     * @param src The buffer to copy from.
     * @param num How many elements to read, should not exceed `size()`.
     *
     * @return The number of elements copied.
     *
     * @throw std::invalid_argument When `num > size()`.
     */
    size_t read_n_from_in_place(const T* src, size_t num) {
        if (num > buffer_.size()) {
            throw std::invalid_argument(
                "num > size() in RingBuffer::copy_n_to()");
        }

        const auto& [num_to_end, num_from_start] =
            split_range_from(current_pos_, num);
        std::copy_n(src, num_to_end, &buffer_[current_pos_]);
        std::copy_n(src + num_to_end, num_from_start, &buffer_[0]);

        return num;
    }

    /**
     * Copy the _last_ `num` samples (going backwards at `pos()`) written to
     * this ring buffer to `dst`. In our case we'll likely read the entire ring
     * buffer at once (i.e. `num == size()`).
     *
     * This does not advance the current position.
     *
     * @param dst Where to write the data to.
     * @param num How many elements to write, should not exceed `size()`.
     *
     * @return The number of elements copied.
     *
     * @throw std::invalid_argument When `num > size()`.
     */
    size_t copy_last_n_to(T* dst, size_t num) {
        if (num > buffer_.size()) {
            throw std::invalid_argument(
                "num > size() in RingBuffer::copy_n_to()");
        }

        // Like all other C-family languages, the % operator is a remainder
        // operator instead of a modulus so you need this abomination when
        // dealing with negative numbers
        const size_t start_pos =
            (current_pos_ - num + buffer_.size()) % buffer_.size();

        const auto& [num_to_end, num_from_start] =
            split_range_from(start_pos, num);
        std::copy_n(&buffer_[start_pos], num_to_end, dst);
        std::copy_n(&buffer_[0], num_from_start, dst + num_to_end);

        return num;
    }

   private:
    /**
     * Returns how to split the range when reading or writing `num` elements
     * starting at `from` in this buffer:
     *
     * ```cpp
     * const auto& [num_to_end, num_from_start] =
     *     split_range_from(pos, num);
     * // Do something with buffer[pos, pos + num_to_end]
     * // Do something with buffer[0, num_from_start] if num_from_start > 0
     * ```
     */
    std::pair<size_t, size_t> split_range_from(size_t from, size_t num) {
        const size_t num_to_end = std::min(num, buffer_.size() - from);
        return std::pair(num_to_end, num - num_to_end);
    }

    std::vector<T> buffer_;
    size_t current_pos_ = 0;
};
