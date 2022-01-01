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

#include <juce_audio_processors/juce_audio_processors.h>
#include <function2/function2.hpp>

/**
 * Run some function on the message thread. This function will be executed
 * synchronously and should thus run in constant time.
 */
class LambdaAsyncUpdater : public juce::AsyncUpdater {
   public:
    LambdaAsyncUpdater(fu2::unique_function<void()> callback);

    void handleAsyncUpdate() override;

   private:
    fu2::unique_function<void()> callback_;

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
    fu2::unique_function<void(const juce::String&, float)> callback_;

    JUCE_DECLARE_NON_COPYABLE_WITH_LEAK_DETECTOR(LambdaParameterListener)
};

/**
 * A wrapper around some `T` that contains an active `T` and an inactive `T`,
 * with a pointer pointing to the currently active object. When some plugin
 * parameter changes that would require us to resize the object, we can resize
 * the inactive object and then swap the two pointers on the next time we fetch
 * the object from the audio processing loop. This prevents locking and memory
 * allocations on the audio thread. Keep in mind that the active and the
 * inactive objects have no relation to each other, and might thus contain
 * completely different data.
 */
template <typename T>
class AtomicallySwappable {
   public:
    /**
     * Default initalizes the objects.
     */
    AtomicallySwappable()
        : pointers_(Pointers{.active = &primary_, .inactive = &secondary_}),
          primary_(),
          secondary_() {}

    /**
     * Initialize the objects with some default value.
     *
     * @param initial The initial value for the object. This will also be copied
     *   to the inactive slot.
     */
    AtomicallySwappable(T initial)
        : pointers_(Pointers{.active = &primary_, .inactive = &secondary_}),
          primary_(initial),
          secondary_(initial) {}

    /**
     * Return a reference to currently active object. This should be done once
     * at the start of the audio processing function, and the same reference
     * should be reused for the remainder of the function.
     */
    T& get() {
        // We'll swap the pointer on the audio thread so that two resizes in a
        // row in between audio processing calls don't cause weird behaviour
        bool expected = true;
        if (needs_swap_.compare_exchange_strong(expected, false)) {
            // The CaS should be atomic, even though GCC will always return
            // false for the `is_lock_free()`/`is_always_lock_free()` on 128-bit
            // types
            static_assert(sizeof(Pointers) == sizeof(T* [2]));

            Pointers current_pointers, updated_pointers;
            do {
                current_pointers = pointers_;
                updated_pointers =
                    Pointers{.active = current_pointers.inactive,
                             .inactive = current_pointers.active};
            } while (!pointers_.compare_exchange_weak(current_pointers,
                                                      updated_pointers));
        }

        return *pointers_.load().active;
    }

    /**
     * Modify the inactive object using the supplied function, and swap the
     * active and the inactive objects on the next call to `get()`. This may
     * block and should thus never be called from the audio thread.
     *
     * @tparam F A function with the signature `void(T&)`.
     */
    template <typename F>
    void modify_and_swap(F modify_fn) {
        // In case two mutations are performed in a row, we don't want the audio
        // thread swapping the objects while we're modifying that same object
        // from another thread
        num_resizing_threads_.fetch_add(1);
        needs_swap_ = false;

        std::lock_guard lock(resize_mutex_);
        modify_fn(*pointers_.load().inactive);

        // If for whatever reason multiple threads are calling this function at
        // the same time, then only the last one may set the swap flag to
        // prevent (admittedly super rare) data races
        if (num_resizing_threads_.fetch_sub(1) == 1) {
            needs_swap_ = true;
        }
    }

    /**
     * Resize both objects down to their smallest size using the supplied
     * function. This should only ever be called from
     * `AudioProcessor::releaseResources()`.
     *
     * @tparam F A function with the signature `void(T&)`.
     */
    template <typename F>
    void clear(F clear_fn) {
        std::lock_guard lock(resize_mutex_);

        clear_fn(primary_);
        clear_fn(secondary_);
    }

   private:
    /**
     * In the unlikely situation that two threads are calling resize at the same
     * time, we'll use a mutex to make sure that those two resizes aren't
     * happening at the same time and we use this `num_resizing_threads` to make
     * sure that `needs_swap` only gets set to `true` when both threads are
     * done. This is to prevent a (super rare) race condition where the audio
     * thread will CaS `needs_swap` to false and swap the active pointer while
     * at the same time another who just got access to the resize mutex is
     * working on the now active object.
     */
    std::atomic_int num_resizing_threads_ = 0;
    std::mutex resize_mutex_;

    struct Pointers {
        T* active;
        T* inactive;
    };
    std::atomic_bool needs_swap_ = false;
    std::atomic<Pointers> pointers_;

    T primary_;
    T secondary_;
};
