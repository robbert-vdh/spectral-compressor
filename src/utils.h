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

#include <atomic>
#include <mutex>

#include <function2/function2.hpp>

/**
 * A wrapper around some resizeable type `T` that contains an active `T` and an
 * inactive `T`, with a pointer pointing to the currently active object. When
 * resizing, we'll resize the inactive `T`, and then set a flag that will cause
 * the active and the inactive objects to get swapped the next time the audio
 * thread requests a reference to the currently active objects. This prevents
 * locking and memory allocations on the audio thread.
 */
template <typename T>
class AtomicResizable {
   public:
    /**
     * Default initalizes the objects.
     *
     * @param resize_and_clear_fn This function should resize an object of type
     *   `T` and potentially also clear its values. A new size of 0 means that
     *   the object has to release its resources. While not strictly necessary,
     *   clearing may be a good idea to avoid weird pops and other artifacts.
     */
    AtomicResizable(fu2::unique_function<void(T&, size_t)> resize_and_clear_fn)
        : resize_and_clear_fn(std::move(resize_and_clear_fn)),
          active(&primary),
          primary(),
          secondary() {}

    /**
     * Initialize the objects with some default value.
     *
     * @param initial The initial value for the object. This will also be copied
     *   to the inactive slot.
     * @param resize_and_clear_fn This function should resize an object of type
     *   `T` and potentially also clear its values. A new size of 0 means that
     *   the object has to release its resources. While not strictly necessary,
     *   clearing may be a good idea to avoid weird pops and other artifacts.
     */
    AtomicResizable(T initial,
                    fu2::unique_function<void(T&, size_t)> resize_and_clear_fn)
        : resize_and_clear_fn(std::move(resize_and_clear_fn)),
          active(&primary),
          primary(initial),
          secondary(initial) {}

    /**
     * Return a reference to currently active object. This should be done at the
     * start of the audio processing function, and the same reference should be
     * reused for the remainder of the function.
     */
    T& get() {
        // We'll swap the pointer on the audio thread so that two resizes in a
        // row in between audio processing calls don't cause weird behaviour
        bool expected = true;
        if (needs_swap.compare_exchange_strong(expected, false)) {
            active = active == &primary ? &secondary : &primary;
        }

        return *active;
    }

    /**
     * Resize and clear the object. This may block and should never be called
     * from the audio thread.
     */
    void resize_and_clear(size_t new_size) {
        std::lock_guard lock(resize_mutex);

        // In case two resizes are performed in a row, we don't want the audio
        // thread swapping the objects while we're performing a second resize
        // TODO: This also isn't entirely safe, although it does require all the
        //       stars to be aligned for this to cause an issue. If a background
        //       thread calls this function once, and then calls it again, then
        //       it can happen that the audio thread does a CaS on `needs_swap`
        //       and starts swapping the pointers in the small time span that
        //       `needs_swap` is set to true. This could cause us to call
        //       `resize_and_clear_fn()` on the same object that gets passed to
        //       the audio thread.
        needs_swap = false;
        resize_and_clear_fn(active == &primary ? secondary : primary, new_size);
        needs_swap = true;
    }

    /**
     * Resize both objects down to their smallest size. This should only ever be
     * called from `AudioProcessor::releaseResources()`.
     */
    void clear() {
        std::lock_guard lock(resize_mutex);

        resize_and_clear_fn(primary, 0);
        resize_and_clear_fn(secondary, 0);
    }

   private:
    fu2::unique_function<void(T&, size_t)> resize_and_clear_fn;

    std::atomic_bool needs_swap = false;
    std::mutex resize_mutex;

    std::atomic<T*> active;
    T primary;
    T secondary;
};
