#pragma once

// We use <oktal/geometry/Vec.hpp> to ensure the compiler finds the file
// from the project's include root, preventing file-not-found errors.
#include <oktal/geometry/Vec.hpp>

namespace oktal {

    template <typename T = double>
    class Box {
    public:
        // 1. Required public type alias
        // This must match Vec<T, 3> exactly for the static_asserts in TestBox.cpp
        using vector_type = Vec<T, 3>;

    private:
        // 2. Backing private data members
        vector_type m_minCorner;
        vector_type m_maxCorner;

    public:
        // 3. Default constructor
        // Your Vec() defaults to zero-initialization, so we rely on that.
        constexpr Box() noexcept = default;

        // 4. Constructor taking two vectors
        constexpr Box(const vector_type& min, const vector_type& max) noexcept
            : m_minCorner(min), m_maxCorner(max) {
        }

        // 5. Accessors (Mutable) - REQUIRED for testSetters
        // These MUST return a reference (vector_type&) so 'b.minCorner() = ...' modifies the box.
        [[nodiscard]] constexpr vector_type& minCorner() noexcept {
            return m_minCorner;
        }

        [[nodiscard]] constexpr vector_type& maxCorner() noexcept {
            return m_maxCorner;
        }

        // 5b. Accessors (Const) - REQUIRED for testConstructorsAndGetters
        // These MUST be present for 'const Box b' usage.
        [[nodiscard]] constexpr const vector_type& minCorner() const noexcept {
            return m_minCorner;
        }

        [[nodiscard]] constexpr const vector_type& maxCorner() const noexcept {
            return m_maxCorner;
        }

        // --- Observers ---

        // 6. extents(): returns (max - min)
        // Uses your Vec::operator-
        [[nodiscard]] constexpr vector_type extents() const noexcept {
            return m_maxCorner - m_minCorner;
        }

        // 7. center(): returns center point
        // Uses Vec::operator+ and Vec::operator/
        // We use static_cast<T>(2) to ensure the scalar type matches T exactly
        // (e.g., 2.0f for float, 2.0 for double)
        [[nodiscard]] constexpr vector_type center() const noexcept {
            return (m_minCorner + m_maxCorner) / static_cast<T>(2);
        }

        // 8. volume(): returns the volume of the box
        [[nodiscard]] constexpr T volume() const noexcept {
            const vector_type e = extents();
            // Uses your Vec::operator[]
            return e[0] * e[1] * e[2];
        }
    };

} // namespace oktal