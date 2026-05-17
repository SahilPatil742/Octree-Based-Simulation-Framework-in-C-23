#pragma once

#include <array>
#include <cstddef>
#include <initializer_list>
#include <type_traits>
#include <cmath>
#include <ostream>
#include <concepts>

namespace oktal {

template <typename T, std::size_t DIM>
class Vec {
public:
    using value_type = T;
    using size_type  = std::size_t;

private:
    std::array<T, DIM> v_{}; // zero-initialized

public:
    // --- Constructors ---
    constexpr Vec() noexcept = default;

    constexpr explicit Vec(const T& value) noexcept {
        for (size_type i = 0; i < DIM; ++i) {
            v_.at(i) = value;
        }
    }

    constexpr Vec(std::initializer_list<T> il) noexcept {
        size_type i = 0;
        for (auto it = il.begin(); it != il.end() && i < DIM; ++it, ++i) {
            v_.at(i) = *it;
        }
        for (; i < DIM; ++i) {
            v_.at(i) = T{};
        }
    }

    template <typename S>
        requires std::convertible_to<S, T>
    constexpr explicit Vec(const Vec<S, DIM>& other) noexcept {
        for (size_type i = 0; i < DIM; ++i) {
            v_.at(i) = static_cast<T>(other[i]);
        }
    }

    // --- Observers ---
    static constexpr size_type dimension() noexcept { return DIM; }

    [[nodiscard]] constexpr size_type size() const noexcept { return DIM; }

    [[nodiscard]] constexpr T* data() noexcept { return v_.data(); }
    [[nodiscard]] constexpr const T* data() const noexcept { return v_.data(); }

    // --- Iterators ---
    [[nodiscard]] constexpr auto begin() noexcept { return v_.begin(); }
    [[nodiscard]] constexpr auto end() noexcept { return v_.end(); }
    [[nodiscard]] constexpr auto begin() const noexcept { return v_.begin(); }
    [[nodiscard]] constexpr auto end() const noexcept { return v_.end(); }
    [[nodiscard]] constexpr auto cbegin() const noexcept { return v_.cbegin(); }
    [[nodiscard]] constexpr auto cend() const noexcept { return v_.cend(); }

    // --- Element access ---
    [[nodiscard]] constexpr T& operator[](size_type i) noexcept {
        return v_.at(i);
    }

    [[nodiscard]] constexpr const T& operator[](size_type i) const noexcept {
        return v_.at(i);
    }

    // --- Equality ---
    constexpr bool operator==(const Vec& other) const noexcept {
        for (size_type i = 0; i < DIM; ++i) {
            if (!(v_.at(i) == other.v_.at(i))) {
                return false;
            }
        }
        return true;
    }

    constexpr bool operator!=(const Vec& other) const noexcept {
        return !(*this == other);
    }

    // --- Unary negation ---
    constexpr Vec operator-() const noexcept {
        Vec out;
        for (size_type i = 0; i < DIM; ++i) {
            out.v_.at(i) = -v_.at(i);
        }
        return out;
    }

    // --- Augmented assignments (vector) ---
    constexpr Vec& operator+=(const Vec& rhs) noexcept {
        for (size_type i = 0; i < DIM; ++i) {
            v_.at(i) += rhs.v_.at(i);
        }
        return *this;
    }

    constexpr Vec& operator-=(const Vec& rhs) noexcept {
        for (size_type i = 0; i < DIM; ++i) {
            v_.at(i) -= rhs.v_.at(i);
        }
        return *this;
    }

    // --- Augmented assignments (scalar) ---
    template <typename S>
        requires std::convertible_to<S, T>
    constexpr Vec& operator*=(S scalar) noexcept {
        const T s = static_cast<T>(scalar);
        for (size_type i = 0; i < DIM; ++i) {
            v_.at(i) *= s;
        }
        return *this;
    }

    template <typename S>
        requires std::convertible_to<S, T>
    constexpr Vec& operator/=(S scalar) noexcept {
        const T s = static_cast<T>(scalar);
        for (size_type i = 0; i < DIM; ++i) {
            v_.at(i) /= s;
        }
        return *this;
    }

    // --- Norms ---
    [[nodiscard]] constexpr T sqrMagnitude() const noexcept {
        long double acc = 0.0L;
        for (size_type i = 0; i < DIM; ++i) {
            const auto x = static_cast<long double>(v_.at(i));
            acc += x * x;
        }
        return static_cast<T>(acc);
    }

    [[nodiscard]] T magnitude() const noexcept {
        const auto sq = static_cast<long double>(sqrMagnitude());
        return static_cast<T>(std::sqrt(sq));
    }
};

// ----------------------------- Non-member operators -----------------------------

template <typename T, std::size_t DIM>
constexpr Vec<T, DIM> operator+(Vec<T, DIM> lhs, const Vec<T, DIM>& rhs) noexcept {
    lhs += rhs;
    return lhs;
}

template <typename T, std::size_t DIM>
constexpr Vec<T, DIM> operator-(Vec<T, DIM> lhs, const Vec<T, DIM>& rhs) noexcept {
    lhs -= rhs;
    return lhs;
}

template <typename T, std::size_t DIM, typename S>
    requires std::convertible_to<S, T>
constexpr Vec<T, DIM> operator*(Vec<T, DIM> v, S scalar) noexcept {
    v *= scalar;
    return v;
}

template <typename T, std::size_t DIM, typename S>
    requires std::convertible_to<S, T>
constexpr Vec<T, DIM> operator*(S scalar, Vec<T, DIM> v) noexcept {
    v *= scalar;
    return v;
}

template <typename T, std::size_t DIM, typename S>
    requires std::convertible_to<S, T>
constexpr Vec<T, DIM> operator/(Vec<T, DIM> v, S scalar) noexcept {
    v /= scalar;
    return v;
}

// ----------------------------- Streaming -----------------------------

template <typename T, std::size_t DIM>
inline std::ostream& operator<<(std::ostream& os, const Vec<T, DIM>& v) {
    os << "{";
    for (std::size_t i = 0; i < DIM; ++i) {
        if (i) {
            os << ", ";
        }
        os << v[i];
    }
    os << "}";
    return os;
}

// ----------------------------- Type aliases -----------------------------

using Vec3F = Vec<float, 3>;
using Vec3D = Vec<double, 3>;

} // namespace oktal
