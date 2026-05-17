#pragma once

#include <array>
#include <cmath>
#include <cstdint>
#include <cstddef>
#include <vector>
#include <stdexcept>
#include <bit>

#include <oktal/geometry/Vec.hpp>

namespace oktal {

using morton_bits_t = std::uint64_t;

class MortonIndex {
public:
    // Maximum representable octree depth
    static constexpr std::size_t MAX_DEPTH =
        (sizeof(morton_bits_t) * 8) / 3;

    // Root cell constructor
    MortonIndex() noexcept : bits_{morton_bits_t{1}} {}

    // Construct from raw Morton bit pattern
    MortonIndex(morton_bits_t bits);

    // Access underlying bit pattern
    [[nodiscard]] morton_bits_t getBits() const;

    [[nodiscard]] static MortonIndex
    fromPath(const std::vector<morton_bits_t>& path);

    [[nodiscard]] std::vector<morton_bits_t> getPath() const;

    [[nodiscard]] std::size_t level() const;
    [[nodiscard]] bool isRoot() const;

    [[nodiscard]] std::size_t siblingIndex() const;
    [[nodiscard]] bool isFirstSibling() const;
    [[nodiscard]] bool isLastSibling() const;

    [[nodiscard]] MortonIndex parent() const;
    [[nodiscard]] MortonIndex safeParent() const;

    [[nodiscard]] MortonIndex child(morton_bits_t c) const;
    [[nodiscard]] MortonIndex safeChild(morton_bits_t c) const;

    [[nodiscard]] bool operator==(const MortonIndex& other) const;
    [[nodiscard]] bool operator!=(const MortonIndex& other) const;

    [[nodiscard]] bool operator<(const MortonIndex& other) const;
    [[nodiscard]] bool operator>(const MortonIndex& other) const;

    [[nodiscard]] bool operator<=(const MortonIndex& other) const;
    [[nodiscard]] bool operator>=(const MortonIndex& other) const;

    [[nodiscard]] Vec<std::size_t, 3> gridCoordinates() const;

    [[nodiscard]] static MortonIndex fromGridCoordinates(
        std::size_t level,
        const Vec<std::size_t, 3>& coords);

private:
    morton_bits_t bits_{1};   // root by default
};

} // namespace oktal
