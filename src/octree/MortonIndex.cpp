#include "oktal/octree/MortonIndex.hpp"

#include <algorithm>   // std::ranges::reverse
#include <bit>
#include <ranges>
#include <stdexcept>
#include <vector>

namespace oktal {

// ------------------------------------------------------------
// Construct from raw Morton bit pattern
// ------------------------------------------------------------
MortonIndex::MortonIndex(morton_bits_t bits)
    : bits_{bits} {}

// ------------------------------------------------------------
// Return underlying bit pattern
// ------------------------------------------------------------
morton_bits_t MortonIndex::getBits() const {
    return bits_;
}

// ------------------------------------------------------------
// Create MortonIndex from a path
// ------------------------------------------------------------
MortonIndex MortonIndex::fromPath(
    const std::vector<morton_bits_t>& path)
{
    auto bits = morton_bits_t{1}; // root = 1
    for (const morton_bits_t c : path) {
        bits = (bits << 3) | (c & 0b111);
    }
    return {bits};
}

// ------------------------------------------------------------
// Convert MortonIndex back to a path
// ------------------------------------------------------------
std::vector<morton_bits_t> MortonIndex::getPath() const {
    std::vector<morton_bits_t> path;

    if (bits_ == 1) {
        return path; // root
    }

    auto temp = bits_;
    while (temp > 1) {
        path.push_back(temp & 0b111);
        temp >>= 3;
    }

    std::ranges::reverse(path);
    return path;
}

// ------------------------------------------------------------
// Level (depth) in the tree
// ------------------------------------------------------------
std::size_t MortonIndex::level() const {
    auto temp = bits_;
    std::size_t lvl = 0;

    while (temp > 1) {
        temp >>= 3;
        ++lvl;
    }

    return lvl;
}

bool MortonIndex::isRoot() const {
    return bits_ == 1;
}

// ------------------------------------------------------------
// Sibling queries
// ------------------------------------------------------------
std::size_t MortonIndex::siblingIndex() const {
    return isRoot() ? 0U : static_cast<std::size_t>(bits_ & 0b111);
}

bool MortonIndex::isFirstSibling() const {
    return isRoot() || siblingIndex() == 0;
}

bool MortonIndex::isLastSibling() const {
    return isRoot() || siblingIndex() == 7;
}

// ------------------------------------------------------------
// Parent / child navigation
// ------------------------------------------------------------
MortonIndex MortonIndex::parent() const {
    return MortonIndex{bits_ >> 3};
}

MortonIndex MortonIndex::safeParent() const {
    if (isRoot()) {
        throw std::logic_error("Root cell has no parent");
    }
    return parent();
}

MortonIndex MortonIndex::child(morton_bits_t c) const {
    return MortonIndex{(bits_ << 3) | (c & 0b111)};
}

MortonIndex MortonIndex::safeChild(morton_bits_t c) const {
    if (c >= 8) {
        throw std::logic_error("Invalid octal child index");
    }
    if (level() >= MAX_DEPTH) {
        throw std::logic_error("Maximum Morton depth exceeded");
    }
    return child(c);
}

// ------------------------------------------------------------
// Comparisons
// ------------------------------------------------------------
bool MortonIndex::operator==(const MortonIndex& other) const {
    return bits_ == other.bits_;
}

bool MortonIndex::operator!=(const MortonIndex& other) const {
    return !(*this == other);
}

bool MortonIndex::operator>(const MortonIndex& other) const {
    const std::size_t this_level  = level();
    const std::size_t other_level = other.level();

    if (this_level >= other_level) {
        return false;
    }

    const std::size_t shift = 3 * (other_level - this_level);
    return bits_ == (other.bits_ >> shift);
}

bool MortonIndex::operator<(const MortonIndex& other) const {
    return other > *this;
}

bool MortonIndex::operator<=(const MortonIndex& other) const {
    return (*this < other) || (*this == other);
}

bool MortonIndex::operator>=(const MortonIndex& other) const {
    return (*this > other) || (*this == other);
}

// ------------------------------------------------------------
// Grid coordinates
// ------------------------------------------------------------
Vec<std::size_t, 3> MortonIndex::gridCoordinates() const {
    Vec<std::size_t, 3> coords{0, 0, 0};
    const auto path = getPath();

    for (const auto c : path) {
        coords[0] = (coords[0] << 1) | ((c >> 0) & 1);
        coords[1] = (coords[1] << 1) | ((c >> 1) & 1);
        coords[2] = (coords[2] << 1) | ((c >> 2) & 1);
    }

    return coords;
}

// ------------------------------------------------------------
// Construct from grid coordinates
// ------------------------------------------------------------
MortonIndex MortonIndex::fromGridCoordinates(
    std::size_t level,
    const Vec<std::size_t, 3>& coords)
{
    if (level > MAX_DEPTH) {
        throw std::invalid_argument("Level exceeds MAX_DEPTH");
    }

    const std::size_t maxCoord = std::size_t{1} << level;
    if (coords[0] >= maxCoord ||
        coords[1] >= maxCoord ||
        coords[2] >= maxCoord)
    {
        throw std::invalid_argument("Grid coordinate out of range");
    }

    std::vector<morton_bits_t> path;
    path.reserve(level);

    for (std::size_t i = 0; i < level; ++i) {
        const std::size_t shift = level - 1 - i;

        const morton_bits_t bx = (coords[0] >> shift) & 1;
        const morton_bits_t by = (coords[1] >> shift) & 1;
        const morton_bits_t bz = (coords[2] >> shift) & 1;

        const morton_bits_t c = (bz << 2) | (by << 1) | bx;
        path.push_back(c);
    }

    return MortonIndex::fromPath(path);
}

} // namespace oktal
