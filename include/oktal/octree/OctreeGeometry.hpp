#pragma once

#include <cstddef>
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include "oktal/geometry/Box.hpp"

namespace oktal {

class OctreeGeometry {
public:
  using Vec3 = Vec<double, 3>;
  using size_type = std::size_t;

private:
  Vec3 origin_;     // bottom-south-west corner
  double side_;     // side length of the entire cube

public:
  // Default ctor -> unit cube at origin
  constexpr OctreeGeometry() noexcept
    : origin_{0.0}, side_{1.0} {}

  // Origin + side length ctor
  constexpr OctreeGeometry(const Vec3 &origin, double side) noexcept
    : origin_{origin}, side_{side} {}

  // ---------------- Accessors ----------------

  [[nodiscard]] constexpr const Vec3 &origin() const noexcept {
    return origin_;
  }

  [[nodiscard]] constexpr double sidelength() const noexcept {
    return side_;
  }

  // ---------------- Cell extents per level ----------------

  // dx(level) = side length / 2^level
  [[nodiscard]] double dx(size_type level) const noexcept;

  // cellExtents(level) = vector (dx, dx, dx)
  [[nodiscard]] Vec3 cellExtents(size_type level) const noexcept;

  // ---------------- Cell geometry ----------------

  [[nodiscard]] Vec3 cellMinCorner(MortonIndex m) const noexcept;
  [[nodiscard]] Vec3 cellMaxCorner(MortonIndex m) const noexcept;
  [[nodiscard]] Box<double> cellBoundingBox(MortonIndex m) const noexcept;
  [[nodiscard]] Vec3 cellCenter(MortonIndex m) const noexcept;
};

} // namespace oktal
