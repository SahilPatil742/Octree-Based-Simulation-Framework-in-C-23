#include "oktal/octree/OctreeGeometry.hpp"

namespace oktal {

// ------------------------------------------------------------
// dx(level) = side / (2^level)
// ------------------------------------------------------------
double OctreeGeometry::dx(size_type level) const noexcept {
    return side_ / static_cast<double>(size_type{1} << level);
}

// ------------------------------------------------------------
// cellExtents(level) = (dx, dx, dx)
// ------------------------------------------------------------
OctreeGeometry::Vec3 OctreeGeometry::cellExtents(size_type level) const noexcept {
    const double h = dx(level);
    return Vec3{h, h, h};
}

// ------------------------------------------------------------
// cellMinCorner: origin + gridCoords * dx
// ------------------------------------------------------------
OctreeGeometry::Vec3 OctreeGeometry::cellMinCorner(MortonIndex m) const noexcept {
    const auto gc = m.gridCoordinates();
    const double h = dx(m.level());

    return Vec3{
        origin_[0] + static_cast<double>(gc[0]) * h,
        origin_[1] + static_cast<double>(gc[1]) * h,
        origin_[2] + static_cast<double>(gc[2]) * h
    };
}

// ------------------------------------------------------------
// cellMaxCorner: minCorner + extents
// ------------------------------------------------------------
OctreeGeometry::Vec3 OctreeGeometry::cellMaxCorner(MortonIndex m) const noexcept {
    const auto min = cellMinCorner(m);
    const double h = dx(m.level());

    return Vec3{
        min[0] + h,
        min[1] + h,
        min[2] + h
    };
}

// ------------------------------------------------------------
// cellBoundingBox: construct a Box<double>
// ------------------------------------------------------------
Box<double> OctreeGeometry::cellBoundingBox(MortonIndex m) const noexcept {
    return Box<double>{cellMinCorner(m), cellMaxCorner(m)};
}

// ------------------------------------------------------------
// cellCenter: midpoint = minCorner + (h/2)
// ------------------------------------------------------------
OctreeGeometry::Vec3 OctreeGeometry::cellCenter(MortonIndex m) const noexcept {
    const auto min = cellMinCorner(m);
    const double h = dx(m.level());
    const double mid = h * 0.5;

    return Vec3{
        min[0] + mid,
        min[1] + mid,
        min[2] + mid
    };
}

} // namespace oktal
