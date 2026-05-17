#pragma once

#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"
#include "oktal/data/GridVector.hpp"

#include <advpt/htgfile/VtkHtgFile.hpp>

#include <filesystem>
#include <string>
#include <string_view>
#include <span>
#include <vector>

namespace oktal::io::vtk {

/**
 * @brief Algorithmic object for exporting CellGrid data to VTK HyperTreeGrid format.
 */
class CellGridVtkWriter {
public:
    CellGridVtkWriter(const CellGrid& grid,
                      const std::filesystem::path& filepath);

    // ============================================================
    // Scalar field export (EXISTING, UNCHANGED)
    // ============================================================
    template <typename T>
    CellGridVtkWriter& writeGridVector(
        std::string_view name,
        std::span<const T> values)
    {
        const auto& octree = grid_->octree();
        std::vector<T> denseData(octree.numberOfNodes(), T{0});

        for (std::size_t i = 0; i < grid_->size(); ++i) {
            if (auto cell = octree.getCell(grid_->mortonIndices()[i])) {
                denseData[cell->streamIndex()] = values[i];
            }
        }

        writer_.writeCellData(std::string(name),
                              std::span<const T>(denseData));
        return *this;
    }

    // ============================================================
    // GridVectorView export (FINAL & CORRECT)
    // ============================================================
    template <typename T, std::size_t Q>
    CellGridVtkWriter& writeGridVector(
        std::string_view name,
        GridVectorView<const T, Q> vector)
    {
        const auto& octree = grid_->octree();
        const std::size_t totalNodes = octree.numberOfNodes();

        // VTK expects array-of-structures (AoS)
        std::vector<T> denseData(totalNodes * Q, T{0});

        for (std::size_t i = 0; i < grid_->size(); ++i) {
            if (auto cell = octree.getCell(grid_->mortonIndices()[i])) {
                const std::size_t cIdx = cell->streamIndex();

                if constexpr (Q == 1) {
                    // rank-1 mdspan
                    denseData[cIdx] = vector[i];
                } else {
                    // rank-2 mdspan → all indices at once
                    for (std::size_t q = 0; q < Q; ++q) {
                        denseData[cIdx * Q + q] = vector[i, q];
                    }
                }
            }
        }

        advpt::htgfile::MdCellDataView<T, 2> mdView(
            denseData.data(), totalNodes, Q);

        writer_.writeCellData(std::string(name), mdView);
        return *this;
    }

private:
    const CellGrid* grid_;
    advpt::htgfile::SnapshotHtgFile writer_;
};

// ============================================================
// Factory helpers
// ============================================================
CellGridVtkWriter exportCellGrid(
    const CellGrid& grid,
    const std::filesystem::path& filepath);

void exportOctree(
    const CellOctree& octree,
    const std::filesystem::path& filepath);

void exportOctree(
    const CellOctree& octree,
    const std::filesystem::path& filepath);

} // namespace oktal::io::vtk
