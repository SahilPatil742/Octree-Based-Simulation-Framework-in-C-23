#include "oktal/io/VtkExport.hpp"

#include <advpt/htgfile/VtkHtgFile.hpp>
#include <cstdint>
#include <vector>

namespace oktal::io::vtk {

namespace {

// Prepare HyperTreeGrid structure
advpt::htgfile::HyperTree prepareHyperTree(const CellOctree& octree)
{
    using advpt::htgfile::HyperTree;
    HyperTree tree;

    const auto& geom = octree.geometry();
    const auto origin = geom.origin();
    const double side = geom.sidelength();

    tree.xCoords = {origin[0], origin[0] + side};
    tree.yCoords = {origin[1], origin[1] + side};
    tree.zCoords = {origin[2], origin[2] + side};

    const std::size_t levels = octree.numberOfLevels();
    tree.nodesPerDepth.resize(levels);
    for (std::size_t l = 0; l < levels; ++l) {
        tree.nodesPerDepth[l] =
            static_cast<int64_t>(octree.numberOfNodes(l));
    }

    const std::size_t totalNodes = octree.numberOfNodes();
    const std::size_t lastLevelNodes = octree.numberOfNodes(levels - 1);
    const std::size_t nonLeafNodes = totalNodes - lastLevelNodes;

    tree.descriptor.assign((nonLeafNodes + 7) / 8, uint8_t{0});

    std::size_t bit = 0;
    for (std::size_t l = 0; l + 1 < levels; ++l) {
        for (const auto& node : octree.nodesStream(l)) {
            if (node.isRefined()) {
                tree.descriptor[bit / 8] |=
                    static_cast<uint8_t>(1u << (7 - (bit % 8)));
            }
            ++bit;
        }
    }

    tree.mask.assign((totalNodes + 7) / 8, uint8_t{0});
    bit = 0;
    for (const auto& node : octree.nodesStream()) {
        if (node.isPhantom() && !node.isRefined()) {
            tree.mask[bit / 8] |=
                static_cast<uint8_t>(1u << (7 - (bit % 8)));
        }
        ++bit;
    }

    return tree;
}

} // namespace

// ============================================================
// CellGridVtkWriter implementation
// ============================================================
CellGridVtkWriter::CellGridVtkWriter(
    const CellGrid& grid,
    const std::filesystem::path& filepath)
    : grid_(&grid),
      writer_(advpt::htgfile::SnapshotHtgFile::create(
          filepath, prepareHyperTree(grid.octree())))
{
    const auto& octree = grid_->octree();
    const std::size_t totalNodes = octree.numberOfNodes();

    // Write refinement level as default cell data
    std::vector<std::size_t> levelData(totalNodes);

    for (std::size_t l = 0; l < octree.numberOfLevels(); ++l) {
        const auto stream = octree.nodesStream(l);
        if (stream.empty()) {
            continue;
        }

        const std::size_t offset =
            static_cast<std::size_t>(
                stream.data() - octree.nodesStream().data());

        for (std::size_t i = 0; i < stream.size(); ++i) {
            levelData[offset + i] = l;
        }
    }

    writer_.writeCellData("level", levelData);
}

// ============================================================
// Free helper functions
// ============================================================
CellGridVtkWriter exportCellGrid(
    const CellGrid& grid,
    const std::filesystem::path& filepath)
{
    return {grid, filepath};
}

void exportOctree(
    const CellOctree& octree,
    const std::filesystem::path& filepath)
{
    const auto tree = prepareHyperTree(octree);
    auto file = advpt::htgfile::SnapshotHtgFile::create(filepath, tree);

    std::vector<std::size_t> levelData(octree.numberOfNodes());
    for (std::size_t l = 0; l < octree.numberOfLevels(); ++l) {
        const auto stream = octree.nodesStream(l);
        if (stream.empty()) {
            continue;
        }

        const std::size_t offset =
            static_cast<std::size_t>(
                stream.data() - octree.nodesStream().data());

        for (std::size_t i = 0; i < stream.size(); ++i) {
            levelData[offset + i] = l;
        }
    }

    file.writeCellData("level", levelData);
}

} // namespace oktal::io::vtk
