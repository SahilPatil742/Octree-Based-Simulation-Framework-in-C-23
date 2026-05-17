#include "oktal/octree/CellOctree.hpp"

#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace oktal {

// ============================================================
// Constructors
// ============================================================
CellOctree::CellOctree() {
    nodes_.emplace_back();
    levelStart_.push_back(0);
    levelCount_.push_back(1);
}

CellOctree::CellOctree(const OctreeGeometry& geometry) : geometry_(geometry) {
    nodes_.emplace_back();
    levelStart_.push_back(0);
    levelCount_.push_back(1);
}

// ============================================================
// Basic queries
// ============================================================
std::size_t CellOctree::numberOfNodes() const noexcept {
    return nodes_.size();
}

std::size_t CellOctree::numberOfNodes(std::size_t level) const noexcept {
    if (level >= levelCount_.size()) {
        return 0;
    }
    return levelCount_[level];
}

std::size_t CellOctree::numberOfLevels() const noexcept {
    return levelCount_.size();
}

std::span<const CellOctree::Node> CellOctree::nodesStream() const noexcept {
    return nodes_;
}

std::span<const CellOctree::Node> CellOctree::nodesStream(std::size_t level) const noexcept {
    if (level >= levelStart_.size()) {
        return {};
    }
    return std::span{nodes_}.subspan(levelStart_[level], levelCount_[level]);
}

// ============================================================
// Descriptor factory
// ============================================================
namespace {
void validateDescriptorStructure(const std::vector<std::string_view>& levels) {
    if (levels.empty() || levels.front().size() != 1) {
        throw std::invalid_argument("Invalid descriptor");
    }

    std::size_t expected = 1;
    for (const auto& lvl : levels) {
        if (lvl.size() != expected) {
            throw std::invalid_argument("Invalid descriptor");
        }
        std::size_t refined = 0;
        for (const char c : lvl) {
            if (c == 'R' || c == 'X') {
                refined++;
            } else if (c != '.' && c != 'P') {
                throw std::invalid_argument("Invalid descriptor");
            }
        }
        expected = refined * 8;
    }
}
} // namespace

CellOctree CellOctree::fromDescriptor(std::string_view desc) {
    std::vector<std::string_view> levels;
    std::size_t pos = 0;

    while (true) {
        const auto next = desc.find('|', pos);
        levels.push_back(desc.substr(pos, next - pos));
        if (next == std::string_view::npos) {
            break;
        }
        pos = next + 1;
    }

    validateDescriptorStructure(levels);

    CellOctree tree;
    tree.nodes_.clear();
    tree.levelStart_.clear();
    tree.levelCount_.clear();

    std::size_t global = 0;
    for (const auto& lvl : levels) {
        tree.levelStart_.push_back(global);
        tree.levelCount_.push_back(lvl.size());
        for (const char c : lvl) {
            tree.nodes_.emplace_back(c == 'R' || c == 'X', c == 'P' || c == 'X');
            global++;
        }
    }

    for (std::size_t lvl = 0; lvl + 1 < levels.size(); ++lvl) {
        const std::size_t nextStart = tree.levelStart_[lvl + 1];
        std::size_t cursor = nextStart;
        for (std::size_t i = 0; i < tree.levelCount_[lvl]; ++i) {
            auto& node = tree.nodes_[tree.levelStart_[lvl] + i];
            if (node.isRefined()) {
                node.setChildrenStartIndex(cursor);
                cursor += 8;
            }
        }
    }
    return tree;
}

// ============================================================
// Factory Utilities
// ============================================================
std::shared_ptr<const CellOctree> CellOctree::createUniformGrid(std::size_t level) {
    return createUniformGrid(OctreeGeometry(), level);
}

std::shared_ptr<const CellOctree> CellOctree::createUniformGrid(OctreeGeometry geom, std::size_t level) {
    auto tree = std::make_shared<CellOctree>(geom);

    for (std::size_t l = 0; l < level; ++l) {
        const std::size_t currentLevelStart = tree->levelStart_[l];
        const std::size_t currentLevelCount = tree->levelCount_[l];
        const std::size_t nextLevelStart = tree->nodes_.size();

        for (std::size_t i = 0; i < currentLevelCount; ++i) {
            const std::size_t nodeIdx = currentLevelStart + i;
            tree->nodes_[nodeIdx].setRefined(true);
            tree->nodes_[nodeIdx].setPhantom(true);
            tree->nodes_[nodeIdx].setChildrenStartIndex(tree->nodes_.size());

            for (int j = 0; j < 8; ++j) {
                tree->nodes_.emplace_back(false, false, 0);
            }
        }
        tree->levelStart_.push_back(nextLevelStart);
        tree->levelCount_.push_back(tree->nodes_.size() - nextLevelStart);
    }
    return tree;
}

// ============================================================
// Cell queries
// ============================================================
std::optional<CellOctree::CellView> CellOctree::getCell(MortonIndex index) const {
    const std::size_t lvl = index.level();
    if (lvl >= levelCount_.size()) {
        return std::nullopt;
    }

    std::size_t idx = 0;
    for (std::size_t l = 0; l < lvl; ++l) {
        const std::size_t shift = 3 * (lvl - l - 1);
        const std::size_t branch = (index.getBits() >> shift) & 0b111;
        const Node& parent = nodes_[idx];
        if (!parent.isRefined()) {
            return std::nullopt;
        }
        idx = parent.childIndex(branch);
    }

    const Node& node = nodes_[idx];
    if (node.isPhantom()) {
        return std::nullopt;
    }
    return CellView{this, index, &nodes_[idx], idx};
}

bool CellOctree::cellExists(MortonIndex index) const {
    return getCell(index).has_value();
}

std::optional<CellOctree::CellView> CellOctree::getRootCell() const {
    if (nodes_.empty() || nodes_.front().isPhantom()) {
        return std::nullopt;
    }
    // Fixed: Using .data() per readability-container-data-pointer
    return CellView{this, MortonIndex{}, nodes_.data(), 0};
}

// ============================================================
// CellView geometry
// ============================================================
Vec<double, 3> CellOctree::CellView::center() const {
    return tree_->geometry_.cellCenter(morton_);
}

Box<double> CellOctree::CellView::boundingBox() const {
    return tree_->geometry_.cellBoundingBox(morton_);
}

} // namespace oktal