#pragma once

#include <algorithm>
#include <concepts>
#include <cstddef>
#include <cstdint>
#include <iterator>
#include <memory>
#include <optional>
#include <span>
#include <stdexcept>
#include <string_view>
#include <utility>
#include <vector>

#include "oktal/geometry/Box.hpp"
#include "oktal/geometry/Vec.hpp"
#include "oktal/octree/MortonIndex.hpp"
#include "oktal/octree/OctreeGeometry.hpp"

namespace oktal {

class OctreeCursor;

class CellOctree {
public:
    class Node {
    public:
        constexpr Node() noexcept = default;
        constexpr Node(bool refined, bool phantom, std::size_t childrenIdx = 0) noexcept
            : refined_{refined}, phantom_{phantom}, childrenStart_{childrenIdx} {}

        [[nodiscard]] constexpr bool isRefined() const noexcept { return refined_; }
        [[nodiscard]] constexpr bool isPhantom() const noexcept { return phantom_; }
        constexpr void setRefined(bool v) noexcept { refined_ = v; }
        constexpr void setPhantom(bool v) noexcept { phantom_ = v; }
        [[nodiscard]] constexpr std::size_t childrenStartIndex() const noexcept { return childrenStart_; }
        constexpr void setChildrenStartIndex(std::size_t idx) noexcept { childrenStart_ = idx; }
        [[nodiscard]] constexpr std::size_t childIndex(std::size_t branch) const noexcept {
            return childrenStart_ + branch;
        }

    private:
        bool refined_{false};
        bool phantom_{false};
        std::size_t childrenStart_{0};
    };

    class CellView {
    public:
        [[nodiscard]] MortonIndex mortonIndex() const noexcept { return morton_; }
        [[nodiscard]] bool isRoot() const noexcept { return morton_.isRoot(); }
        [[nodiscard]] bool isRefined() const noexcept { return node_->isRefined(); }
        [[nodiscard]] std::size_t level() const noexcept { return morton_.level(); }
        [[nodiscard]] std::size_t streamIndex() const noexcept { return index_; }
        [[nodiscard]] Vec<double, 3> center() const;
        [[nodiscard]] Box<double> boundingBox() const;

    private:
        friend class CellOctree;
        friend class OctreeCursor;
        CellView(const CellOctree* tree, MortonIndex morton, const Node* node, std::size_t index) noexcept
            : tree_{tree}, morton_{morton}, node_{node}, index_{index} {}

        const CellOctree* tree_;
        MortonIndex morton_;
        const Node* node_;
        std::size_t index_;
    };

    CellOctree();
    explicit CellOctree(const OctreeGeometry& geometry);

    [[nodiscard]] static std::shared_ptr<const CellOctree> createUniformGrid(std::size_t level);
    [[nodiscard]] static std::shared_ptr<const CellOctree> createUniformGrid(OctreeGeometry geom, std::size_t level);

    [[nodiscard]] const OctreeGeometry& geometry() const noexcept { return geometry_; }
    [[nodiscard]] std::size_t numberOfNodes() const noexcept;
    [[nodiscard]] std::size_t numberOfNodes(std::size_t level) const noexcept;
    [[nodiscard]] std::size_t numberOfLevels() const noexcept;
    [[nodiscard]] std::span<const Node> nodesStream() const noexcept;
    [[nodiscard]] std::span<const Node> nodesStream(std::size_t level) const noexcept;
    [[nodiscard]] static CellOctree fromDescriptor(std::string_view descriptor);
    [[nodiscard]] std::optional<CellView> getCell(MortonIndex index) const;
    [[nodiscard]] bool cellExists(MortonIndex index) const;
    [[nodiscard]] std::optional<CellView> getRootCell() const;

    [[nodiscard]] auto preOrderDepthFirstRange() const;
    [[nodiscard]] auto horizontalRange(std::size_t level) const;

private:
    OctreeGeometry geometry_;
    std::vector<Node> nodes_;
    std::vector<std::size_t> levelStart_;
    std::vector<std::size_t> levelCount_;
};

class OctreeCursor {
public:
    OctreeCursor() noexcept = default;
    explicit OctreeCursor(const CellOctree& octree) : tree_(&octree) {
        if (!octree.nodesStream().empty()) {
            path_.reserve(16);
            path_.push_back(0);
        }
    }
    OctreeCursor(const CellOctree& octree, std::span<const std::size_t> path)
        : tree_(&octree), path_(path.begin(), path.end()) {}

    [[nodiscard]] const CellOctree* octree() const noexcept { return tree_; }
    [[nodiscard]] std::span<const std::size_t> path() const noexcept { return path_; }
    [[nodiscard]] bool empty() const noexcept { return tree_ == nullptr; }
    [[nodiscard]] bool end() const noexcept { return tree_ != nullptr && path_.empty(); }
    [[nodiscard]] std::size_t currentLevel() const noexcept { return path_.empty() ? 0 : path_.size() - 1; }
    [[nodiscard]] std::size_t currentStreamIndex() const { return path_.empty() ? 0 : path_.back(); }

    [[nodiscard]] MortonIndex mortonIndex() const {
        MortonIndex m;
        for (std::size_t i = 1; i < path_.size(); ++i) {
            const std::size_t parentStreamIdx = path_[i - 1];
            const std::size_t currentStreamIdx = path_[i];
            const auto& parentNode = tree_->nodesStream()[parentStreamIdx];
            const std::size_t octant = currentStreamIdx - parentNode.childrenStartIndex();
            m = m.child(octant);
        }
        return m;
    }

    [[nodiscard]] std::optional<CellOctree::CellView> currentCell() const {
        if (end() || empty()) {
            return std::nullopt;
        }
        const std::size_t idx = currentStreamIndex();
        if (idx >= tree_->nodesStream().size()) {
            return std::nullopt;
        }
        const auto& node = tree_->nodesStream()[idx];
        if (node.isPhantom()) {
            return std::nullopt;
        }
        return CellOctree::CellView(tree_, mortonIndex(), &node, idx);
    }

    [[nodiscard]] bool firstSibling() const {
        if (path_.size() <= 1) {
            return true;
        }
        const auto& parent = tree_->nodesStream()[path_[path_.size() - 2]];
        return currentStreamIndex() == parent.childrenStartIndex();
    }

    [[nodiscard]] bool lastSibling() const {
        if (path_.size() <= 1) {
            return true;
        }
        const auto& parent = tree_->nodesStream()[path_[path_.size() - 2]];
        return currentStreamIndex() == parent.childrenStartIndex() + 7;
    }

    bool operator==(const OctreeCursor& other) const noexcept {
        if (tree_ != other.tree_) {
            return false;
        }
        if (tree_ == nullptr) {
            return true;
        }
        return path_ == other.path_;
    }
    bool operator!=(const OctreeCursor& other) const noexcept { return !(*this == other); }

    void ascend() {
        if (!path_.empty()) {
            path_.pop_back();
        }
    }
    void descend() {
        if (path_.empty()) {
            return;
        }
        const auto& node = tree_->nodesStream()[currentStreamIndex()];
        if (node.isRefined()) {
            path_.push_back(node.childrenStartIndex());
        }
    }

    void descend(std::size_t childIdx) {
        if (path_.empty()) {
            return;
        }
        if (childIdx > 7) {
            throw std::out_of_range("Child index out of range (0-7)");
        }
        const auto& node = tree_->nodesStream()[currentStreamIndex()];
        if (node.isRefined()) {
            path_.push_back(node.childrenStartIndex() + childIdx);
        }
    }

    void previousSibling() {
        if (path_.size() <= 1) {
            return;
        }
        const auto& parent = tree_->nodesStream()[path_[path_.size() - 2]];
        if (currentStreamIndex() > parent.childrenStartIndex()) {
            path_.back()--;
        }
    }

    void nextSibling() {
        if (path_.size() <= 1) {
            return;
        }
        const auto& parent = tree_->nodesStream()[path_[path_.size() - 2]];
        if (currentStreamIndex() < parent.childrenStartIndex() + 7) {
            path_.back()++;
        }
    }

    void toSibling(std::size_t siblingIdx) {
        if (path_.size() <= 1) {
            if (siblingIdx != 0) {
                throw std::out_of_range("Root has no siblings");
            }
            return;
        }
        if (siblingIdx > 7) {
            throw std::out_of_range("Sibling index out of range (0-7)");
        }
        const auto& parent = tree_->nodesStream()[path_[path_.size() - 2]];
        path_.back() = parent.childrenStartIndex() + siblingIdx;
    }

    void toEnd() { path_.clear(); }

private:
    const CellOctree* tree_ = nullptr;
    std::vector<std::size_t> path_;
};

inline bool advanceCursorLogic(OctreeCursor& c) {
    if (!c.lastSibling()) {
        c.nextSibling();
    } else {
        while (!c.end() && c.lastSibling()) {
            c.ascend();
        }
        if (c.end()) {
            return false;
        }
        c.nextSibling();
    }
    return true;
}

template <typename T>
concept OctreeIteratorPolicy = std::semiregular<T> && requires(const T t, OctreeCursor& c) {
    { t.advance(c) } -> std::same_as<void>;
};

struct DFSPolicy {
    static void advance(OctreeCursor& c) {
        if (c.empty() || c.end()) {
            return;
        }
        do {
            const auto& node = c.octree()->nodesStream()[c.currentStreamIndex()];
            if (node.isRefined()) {
                c.descend();
            } else {
                while (!c.end() && c.lastSibling()) {
                    c.ascend();
                }
                if (c.end()) {
                    return;
                }
                c.nextSibling();
            }
        } while (!c.end() && !c.currentCell().has_value());
    }
};

class HorizontalPolicy {
public:
    constexpr HorizontalPolicy() noexcept = default;
    explicit constexpr HorizontalPolicy(std::size_t level) noexcept : level_(level) {}

    void advance(OctreeCursor& c) const {
        if (c.empty() || c.end()) {
            return;
        }
        do {
            if (c.currentLevel() == level_) {
                if (!advanceCursorLogic(c)) {
                    return;
                }
            }
            while (!c.end() && c.currentLevel() < level_) {
                const auto& node = c.octree()->nodesStream()[c.currentStreamIndex()];
                if (node.isRefined()) {
                    c.descend();
                } else {
                    if (!advanceCursorLogic(c)) {
                        return;
                    }
                }
            }
        } while (!c.end() && (c.currentLevel() != level_ || !c.currentCell().has_value()));
    }

private:
    std::size_t level_ = 0;
};

template <OctreeIteratorPolicy TPolicy>
class OctreeCellsIterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = CellOctree::CellView;
    using difference_type = std::ptrdiff_t;
    using pointer = void;
    using reference = CellOctree::CellView;

    OctreeCellsIterator() = default;
    OctreeCellsIterator(OctreeCursor cursor, TPolicy policy)
        : cursor_(std::move(cursor)), policy_(std::move(policy)) {}

    [[nodiscard]] reference operator*() const {
        const auto opt = cursor_.currentCell();
        if (!opt.has_value()) {
            throw std::runtime_error("Iterator dereference failed");
        }
        return *opt;
    }

    OctreeCellsIterator& operator++() {
        policy_.advance(cursor_);
        return *this;
    }
    OctreeCellsIterator operator++(int) {
        auto tmp = *this;
        ++(*this);
        return tmp;
    }
    [[nodiscard]] bool operator==(const OctreeCellsIterator& other) const noexcept { return cursor_ == other.cursor_; }
    [[nodiscard]] bool operator!=(const OctreeCellsIterator& other) const noexcept { return !(*this == other); }

private:
    OctreeCursor cursor_;
    TPolicy policy_;
};

template <OctreeIteratorPolicy TPolicy>
using OctreeIterator = OctreeCellsIterator<TPolicy>;

template <OctreeIteratorPolicy TPolicy>
class OctreeCellsRange {
public:
    OctreeCellsRange(OctreeCursor start, OctreeCursor end, TPolicy policy)
        : start_(std::move(start)), end_(std::move(end)), policy_(std::move(policy)) {}

    [[nodiscard]] OctreeCellsIterator<TPolicy> begin() const { return OctreeCellsIterator<TPolicy>(start_, policy_); }
    [[nodiscard]] OctreeCellsIterator<TPolicy> end() const { return OctreeCellsIterator<TPolicy>(end_, policy_); }

private:
    OctreeCursor start_;
    OctreeCursor end_;
    TPolicy policy_;
};

inline auto CellOctree::preOrderDepthFirstRange() const {
    OctreeCursor start(*this);
    OctreeCursor endCursor(*this);
    endCursor.toEnd();
    if (!start.end() && !start.currentCell().has_value()) {
        DFSPolicy::advance(start);
    }
    return OctreeCellsRange<DFSPolicy>(std::move(start), std::move(endCursor), DFSPolicy{});
}

inline auto CellOctree::horizontalRange(std::size_t level) const {
    OctreeCursor start(*this);
    OctreeCursor endCursor(*this);
    endCursor.toEnd();
    const HorizontalPolicy policy(level);

    if (level > 0) {
        while (!start.end() && start.currentLevel() < level) {
            const auto& node = nodesStream()[start.currentStreamIndex()];
            if (node.isRefined()) {
                start.descend();
            } else {
                if (!advanceCursorLogic(start)) {
                    break;
                }
            }
        }
    }

    if (start.currentLevel() != level) {
        start.toEnd();
    } else if (!start.currentCell().has_value()) {
        policy.advance(start);
    }
    return OctreeCellsRange<HorizontalPolicy>(std::move(start), std::move(endCursor), policy);
}

} // namespace oktal