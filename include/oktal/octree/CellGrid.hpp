#pragma once

#include <vector>
#include <memory>
#include <span>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <algorithm>
#include <iterator>
#include <array>

#include "oktal/geometry/Vec.hpp"
#include "oktal/geometry/Box.hpp"
#include "oktal/octree/CellOctree.hpp"
#include "oktal/octree/MortonIndex.hpp"

namespace oktal {

class CellGrid;

/* ---------------- Periodicity ---------------- */

class PeriodicityMapper {
public:
    PeriodicityMapper() = default;
    PeriodicityMapper(const PeriodicityMapper&) = default;
    PeriodicityMapper& operator=(const PeriodicityMapper&) = default;
    PeriodicityMapper(PeriodicityMapper&&) = default;
    PeriodicityMapper& operator=(PeriodicityMapper&&) = default;
    virtual ~PeriodicityMapper() = default;

    virtual bool operator()(Vec<long, 3>& coords,
                            const Vec<long, 3>& domainSize) const = 0;
};

class NoPeriodicity final : public PeriodicityMapper {
public:
    bool operator()(Vec<long, 3>& coords,
                    const Vec<long, 3>& domainSize) const override {
        for (std::size_t i = 0; i < 3; ++i) {
            if (coords[i] < 0 || coords[i] >= domainSize[i]) {
                return false;
            }
        }
        return true;
    }
};

class Torus final : public PeriodicityMapper {
public:
    explicit Torus(std::initializer_list<bool> periodic);

    bool operator()(Vec<long, 3>& coords,
                    const Vec<long, 3>& domainSize) const override;

private:
    std::array<bool, 3> m_periodic{};
};

/* ---------------- Builder ---------------- */

class CellGridBuilder {
    std::shared_ptr<const CellOctree> m_octree;
    std::vector<std::size_t> m_levels;
    std::vector<Vec<std::ptrdiff_t, 3>> m_neighborhood;
    std::shared_ptr<PeriodicityMapper> m_periodicity;

public:
    explicit CellGridBuilder(std::shared_ptr<const CellOctree> octree);

    CellGridBuilder& levels(std::span<const std::size_t> levels);
    CellGridBuilder& levels(std::initializer_list<std::size_t> levels);

    CellGridBuilder& neighborhood(
        std::span<const Vec<std::ptrdiff_t, 3>> offsets);
    CellGridBuilder& neighborhood(
        std::initializer_list<Vec<std::ptrdiff_t, 3>> offsets);

    template <typename T>
    CellGridBuilder& periodicityMapper(T mapper) {
        m_periodicity = std::make_shared<T>(std::move(mapper));
        return *this;
    }

    CellGrid build();

private:
    // --- helper functions to keep build() simple ---
    std::vector<std::size_t> collectTargetLevels() const;

    void enumerateCells(
        CellGrid& grid,
        const std::vector<std::size_t>& levels) const;

    void buildAdjacency(CellGrid& grid) const;
};

/* ---------------- CellGrid ---------------- */

class CellGrid {
public:
    static constexpr std::size_t NOT_ENUMERATED =
        std::numeric_limits<std::size_t>::max();
    static constexpr std::size_t NO_NEIGHBOR = NOT_ENUMERATED;

    class CellView;

    class Iterator {
        const CellGrid* m_grid = nullptr;
        std::size_t m_current = 0;

    public:
        using iterator_category = std::forward_iterator_tag;
        using iterator_concept  = std::forward_iterator_tag;
        using value_type        = CellView;
        using difference_type   = std::ptrdiff_t;

        Iterator() = default;
        Iterator(const CellGrid* grid, std::size_t index)
            : m_grid(grid), m_current(index) {}

        [[nodiscard]] CellView operator*() const;
        Iterator& operator++() { ++m_current; return *this; }
        Iterator operator++(int) { Iterator tmp{*this}; ++(*this); return tmp; }

        bool operator==(const Iterator& other) const {
            return m_current == other.m_current;
        }
        bool operator!=(const Iterator& other) const {
            return !(*this == other);
        }
    };

    class CellView {
        const CellGrid* m_grid = nullptr;
        std::size_t m_index = NOT_ENUMERATED;

    public:
        CellView() = default;
        CellView(const CellGrid* grid, std::size_t index)
            : m_grid(grid), m_index(index) {}

        [[nodiscard]] std::size_t enumerationIndex() const { return m_index; }

        // REQUIRED by tests (implicit!)
        operator std::size_t() const { return m_index; }
        operator bool() const { return m_index != NOT_ENUMERATED; }

        [[nodiscard]] bool isValid() const {
            return m_index != NOT_ENUMERATED;
        }

        [[nodiscard]] Vec<double, 3> center() const;
        [[nodiscard]] std::size_t level() const;
        [[nodiscard]] Box<double> boundingBox() const;
        [[nodiscard]] MortonIndex mortonIndex() const;

        [[nodiscard]] CellView neighbor(Vec<std::ptrdiff_t, 3> offset) const;
    };

    static CellGridBuilder create(std::shared_ptr<const CellOctree> octree);

    [[nodiscard]] const CellOctree& octree() const { return *m_octree; }
    [[nodiscard]] std::span<const MortonIndex> mortonIndices() const {
        return m_mortonIndices;
    }

    [[nodiscard]] std::size_t getEnumerationIndex(
        std::size_t streamIndex) const;
    [[nodiscard]] std::size_t getEnumerationIndex(
        const CellOctree::CellView& cell) const;

    [[nodiscard]] std::span<const std::size_t> neighborIndices(
        Vec<std::ptrdiff_t, 3> offset) const;

    [[nodiscard]] Iterator begin() const { return Iterator{this, 0}; }
    [[nodiscard]] Iterator end() const {
        return Iterator{this, m_mortonIndices.size()};
    }
    [[nodiscard]] std::size_t size() const { return m_mortonIndices.size(); }

    [[nodiscard]] CellView operator[](std::size_t index) const {
        return CellView{this, index};
    }

private:
    friend class CellGridBuilder;
    CellGrid() = default;

    std::shared_ptr<const CellOctree> m_octree;
    std::vector<MortonIndex> m_mortonIndices;
    std::vector<std::pair<std::size_t, std::size_t>> m_streamToEnum;
    std::vector<Vec<std::ptrdiff_t, 3>> m_neighborOffsets;
    std::vector<std::vector<std::size_t>> m_adjacencyLists;
    std::vector<CellOctree::CellView> m_octreeCells;
};

inline CellGrid::CellView CellGrid::Iterator::operator*() const {
    return CellView{m_grid, m_current};
}

} // namespace oktal
