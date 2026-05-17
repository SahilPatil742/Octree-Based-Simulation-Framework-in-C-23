#include "oktal/octree/CellGrid.hpp"

#include <algorithm>
#include <ranges>
#include <map>
#include <stdexcept>

namespace oktal {

/* ---------- Torus Periodicity ---------- */

Torus::Torus(std::initializer_list<bool> periodic) {
    const std::span<const bool> periodic_span{periodic};

    for (std::size_t i = 0; i < 3; ++i) {
        if (i < periodic_span.size()) {
            m_periodic.at(i) = periodic_span[i];
        } else {
            m_periodic.at(i) = false;
        }
    }
}

bool Torus::operator()(Vec<long, 3>& coords,
                       const Vec<long, 3>& domainSize) const {
    for (std::size_t i = 0; i < 3; ++i) {
        if (m_periodic.at(i)) {
            coords[i] =
                ((coords[i] % domainSize[i]) + domainSize[i]) % domainSize[i];
        } else {
            if (coords[i] < 0 || coords[i] >= domainSize[i]) {
                return false;
            }
        }
    }
    return true;
}

/* ---------- CellGridBuilder ---------- */

CellGridBuilder::CellGridBuilder(std::shared_ptr<const CellOctree> octree)
    : m_octree(std::move(octree)) {}

CellGridBuilder&
CellGridBuilder::levels(std::span<const std::size_t> levels) {
    m_levels.assign(levels.begin(), levels.end());
    return *this;
}

CellGridBuilder&
CellGridBuilder::levels(std::initializer_list<std::size_t> levels) {
    m_levels.assign(levels.begin(), levels.end());
    return *this;
}

CellGridBuilder&
CellGridBuilder::neighborhood(
    std::span<const Vec<std::ptrdiff_t, 3>> offsets) {
    m_neighborhood.assign(offsets.begin(), offsets.end());
    return *this;
}

CellGridBuilder&
CellGridBuilder::neighborhood(
    std::initializer_list<Vec<std::ptrdiff_t, 3>> offsets) {
    m_neighborhood = offsets;
    return *this;
}

/* ---------- Helper Member Functions ---------- */

std::vector<std::size_t>
CellGridBuilder::collectTargetLevels() const {
    if (!m_levels.empty()) {
        return m_levels;
    }

    const auto nLevels = m_octree->numberOfLevels();
    std::vector<std::size_t> levels;
    levels.reserve(nLevels);  // ✅ clang-tidy fix

    for (std::size_t l = 0; l < nLevels; ++l) {
        levels.push_back(l);
    }
    return levels;
}


void CellGridBuilder::enumerateCells(
    CellGrid& grid,
    const std::vector<std::size_t>& levels) const {

    for (const std::size_t lvl : levels) {
        if (lvl >= m_octree->numberOfLevels()) {
            continue;
        }

        for (const auto cell : m_octree->horizontalRange(lvl)) {
            if (cell.level() != lvl) {
                continue;
            }

            if (grid.getEnumerationIndex(cell.streamIndex()) !=
                CellGrid::NOT_ENUMERATED) {
                continue;
            }

            const auto enumIdx = grid.m_mortonIndices.size();
            grid.m_mortonIndices.push_back(cell.mortonIndex());
            grid.m_octreeCells.push_back(cell);
            grid.m_streamToEnum.emplace_back(
                cell.streamIndex(), enumIdx);
        }
    }

    std::ranges::sort(
        grid.m_streamToEnum,
        [](const auto& a, const auto& b) {
            return a.first < b.first;
        });
}

void CellGridBuilder::buildAdjacency(CellGrid& grid) const {
    if (m_neighborhood.empty()) {
        return;
    }

    const auto periodicity =
        m_periodicity ? m_periodicity
                      : std::make_shared<NoPeriodicity>();

    grid.m_neighborOffsets = m_neighborhood;
    grid.m_adjacencyLists.resize(m_neighborhood.size());

    for (auto& list : grid.m_adjacencyLists) {
        list.assign(
            grid.m_mortonIndices.size(),
            CellGrid::NO_NEIGHBOR);
    }

    for (std::size_t i = 0; i < grid.m_mortonIndices.size(); ++i) {
        const auto& cell = grid.m_octreeCells[i];
        const auto lvl = cell.level();
        const auto coords =
            cell.mortonIndex().gridCoordinates();
        const long dim = 1L << lvl;
        const Vec<long, 3> domain{dim, dim, dim};

        for (std::size_t n = 0; n < m_neighborhood.size(); ++n) {
            Vec<long, 3> target{
                static_cast<long>(coords[0]) +
                    static_cast<long>(m_neighborhood[n][0]),
                static_cast<long>(coords[1]) +
                    static_cast<long>(m_neighborhood[n][1]),
                static_cast<long>(coords[2]) +
                    static_cast<long>(m_neighborhood[n][2])};

            if (!(*periodicity)(target, domain)) {
                continue;
            }

            const Vec<std::size_t, 3> unsignedTarget{
                static_cast<std::size_t>(target[0]),
                static_cast<std::size_t>(target[1]),
                static_cast<std::size_t>(target[2])};

            const auto neighborMorton =
                MortonIndex::fromGridCoordinates(
                    lvl, unsignedTarget);

            if (const auto neighbor =
                    m_octree->getCell(neighborMorton)) {
                grid.m_adjacencyLists[n][i] =
                    grid.getEnumerationIndex(*neighbor);
            }
        }
    }
}

/* ---------- build() ---------- */

CellGrid CellGridBuilder::build() {
    CellGrid grid;
    grid.m_octree = m_octree;

    const auto levels = collectTargetLevels();
    enumerateCells(grid, levels);
    buildAdjacency(grid);

    return grid;
}

/* ---------- CellGrid Queries ---------- */

CellGridBuilder
CellGrid::create(std::shared_ptr<const CellOctree> octree) {
    return CellGridBuilder{std::move(octree)};
}

std::size_t
CellGrid::getEnumerationIndex(std::size_t streamIndex) const {
    auto it = std::ranges::lower_bound(
        m_streamToEnum,
        streamIndex,
        {},
        &std::pair<std::size_t, std::size_t>::first);

    if (it != m_streamToEnum.end() &&
        it->first == streamIndex) {
        return it->second;
    }
    return NOT_ENUMERATED;
}

std::size_t
CellGrid::getEnumerationIndex(
    const CellOctree::CellView& cell) const {
    return getEnumerationIndex(cell.streamIndex());
}

std::span<const std::size_t>
CellGrid::neighborIndices(
    Vec<std::ptrdiff_t, 3> offset) const {
    auto it = std::ranges::find(m_neighborOffsets, offset);
    if (it == m_neighborOffsets.end()) {
        throw std::out_of_range(
            "Neighborhood offset not found.");
    }

    const auto index = static_cast<std::size_t>(
        std::distance(m_neighborOffsets.begin(), it));

    return m_adjacencyLists[index];
}

/* ---------- CellView ---------- */

MortonIndex CellGrid::CellView::mortonIndex() const {
    return isValid()
        ? m_grid->m_mortonIndices[m_index]
        : MortonIndex{0};
}

Vec<double, 3> CellGrid::CellView::center() const {
    return isValid()
        ? m_grid->m_octreeCells[m_index].center()
        : Vec<double, 3>{0, 0, 0};
}

std::size_t CellGrid::CellView::level() const {
    return isValid()
        ? m_grid->m_octreeCells[m_index].level()
        : 0;
}

Box<double>
CellGrid::CellView::boundingBox() const {
    return isValid()
        ? m_grid->m_octreeCells[m_index].boundingBox()
        : Box<double>{{0, 0, 0}, {0, 0, 0}};
}

CellGrid::CellView
CellGrid::CellView::neighbor(
    Vec<std::ptrdiff_t, 3> offset) const {
    if (!isValid()) {
        return CellView{m_grid, NO_NEIGHBOR};
    }

    auto it = std::ranges::find(
        m_grid->m_neighborOffsets, offset);
    if (it == m_grid->m_neighborOffsets.end()) {
        return CellView{m_grid, NO_NEIGHBOR};
    }

    const auto dir = static_cast<std::size_t>(
        std::distance(
            m_grid->m_neighborOffsets.begin(), it));

    return CellView{
        m_grid,
        m_grid->m_adjacencyLists[dir][m_index]};
}

} // namespace oktal
