#include <cmath>
#include <iostream>
#include <memory>
#include <span>
#include <string>
#include <vector>
#include <stdexcept>

#include "oktal/io/VtkExport.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/octree/CellOctree.hpp"

using namespace oktal;

namespace { // Internal linkage for misc-use-internal-linkage

/**
 * @brief Jacobi Solver extracted to keep main's cognitive complexity below 25.
 */
void solvePoisson(const CellGrid& grid,
                  std::vector<double>& u,
                  const std::vector<double>& f,
                  const std::vector<bool>& isBoundary,
                  const std::vector<Vec<std::ptrdiff_t, 3>>& offsets,
                  double h2,
                  std::size_t maxIter,
                  double epsilon) {
    const std::size_t nCells = grid.size();
    std::vector<double> u_next = u;
    std::vector<double> residual(nCells, 0.0);

    double l2_norm = epsilon + 1.0;
    std::size_t iter = 0;

    while (iter < maxIter && l2_norm > epsilon) {
        // Step 1: Update values for non-boundary cells
        for (std::size_t i = 0; i < nCells; ++i) {
            if (isBoundary.at(i)) {
                continue;
            }

            double neighbors_sum = 0.0;
            for (const auto& off : offsets) {
                const auto nb = grid[i].neighbor(off);
                neighbors_sum += u.at(nb.enumerationIndex());
            }
            u_next.at(i) = (1.0 / 6.0) * (h2 * f.at(i) + neighbors_sum);
        }
        u = u_next;

        // Step 2: Compute L2 norm of the residual
        double sum_r2 = 0.0;
        for (std::size_t i = 0; i < nCells; ++i) {
            if (isBoundary.at(i)) {
                residual.at(i) = 0.0;
            } else {
                double neighbors_sum = 0.0;
                for (const auto& off : offsets) {
                    const auto nb = grid[i].neighbor(off);
                    neighbors_sum += u.at(nb.enumerationIndex());
                }
                residual.at(i) = f.at(i) + (1.0 / h2) * (-6.0 * u.at(i) + neighbors_sum);
            }
            sum_r2 += residual.at(i) * residual.at(i);
        }

        l2_norm = std::sqrt(sum_r2 / static_cast<double>(nCells));
        iter++;
    }

    // Use \n instead of std::endl to satisfy performance-avoid-endl
    std::cout << l2_norm << "\n" << iter << "\n";
}

} // namespace

int main(int argc, char** argv) {
    try {
        // Wrap argv in a span to satisfy pointer arithmetic checks.
        const std::span<char*> args(argv, static_cast<std::size_t>(argc));

        if (argc < 4) {
            std::cerr << "Usage: " << args[0] << " <level> <maxIter> <epsilon> [outFile]\n";
            return 1;
        }

        const std::size_t level = std::stoul(args[1]);
        const std::size_t maxIter = std::stoul(args[2]);
        const double epsilon = std::stod(args[3]);

        auto octree = CellOctree::createUniformGrid(level);
        const double h = 1.0 / std::pow(2.0, static_cast<double>(level));
        const double h2 = h * h;

        const std::vector<Vec<std::ptrdiff_t, 3>> offsets = {
            {-1, 0, 0}, {1, 0, 0}, {0, -1, 0}, {0, 1, 0}, {0, 0, -1}, {0, 0, 1}};

        auto grid = CellGrid::create(octree).levels({level}).neighborhood(offsets).build();

        const std::size_t nCells = grid.size();
        std::vector<double> u(nCells, 0.0);
        std::vector<double> f(nCells, 0.0);
        std::vector<bool> isBoundary(nCells, false);

        for (std::size_t i = 0; i < nCells; ++i) {
            const auto center = grid[i].center();
            f.at(i) = std::sin(M_PI * center[0]) * std::sin(M_PI * center[1]) * std::sin(M_PI * center[2]);

            for (const auto& off : offsets) {
                if (!grid[i].neighbor(off).isValid()) {
                    isBoundary.at(i) = true;
                    break;
                }
            }
        }

        solvePoisson(grid, u, f, isBoundary, offsets, h2, maxIter, epsilon);

        if (argc > 4) {
            const std::string outFile = args[4];
            io::vtk::exportCellGrid(grid, outFile).writeGridVector("u", std::span<const double>{u});
        }

    } catch (const std::exception& e) {
        std::cerr << "Standard exception caught in main: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "Unknown exception caught in main.\n";
        return 1;
    }

    return 0;
}