#include <iostream>
#include <filesystem>
#include <chrono>
#include <fstream>
#include <cmath>
#include <iomanip>
#include <numbers>
#include <vector>
#include <string_view>
#include <span>

#include "oktal/octree/CellOctree.hpp"
#include "oktal/octree/CellGrid.hpp"
#include "oktal/lbm/D3Q19.hpp"
#include "oktal/lbm/LbmKernels.hpp"
#include "oktal/lbm/TaylorGreen.hpp"
#include "oktal/io/VtkExport.hpp"

using namespace oktal;

int main(int argc, char** argv) {
    // Wrap argv to avoid pointer arithmetic warnings
    const std::span<char*> args(argv, static_cast<size_t>(argc));

    try {
        if (argc != 3) {
            std::cerr << "Usage: " << args[0] << " <refinement-level> <output-directory>\n";
            return 1;
        }

        const size_t level = std::stoul(args[1]);
        if (level < 5) {
            std::cerr << "Error: Refinement level must be at least five.\n";
            return 1;
        }

        const std::filesystem::path outDir = args[2];
        std::filesystem::create_directories(outDir);

        const lbm::TaylorGreen tgv(level);
        auto octree = CellOctree::createUniformGrid(tgv.geometry(), level);

        const auto grid = CellGrid::create(octree)
                              .neighborhood(lbm::D3Q19::CS)
                              .periodicityMapper(Torus({true, true, true}))
                              .build();

        GridVector<double, 19> f(grid);
        GridVector<double, 19> fTmp(grid);
        GridVector<double, 1> rho(grid);
        GridVector<double, 3> u(grid);
        GridVector<double, 3> uErr(grid);

        for (const auto cell : grid) {
            const auto pos = cell.center();
            const std::size_t idx = cell.enumerationIndex();
            rho[idx] = tgv.rho(pos, 0.0);
            const auto velocity = tgv.u(pos, 0.0);
            u[idx, 0] = velocity[0];
            u[idx, 1] = velocity[1];
            u[idx, 2] = velocity[2];
        }

        const lbm::InitializePdfs initKernel;
        for (const auto cell : grid) {
            initKernel(cell, f.view(), rho.const_view(), u.const_view());
        }

        const lbm::ComputeMacroscopicQuantities macroscopic;
        const lbm::Collide collide(tgv.omega());
        const lbm::Stream stream;

        const size_t numSteps = tgv.numberOfTimesteps();
        const auto startTime = std::chrono::steady_clock::now();

        for (size_t step = 0; step <= numSteps; ++step) {
            const double tPhys = static_cast<double>(step) * tgv.dt();

            if (step % 50 == 0 || step == numSteps) {
                for (const auto cell : grid) {
                    const auto pos = cell.center();
                    const auto analyticU = tgv.u(pos, tPhys);
                    const std::size_t idx = cell.enumerationIndex();
                    uErr[idx, 0] = u[idx, 0] - analyticU[0];
                    uErr[idx, 1] = u[idx, 1] - analyticU[1];
                    uErr[idx, 2] = u[idx, 2] - analyticU[2];
                }

                auto writer = io::vtk::exportCellGrid(grid, outDir / ("step_" + std::to_string(step)));
                writer.writeGridVector<double, 1>("rho", rho.const_view());
                writer.writeGridVector<double, 3>("u", u.const_view());
                writer.writeGridVector<double, 3>("u_err", uErr.const_view());

                const auto now = std::chrono::steady_clock::now();
                const std::chrono::duration<double> elapsed = now - startTime;
                std::cout << "Step: " << step << "/" << numSteps << " (Elapsed: " << elapsed.count() << "s)\n";
            }

            if (step < numSteps) {
                for (const auto cell : grid) {
                    macroscopic(cell, f.const_view(), rho.view(), u.view());
                }
                for (const auto cell : grid) {
                    collide(cell, f.view(), rho.const_view(), u.const_view());
                }
                for (const auto cell : grid) {
                    stream(cell, fTmp.view(), f.const_view());
                }
                std::swap(f, fTmp);
            }
        }

        double sumErrSqX = 0.0;
        double sumErrSqY = 0.0;
        double sumAnalyticSqX = 0.0;
        double sumAnalyticSqY = 0.0;
        const double tFinal = static_cast<double>(numSteps) * tgv.dt();

        for (const auto cell : grid) {
            const auto analyticU = tgv.u(cell.center(), tFinal);
            const std::size_t idx = cell.enumerationIndex();

            sumErrSqX += std::pow(u[idx, 0] - analyticU[0], 2);
            sumErrSqY += std::pow(u[idx, 1] - analyticU[1], 2);
            sumAnalyticSqX += std::pow(analyticU[0], 2);
            sumAnalyticSqY += std::pow(analyticU[1], 2);
        }

        std::ofstream errFile(outDir / "errors.txt");
        errFile << std::sqrt(sumErrSqX / sumAnalyticSqX) << " "
                << std::sqrt(sumErrSqY / sumAnalyticSqY) << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Unhandled exception: " << e.what() << "\n";
        return 1;
    } catch (...) {
        std::cerr << "Unknown exception occurred.\n";
        return 1;
    }

    return 0;
}
