#pragma once

#include "oktal/octree/CellGrid.hpp"
#include "oktal/lbm/D3Q19.hpp"

namespace oktal::lbm {

/**
 * @brief Helper to compute the discrete equilibrium distribution.
 */
inline double computeEquilibrium(std::size_t i, double rho, const std::array<double, 3>& u) {
    const double ci_u = static_cast<double>(D3Q19::CS[i][0]) * u[0] + 
                        static_cast<double>(D3Q19::CS[i][1]) * u[1] + 
                        static_cast<double>(D3Q19::CS[i][2]) * u[2];
    
    const double u_u = u[0] * u[0] + u[1] * u[1] + u[2] * u[2];
    
    return D3Q19::W[i] * rho * (1.0 + 3.0 * ci_u + 4.5 * ci_u * ci_u - 1.5 * u_u);
}

class InitializePdfs {
public:
    void operator()(CellGrid::CellView cell,
                    D3Q19LatticeView pdfs,
                    GridVectorView<const double, 1> rho,
                    GridVectorView<const double, 3> u) const {
        const std::size_t cellIdx{cell};
        
        // rho has rank 1: access with one index
        const double r = rho[cellIdx];
        // u has rank 2 (cell index, dimension): access with comma
        const std::array<double, 3> vel = {u[cellIdx, 0], u[cellIdx, 1], u[cellIdx, 2]};

        for (std::size_t i = 0; i < D3Q19::Q; ++i) {
            pdfs[cellIdx, i] = computeEquilibrium(i, r, vel);
        }
    }
};

class ComputeMacroscopicQuantities {
public:
    void operator()(CellGrid::CellView cell,
                    D3Q19LatticeConstView pdfs,
                    GridVectorView<double, 1> rho,
                    GridVectorView<double, 3> u) const {
        const std::size_t cellIdx{cell};
        
        double r = 0.0;
        for (std::size_t i = 0; i < D3Q19::Q; ++i) {
            r += pdfs[cellIdx, i];
        }
        rho[cellIdx] = r;

        const double invRho = 1.0 / r;
        for (std::size_t d = 0; d < 3; ++d) {
            double momentum_d = 0.0;
            for (std::size_t i = 0; i < D3Q19::Q; ++i) {
                momentum_d += static_cast<double>(D3Q19::CS[i][d]) * pdfs[cellIdx, i];
            }
            u[cellIdx, d] = invRho * momentum_d;
        }
    }
};

class Collide {
public:
    explicit Collide(double omega) : m_omega(omega) {}

    void operator()(CellGrid::CellView cell,
                    D3Q19LatticeView pdfs,
                    GridVectorView<const double, 1> rho,
                    GridVectorView<const double, 3> u) const {
        const std::size_t cellIdx{cell};
        const double r = rho[cellIdx];
        const std::array<double, 3> vel = {u[cellIdx, 0], u[cellIdx, 1], u[cellIdx, 2]};

        for (std::size_t i = 0; i < D3Q19::Q; ++i) {
            const double feq = computeEquilibrium(i, r, vel);
            pdfs[cellIdx, i] += m_omega * (feq - pdfs[cellIdx, i]);
        }
    }

private:
    double m_omega;
};

class Stream {
public:
    void operator()(CellGrid::CellView cell,
                    D3Q19LatticeView pdfsDst,
                    D3Q19LatticeConstView pdfsSrc) const {
        const std::size_t cellIdx{cell};

        // Handle stationary population (i=0) explicitly
        // It stays in the same cell: x + (0,0,0) = x
        pdfsDst[cellIdx, 0] = pdfsSrc[cellIdx, 0];

        // Handle moving populations (i=1 to 18)
        for (std::size_t i = 1; i < D3Q19::Q; ++i) {
            auto neighbor = cell.neighbor(D3Q19::CS[i]);
            
            if (neighbor) {
                const std::size_t neighborIdx{neighbor};
                // Equation 22: f_i(x + c_i, t+1) = f_i(x, t)
                pdfsDst[neighborIdx, i] = pdfsSrc[cellIdx, i];
            }
        }
    }
};

} // namespace oktal::lbm
