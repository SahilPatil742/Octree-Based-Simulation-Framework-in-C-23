#pragma once

#include <vector>
#include <concepts>
#include <algorithm> 
#include <utility>   
#include <type_traits>
#include <cassert>
#include <span>      


#include <experimental/mdspan>


namespace oktal {


    template<typename T, size_t Q>
    requires (std::semiregular<std::remove_const_t<T>> && Q > 0)
    using GridVectorView = std::experimental::mdspan<
        T,
        std::conditional_t<
            Q == 1,
            std::experimental::extents<size_t, std::experimental::dynamic_extent>,
            std::experimental::extents<size_t, std::experimental::dynamic_extent, Q>
        >,
        std::experimental::layout_left
    >;


    template <typename T, size_t Q>
    requires (std::semiregular<T> && Q > 0)
    class GridVector {
    public:
        using value_type = T;
        static constexpr size_t NUM_COMPONENTS = Q;

        explicit GridVector(const auto& cells) 
            : m_numCells(cells.size())
            , m_capacity(cells.size() * Q)
        {
            if (m_capacity > 0) {
                m_data = new T[m_capacity](); 
            }
        }


        ~GridVector() {
            delete[] m_data;
        }

        GridVector(const GridVector& other) 
            : m_numCells(other.m_numCells)
            , m_capacity(other.m_capacity) 
        {
            if (m_capacity > 0) {
                m_data = new T[m_capacity];
                std::copy(other.m_data, other.m_data + m_capacity, m_data);
            }
        }


        GridVector& operator=(const GridVector& other) {
            if (this != &other) {
                if (m_capacity != other.m_capacity) {
                    delete[] m_data;
                    m_capacity = other.m_capacity;
                    m_data = (m_capacity > 0) ? new T[m_capacity] : nullptr;
                }
                m_numCells = other.m_numCells;
                if (m_capacity > 0) {
                    std::copy(other.m_data, other.m_data + m_capacity, m_data);
                }
            }
            return *this;
        }


        GridVector(GridVector&& other) noexcept 
            : m_data(std::exchange(other.m_data, nullptr))
            , m_numCells(std::exchange(other.m_numCells, 0))
            , m_capacity(std::exchange(other.m_capacity, 0)) 
        {}


        GridVector& operator=(GridVector&& other) noexcept {
            if (this != &other) {
                delete[] m_data;
                m_data = std::exchange(other.m_data, nullptr);
                m_numCells = std::exchange(other.m_numCells, 0);
                m_capacity = std::exchange(other.m_capacity, 0);
            }
            return *this;
        }



        size_t allocSize() const {
            return m_capacity;
        }

        T* data() { return m_data; }
        const T* data() const { return m_data; }


        GridVectorView<T, Q> view() {
            if constexpr (Q == 1) {
                return GridVectorView<T, Q>(m_data, m_numCells);
            } else {
                return GridVectorView<T, Q>(m_data, m_numCells, Q);
            }
        }

        GridVectorView<const T, Q> view() const {
            if constexpr (Q == 1) {
                return GridVectorView<const T, Q>(m_data, m_numCells);
            } else {
                return GridVectorView<const T, Q>(m_data, m_numCells, Q);
            }
        }

        GridVectorView<const T, Q> const_view() const {
            return view();
        }


        
        operator GridVectorView<T, Q>() {
            return view();
        }

        operator GridVectorView<const T, Q>() const {
            return view();
        }



        // Scalar case (Q == 1)
        decltype(auto) operator[](size_t cellIdx) requires (Q == 1) {
            return view()[cellIdx];
        }

        decltype(auto) operator[](size_t cellIdx) const requires (Q == 1) {
            return view()[cellIdx];
        }

        // Vector case (Q > 1)
        decltype(auto) operator[](size_t cellIdx, size_t compIdx) requires (Q > 1) {
            return view()[cellIdx, compIdx];
        }

        decltype(auto) operator[](size_t cellIdx, size_t compIdx) const requires (Q > 1) {
            return view()[cellIdx, compIdx];
        }

    private:
        T* m_data = nullptr;
        size_t m_numCells = 0;
        size_t m_capacity = 0;
    };


    namespace data {
        using ::oktal::GridVector;
        using ::oktal::GridVectorView;
    }

} // namespace oktal