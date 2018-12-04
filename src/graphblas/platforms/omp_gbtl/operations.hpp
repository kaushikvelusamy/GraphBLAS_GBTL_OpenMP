#ifndef GB_SEQUENTIAL_OPERATIONS_HPP

#define GB_SEQUENTIAL_OPERATIONS_HPP
#pragma once
#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <graphblas/algebra.hpp>
// #include <graphblas/platforms/sequential/TransposeView.hpp>

#include <graphblas/platforms/omp_gbtl/sparse_mxm.hpp>
#include <graphblas/platforms/omp_gbtl/sparse_reduce.hpp>
#include <graphblas/platforms/omp_gbtl/sparse_transpose.hpp>


namespace GraphBLAS
{
    namespace backend
    {

        template<typename MatrixT>
        inline TransposeView<MatrixT> transpose(MatrixT const &A)
        {
            return TransposeView<MatrixT>(A);
        }

    } // backend
} // GraphBLAS

#endif // GB_SEQUENTIAL_OPERATIONS_HPP
