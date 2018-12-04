#ifndef ALGORITHMS_TRIANGLE_COUNT_HPP
#define ALGORITHMS_TRIANGLE_COUNT_HPP

#include <iostream>
#include <chrono>

#include <graphblas/graphblas.hpp>

namespace algorithms

{
   //************************************************************************
    template<typename MatrixT>
    typename MatrixT::ScalarType triangle_count_masked(MatrixT const &L)
    {
        using T = typename MatrixT::ScalarType;
        GraphBLAS::IndexType rows(L.nrows());
        GraphBLAS::IndexType cols(L.ncols());

        std::cout << "DEBUG: allocating temporary matrix." << std::endl;
        MatrixT B(rows, cols);

        std::cout << "DEBUG: before mxm function"<< std::endl;
        GraphBLAS::mxm(B,
                       L, GraphBLAS::NoAccumulate(),
                       GraphBLAS::ArithmeticSemiring<T>(),
                       L, GraphBLAS::transpose(L));

        std::cout << "DEBUG: before reduce function"<< std::endl;
        T sum = 0;
        GraphBLAS::reduce(sum,
                          GraphBLAS::NoAccumulate(),
                          GraphBLAS::PlusMonoid<T>(),
                          B);
        return sum;
    }
/*
   ************************************************************************
    template<typename LMatrixT, typename MatrixT>
    typename MatrixT::ScalarType triangle_count_newGBTL(LMatrixT const &L,
                                                        MatrixT  const &U)
    {
        auto start = std::chrono::steady_clock::now();

        using T = typename MatrixT::ScalarType;

	GraphBLAS::IndexType rows(L.nrows());
        GraphBLAS::IndexType cols(L.ncols());

        MatrixT B(rows, cols);
        GraphBLAS::mxm(B, GraphBLAS::NoMask(), GraphBLAS::NoAccumulate(),
                       GraphBLAS::ArithmeticSemiring<T>(),
                       L,
                       U);  /// @todo can't use transpose(L) here as LMatrix may
                            /// already be a TransposeView (nesting not supported)

        auto finish = std::chrono::steady_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::milliseconds>
            (finish - start);
        start = finish;
        //std::cout << "mxm elapsed time: " << duration.count() << " msec." << std::endl;


        T sum = 0;
        MatrixT C(rows, cols);
        GraphBLAS::eWiseMult(C, GraphBLAS::NoMask(), GraphBLAS::NoAccumulate(),
                             GraphBLAS::Times<T>(),
                             L, B, true);

        GraphBLAS::reduce(sum, GraphBLAS::NoAccumulate(),
                          GraphBLAS::PlusMonoid<T>(), C);
        finish = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>
            (finish - start);
        start = finish;
        //std::cout << "count1 elapsed time: " << duration.count() << " msec." << std::endl;

        // for undirected graph you can stop here and return 'sum'

        GraphBLAS::eWiseMult(C, GraphBLAS::NoMask(), GraphBLAS::NoAccumulate(),
                             GraphBLAS::Times<T>(),
                             U, B, true);

        GraphBLAS::reduce(sum, GraphBLAS::Plus<T>(),
                          GraphBLAS::PlusMonoid<T>(), C);
        finish = std::chrono::steady_clock::now();
        duration = std::chrono::duration_cast<std::chrono::milliseconds>
            (finish - start);
        //std::cout << "count2 elapsed time: " << duration.count() << " msec." << std::endl;

        return sum / static_cast<T>(2);
    }
  */
} // algorithms

#endif // ALGORITHMS_TRIANGLE_COUNT_HPP
