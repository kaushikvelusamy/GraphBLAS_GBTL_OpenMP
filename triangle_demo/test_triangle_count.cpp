#include <iostream>
#include <graphblas/graphblas.hpp>
#include <algorithms/triangle_count.hpp>

using namespace GraphBLAS;
using namespace algorithms;

#define BOOST_TEST_MAIN
#define BOOST_TEST_MODULE triangle_count_test_suite

#include <boost/test/included/unit_test.hpp>

BOOST_AUTO_TEST_SUITE(BOOST_TEST_MODULE)

//****************************************************************************
BOOST_AUTO_TEST_CASE(test_triangle_count_masked)
{
    //Matrix<double, DirectedMatrixTag> testtriangle(
    //                       {{0,1,1,1,0},
    //                        {1,0,1,0,1},
    //                        {1,1,0,1,1},
    //                        {1,0,1,0,1},
    //                        {0,1,1,1,0}});

    std::vector<double> ar={0, 0, 0, 1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 4, 4, 4};
    std::vector<double> ac={1, 2, 3, 0, 2, 4, 0, 1, 3, 4, 0, 2, 4, 1, 2, 3};
    std::vector<double> av={1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    Matrix<double, DirectedMatrixTag> testtriangle(5,5);
    testtriangle.build(ar.begin(), ac.begin(), av.begin(), av.size());

    Matrix<double, DirectedMatrixTag> L(5,5), U(5,5);
//    GraphBLAS::split(testtriangle, L, U);

    IndexType result = triangle_count_masked(L);
    BOOST_CHECK_EQUAL(result, 4);
}

BOOST_AUTO_TEST_SUITE_END()
