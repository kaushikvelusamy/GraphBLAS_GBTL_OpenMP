#pragma once

#include <cstddef>
#include <graphblas/platforms/omp_gbtl/LilSparseMatrix.hpp>


namespace GraphBLAS
{
    namespace backend
    {
        // A marker class for when we should have no mask
        // @todo: Find somewhere else to put this
        class NoMask
        {
        public:
            NoMask() {}

            friend std::ostream &operator<<(std::ostream             &os,
                                            NoMask          const    &mask)
            {
                os << "No mask";
                return os;
            }
        };


        template<typename ScalarT, typename... TagsT>
        class Matrix : public LilSparseMatrix<ScalarT>
        {
        public:
            typedef ScalarT ScalarType;

            // construct an empty matrix of fixed dimensions
            Matrix(IndexType   num_rows,
                   IndexType   num_cols)
                : LilSparseMatrix<ScalarT>(num_rows, num_cols)
            {
            }

            // copy construct
            Matrix(Matrix const &rhs)
                : LilSparseMatrix<ScalarT>(rhs)
            {
            }

            // construct a dense matrix from dense data.
            Matrix(std::vector<std::vector<ScalarT> > const &values)
                : LilSparseMatrix<ScalarT>(values)
            {
            }

            // construct a sparse matrix from dense data and a zero val.
            Matrix(std::vector<std::vector<ScalarT> > const &values,
                   ScalarT                                   zero)
                : LilSparseMatrix<ScalarT>(values, zero)
            {
            }

            ~Matrix() {}  // virtual?

       };
    }
}

