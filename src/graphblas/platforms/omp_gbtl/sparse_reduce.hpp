#ifndef GB_SEQUENTIAL_SPARSE_REDUCE_HPP
#define GB_SEQUENTIAL_SPARSE_REDUCE_HPP

#pragma once

#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <iostream>
#include <graphblas/algebra.hpp>

#include "sparse_helpers.hpp"


namespace GraphBLAS
{
    namespace backend
    {
        /// Implementation of 4.3.9.3 reduce: Matrix to scalar variant
        template<typename ValueT,
                 typename AccumT,
                 typename MonoidT, // monoid only
                 typename AMatrixT>
        inline Info reduce_matrix_to_scalar(ValueT         &val,
                                            AccumT          accum,
                                            MonoidT         op,
                                            AMatrixT const &A)
        {
            // Do the basic reduction work with the monoid
            typedef typename MonoidT::result_type D3ScalarType;
            typedef typename AMatrixT::ScalarType AScalarType;
            typedef std::vector<std::tuple<IndexType,AScalarType> >  ARowType;

            D3ScalarType t = op.identity();

            if (A.nvals() > 0)
            {
                for (IndexType row_idx = 0; row_idx < A.nrows(); ++row_idx)
                {
                    /// @todo Can't be a reference because A might be transpose
                    /// view.  Need to specialize on TransposeView and getCol()
                    ARowType const A_row(A.getRow(row_idx));

                    /// @todo There is something hinky with domains here.  How
                    /// does one perform the reduction in A domain but produce
                    /// partial results in D3(op)?
                    D3ScalarType tmp;
                    if (reduction(tmp, A_row, op))
                    {
                        t = op(t, tmp); // reduce each row
                    }
                }
            }

            // =================================================================
            // Accumulate into Z
            /// @todo Do we need a type generator for z: D(w) if no accum,
            /// or D3(accum). I think that D(z) := D(val) should be equivalent, but
            /// still need to work the proof.
            ValueT z;
            opt_accum_scalar(z, val, t, accum);

            // Copy Z into the final output
            val = z;
            return SUCCESS;
        }

    } // backend
} // GraphBLAS

#endif
