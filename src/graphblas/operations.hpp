#ifndef GB_OPERATIONS_HPP
#define GB_OPERATIONS_HPP

#pragma once

#include <cstddef>
#include <vector>

#include <graphblas/algebra.hpp>
#include <graphblas/TransposeView.hpp>
#include <graphblas/Matrix.hpp>
#define GB_INCLUDE_BACKEND_TRANSPOSE_VIEW 1

#define GB_INCLUDE_BACKEND_OPERATIONS 1
#include <backend_include.hpp>

namespace GraphBLAS
{
	// 4.3.1: Matrix-matrix multiply
	template<typename CMatrixT,
		typename MaskT,
		typename AccumT,
		typename SemiringT,
		typename AMatrixT,
		typename BMatrixT>

			inline Info mxm(CMatrixT         &C,
					MaskT      const &Mask,
					AccumT            accum,
					SemiringT         op,
					AMatrixT   const &A,
					BMatrixT   const &B,
					bool              replace_flag = false)
			{
				CHECK_STATUS(
				backend::mxm(C.m_mat, Mask.m_mat, accum, op, A.m_mat, B.m_mat,
						replace_flag)
						);
return SUCCESS;
			}

	// 4.3.9.3: reduce - matrix-scalar variant
	/// @note We aren't supporting transpose of matrix here. The spec does not
	/// require support.
	template<typename ValueT,
		typename AccumT,
		typename MonoidT, // monoid only
		typename AScalarT,
		typename ...ATagsT>
			inline Info reduce(
					ValueT                            &val,
					AccumT                             accum,
					MonoidT                            op,
					Matrix<AScalarT, ATagsT...> const &A)
			{
 CHECK_STATUS(
				backend::reduce_matrix_to_scalar(val, accum, op, A.m_mat)
				);
return SUCCESS;
			}

    //************************************************************************
    // Transpose
    //************************************************************************

    // 4.3.10: transpose
    template<typename CMatrixT,
             typename MaskT,
             typename AccumT,
             typename AMatrixT>
    inline Info transpose(CMatrixT       &C,
                          MaskT    const &Mask,
                          AccumT          accum,
                          AMatrixT const &A,
                          bool            replace_flag = false)
    {


        CHECK_STATUS(check_nrows_nrows(C, Mask, "transpose: C.nrows != Mask.nrows"));
        CHECK_STATUS(check_ncols_ncols(C, Mask, "transpose: C.ncols != Mask.ncols"));
        CHECK_STATUS(check_ncols_nrows(C, A, "transpose: C.ncols != A.nrows"));
        CHECK_STATUS(check_ncols_nrows(A, C, "transpose: A.ncols != C.nrows"));

        CHECK_STATUS(
            backend::transpose(C.m_mat, Mask.m_mat, accum, A.m_mat, replace_flag)
            );


        return SUCCESS;
    }

    //************************************************************************
    // Views
    //************************************************************************

    /**
     * @brief  "Flip" the rows and columns of a matrix
     * @param[in]  A  The matrix to transpose
     *
     */
    template<typename MatrixT>
    inline TransposeView<MatrixT> transpose(MatrixT const &A)
    {
        return TransposeView<MatrixT>(backend::transpose(A.m_mat));
    }



} // GraphBLAS


#endif
