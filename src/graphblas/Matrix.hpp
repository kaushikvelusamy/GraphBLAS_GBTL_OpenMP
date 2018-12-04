#pragma once

#include <cstddef>
#include <type_traits>
#include <graphblas/detail/param_unpack.hpp>

#define GB_INCLUDE_BACKEND_MATRIX 1
#include <backend_include.hpp>

namespace GraphBLAS
{
	// We need to declare class so we can include later.
	template<typename ScalarT, typename... TagsT>
		class Vector;

	template<typename ScalarT, typename... TagsT>
		class Matrix;

	template<typename MatrixT>
		class TransposeView;


	template<typename ScalarT, typename... TagsT>
		class Matrix
		{
			public:
				typedef matrix_tag  tag_type;

				typedef ScalarT     ScalarType;
				typedef typename detail::matrix_generator::result<
					ScalarT,
					detail::SparsenessCategoryTag,
					detail::DirectednessCategoryTag,
					TagsT... ,
					detail::NullTag,
					detail::NullTag >::type BackendType;

				/**
				 * @brief Construct an empty matrix with the specified shape.
				 *
				 * @note The backend should be able to decide when to ignore any of the
				 *       tags and/or arguments.
				 *
				 * @param[in] num_rows  Number of rows in the matrix
				 * @param[in] num_cols  Number of columns in the matrix
				 * @param[in] zero      The "zero" value, additive identity, and
				 *                      the structural zero.
				 */
				Matrix(IndexType num_rows, IndexType num_cols)
					: m_mat(num_rows, num_cols)
				{
				}

				/**
				 * @brief Copy constructor.
				 *
				 * @param[in] rhs   The matrix to copy.
				 */
				Matrix(Matrix<ScalarT, TagsT...> const &rhs)
					: m_mat(rhs.m_mat)
				{
				}

				/**
				 * @brief Construct a dense matrix from dense data
				 *
				 * @param[in] values The dense matrix from which to construct a
				 *                   sparse matrix from.
				 *
				 * @todo Should we really support this interface?
				 */
				Matrix(std::vector<std::vector<ScalarT> > const &values)
					: m_mat(values)
					{
					}

				/**
				 * @brief Construct a sparse matrix from dense data and a sentinel zero value.
				 *
				 * @param[in] values The dense matrix from which to construct a
				 *                   sparse matrix from.
				 * @param[in] zero   The "zero" value used to determine implied
				 *                   zeroes (no stored value) in the sparse structure
				 *
				 * @todo Should we really support this interface?
				 */
				Matrix(std::vector<std::vector<ScalarT> > const &values, ScalarT zero)
					: m_mat(values, zero)
					{
					}

				~Matrix() { }

				/// @todo Should assignment work only if dimensions are same?
				Matrix<ScalarT, TagsT...> &
					operator=(Matrix<ScalarT, TagsT...> const &rhs)
					{
						if (this != &rhs)
						{
							// backend currently doing dimension check.
							m_mat = rhs.m_mat;
						}
						return *this;
					}


				/// @todo need to change to mix and match internal types
				bool operator==(Matrix<ScalarT, TagsT...> const &rhs) const
				{
					return (m_mat == rhs.m_mat);
				}

				bool operator!=(Matrix<ScalarT, TagsT...> const &rhs) const
				{
					//return !(m_mat == rhs.m_mat);
					return !(*this == rhs);
				}

				/**
				 * Populate the matrix with stored values (using iterators).
				 *
				 * @param[in]  i_it      Row index iterator
				 * @param[in]  j_it      Column index iterator
				 * @param[in]  v_it      Value (scalar) iterator
				 * @param[in]  num_vals  Number of elements to store
				 * @param[in]  dup       Binary function to call when value is being stored
				 *                       in a location that already has a stored value.
				 *                       stored_val = dup(stored_val, *v_it)
				 *
				 * @todo The C spec says it is an error to call build on a non-empty
				 *       matrix.  Unclear if the C++ should.
				 */
				template<typename RAIteratorI,
					typename RAIteratorJ,
					typename RAIteratorV,
					typename BinaryOpT = GraphBLAS::Second<ScalarType> >
						void build(RAIteratorI  i_it,
								RAIteratorJ  j_it,
								RAIteratorV  v_it,
								IndexType    num_vals,
								BinaryOpT    dup = BinaryOpT())
						{

							m_mat.build(i_it, j_it, v_it, num_vals, dup);
						}

				/**
				 * Populate the matrix with stored values (using iterators).
				 *
				 * @param[in]  row_indices  Array of row indices
				 * @param[in]  col_indices  Array of column indices
				 * @param[in]  values       Array of values
				 * @param[in]  dup          binary function to call when value is being stored
				 *                          in a location that already has a stored value.
				 *                          stored_val = dup(stored_val, *v_it)
				 *
				 * @todo The C spec says it is an error to call build on a non-empty
				 *       matrix.  Unclear if the C++ should.
				 */
				template<typename ValueT,
					typename BinaryOpT = GraphBLAS::Second<ScalarType> >
						inline Info build(IndexArrayType       const &row_indices,
								IndexArrayType       const &col_indices,
								std::vector<ValueT>  const &values,
								BinaryOpT                   dup = BinaryOpT())
						{

							if ((row_indices.size() != col_indices.size()) ||
									(row_indices.size() != values.size()))
							{
								std::cerr << "throw DimensionException(\"Matrix::build\")\n";
							}

							return m_mat.build(row_indices.begin(), col_indices.begin(),
									values.begin(), values.size(), dup);
						}

				void clear() { m_mat.clear(); }

				IndexType nrows() const  { return m_mat.nrows(); }
				IndexType ncols() const  { return m_mat.ncols(); }
				IndexType nvals() const  { return m_mat.nvals(); }

				bool hasElement(IndexType row, IndexType col) const
				{
					return m_mat.hasElement(row, col);
				}

				/// @todo I don't think this is a valid interface for sparse
				void setElement(IndexType row, IndexType col, ScalarT const &val)
				{
					m_mat.setElement(row, col, val);
				}

				/// @throw NoValueException if there is no value stored at (row,col)
				ScalarT extractElement(IndexType row, IndexType col) const
				{
					return m_mat.extractElement(row, col);
				}

				template<typename RAIteratorIT,
					typename RAIteratorJT,
					typename RAIteratorVT>
						inline void extractTuples(RAIteratorIT        row_it,
								RAIteratorJT        col_it,
								RAIteratorVT        values) const
						{
							m_mat.extractTuples(row_it, col_it, values);
						}

				template <typename RowSequenceT,
					 typename ColSequenceT>
						 inline void extractTuples(RowSequenceT            &row_indices,
								 ColSequenceT            &col_indices,
								 std::vector<ScalarT>    &values) const
						 {
							 m_mat.extractTuples(row_indices.begin(),
									 col_indices.begin(),
									 values.begin());
						 }

			private:

				// 4.3.1:
				template<typename CMatrixT,
					typename MaskT,
					typename AccumT,
					typename SemiringT,
					typename AMatrixT,
					typename BMatrixT>
						friend inline Info mxm(CMatrixT         &C,
								MaskT      const &Mask,
								AccumT            accum,
								SemiringT         op,
								AMatrixT   const &A,
								BMatrixT   const &B,
								bool              replace_flag);

				// 4.3.9.1
				template<typename WVectorT,
					typename MaskT,
					typename AccumT,
					typename BinaryOpT,  // monoid or binary op only
					typename AMatrixT>
						friend inline Info reduce(WVectorT        &u,
								MaskT     const &mask,
								AccumT           accum,
								BinaryOpT        op,
								AMatrixT  const &A,
								bool             replace_flag);

				// 4.3.9.3
				template<typename ValueT,
					typename AccumT,
					typename MonoidT, // monoid only
					typename AScalarT,
					typename... ATagsT>
						friend inline Info reduce(
								ValueT                                       &dst,
								AccumT                                        accum,
								MonoidT                                       op,
								GraphBLAS::Matrix<AScalarT, ATagsT...> const &A);

				//--------------------------------------------------------------------

				// 4.3.10
				template<typename CMatrixT,
					typename MaskT,
					typename AccumT,
					typename AMatrixT>
						friend inline Info transpose(CMatrixT       &C,
								MaskT    const &Mask,
								AccumT          accum,
								AMatrixT const &A,
								bool            replace_flag);

				//--------------------------------------------------------------------

				template<typename MatrixT>
					friend inline GraphBLAS::TransposeView<MatrixT> transpose(MatrixT const &A);




			private:
				BackendType m_mat;
		};

	//**************************************************************************
	// GrB_NULL mask: should be GrB_FULL
	class NoMask
	{
		public:
			typedef bool ScalarType; // not necessary?
			typedef backend::NoMask BackendType; // not necessary?

			backend::NoMask m_mat;  // can be const?
			backend::NoMask m_vec;
	};



} // end namespace GraphBLAS


