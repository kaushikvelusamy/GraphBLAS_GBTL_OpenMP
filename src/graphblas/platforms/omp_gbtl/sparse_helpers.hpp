#ifndef GB_SEQUENTIAL_SPARSE_HELPERS_HPP
#define GB_SEQUENTIAL_SPARSE_HELPERS_HPP

#pragma once

#include <functional>
#include <utility>
#include <vector>
#include <iterator>
#include <iostream>
#include <string>
#include <graphblas/algebra.hpp>

namespace GraphBLAS
{
    namespace backend
    {
        template <typename DstMatrixT,
                  typename SrcMatrixT>
        void sparse_copy(DstMatrixT &dstMatrix,
                         SrcMatrixT const &srcMatrix)
        {
            typedef typename SrcMatrixT::ScalarType SrcScalarType;
            typedef typename DstMatrixT::ScalarType DstScalarType;

            typedef std::vector<std::tuple<IndexType, SrcScalarType> > SrcRowType;
            typedef std::vector<std::tuple<IndexType, DstScalarType> > DstRowType;

            // Copying removes the contents of the other matrix so clear it first.
            std::cout << "DEBUG: sparse_copy clear dst" << std::endl;
            dstMatrix.clear();
            dstMatrix.printInfo(std::cout, "dstMatrix");

            srcMatrix.printInfo(std::cout, "srcMatrix");

            IndexType nrows(dstMatrix.nrows());
            for (IndexType row_idx = 0; row_idx < nrows; ++row_idx)
            {
                std::cout << "DEBUG: sparse_copy copy row " << row_idx << std::endl;
                SrcRowType srcRow = srcMatrix.getRow(row_idx);
                DstRowType dstRow;

                // We need to construct a new row with the appropriate cast!
                for (auto it = srcRow.begin(); it != srcRow.end(); ++it)
                {
                    IndexType idx;
                    SrcScalarType srcVal;
                    std::tie(idx, srcVal) = *it;
                    dstRow.push_back(std::make_tuple(idx, static_cast<DstScalarType>(srcVal)));
                }

                if (!dstRow.empty())
                    dstMatrix.setRow(row_idx, dstRow);
            }
        }

        // @todo: Make a sparse copy where they are the same type for efficiency

        //**********************************************************************
        /// Increments the provided iterate while the value is less
        /// than the provided index
        template <typename I>
        void increment_until_true(
            // std::vector<std::tuple<GraphBLAS::IndexType,D> >::const_iterator    &iter,
            // std::vector<std::tuple<GraphBLAS::IndexType,D> >::const_iterator    &iter_end,
            I          &iter,
            I const    &iter_end)
        {
            using value_type = typename std::iterator_traits<I>::value_type;
            using data_type = typename std::tuple_element<1, value_type>::type;

            GraphBLAS::IndexType tmp_idx;
            data_type tmp_val;

            if (iter == iter_end)
                return;

            std::tie(tmp_idx, tmp_val) = *iter;
            while (!tmp_val && iter != iter_end)
            {
                //std::cout << "Iter not true. index: " + std::to_string(tmp_idx) +
                //    " incrementing."  << std::endl;;
                ++iter;
                if (iter == iter_end)
                    return;

                // Reload
                std::tie(tmp_idx, tmp_val) = *iter;
            }
        }

        //**********************************************************************

        /// Increments the provided iterate while the value is less than the provided index
        template <typename I>
        void increment_while_below(
                // std::vector<std::tuple<GraphBLAS::IndexType,D> >::const_iterator    &iter,
                // std::vector<std::tuple<GraphBLAS::IndexType,D> >::const_iterator    &iter_end,
                I         &iter,
                I const   &iter_end,
                GraphBLAS::IndexType                                                 idx)
        {
            using value_type = typename std::iterator_traits<I>::value_type;
            using data_type = typename std::tuple_element<1, value_type>::type;

            GraphBLAS::IndexType tmp_idx;
            data_type tmp_val;

            if (iter == iter_end)
                return;

            std::tie(tmp_idx, tmp_val) = *iter;
            while (tmp_idx < idx && iter != iter_end)
            {
                //std::cerr << "Iter at: " + std::to_string(tmp_idx) +  ", below index: " +
                //    std::to_string(idx) + " incrementing."  << std::endl;;
                ++iter;
                if (iter == iter_end)
                    return;

                // Reload
                std::tie(tmp_idx, tmp_val) = *iter;
            }
        }

        //**********************************************************************
        template <typename I, typename V>
        void increment_and_add_while_below(
                I                              &iter,
                I                       const  &iter_end,
                GraphBLAS::IndexType            idx,
                V                               &vec)
        {
            // @todo: Tighten up the template types above to bind the vector
            // type to the iter type explicitly
            using value_type = typename std::iterator_traits<I>::value_type;
            using data_type = typename std::tuple_element<1, value_type>::type;

            GraphBLAS::IndexType tmp_idx;
            data_type tmp_val;

            if (iter == iter_end)
                return;

            std::tie(tmp_idx, tmp_val) = *iter;
            while (iter != iter_end && tmp_idx < idx )
            {
                vec.push_back(*iter);

                ++iter;
                if (iter == iter_end)
                    return;

                // Reload
                std::tie(tmp_idx, tmp_val) = *iter;
            }
        }

        //************************************************************************
        /// A dot product of two sparse vectors (vectors<tuple(index,value)>)
        template <typename D1, typename D2, typename D3, typename SemiringT>
        bool dot(D3                                                      &ans,
                 std::vector<std::tuple<GraphBLAS::IndexType,D1> > const &vec1,
                 std::vector<std::tuple<GraphBLAS::IndexType,D2> > const &vec2,
                 SemiringT                                                op)
        {
            bool value_set(false);
            ans = op.zero();

            if (vec2.empty() || vec1.empty())
            {
                return value_set;
            }

            auto v1_it = vec1.begin();
            auto v2_it = vec2.begin();

            // pull first value out of the row
            D1 a_val;
            D2 u_val;
            GraphBLAS::IndexType a_idx, u_idx;

            // loop through both ordered sets to compute sparse dot prod
            while ((v1_it != vec1.end()) &&
                   (v2_it != vec2.end()))
            {
                std::tie(a_idx, a_val) = *v1_it;
                std::tie(u_idx, u_val) = *v2_it;

                //std::cerr << "Examine u idx,val = " << u_idx << "," << u_val
                //          << "; A col_idx,val = " << a_idx << "," << a_val << std::endl;

                if (u_idx == a_idx)
                {
                    //std::cerr << ans << " + " << a_val << " * " << u_val << " = ";
                    ans = op.add(ans, op.mult(a_val, u_val));
                    value_set = true;
                    //std::cerr << ans << std::endl;

                    ++v2_it;
                    ++v1_it;
                }
                else if (u_idx > a_idx)
                {
                    //std::cerr << "Advancing v1_it" << std::endl;
                    ++v1_it;
                }
                else
                {
                    //std::cerr << "Advancing v2_it" << std::endl;
                    ++v2_it;
                }
            }

            return value_set;
        }

        //************************************************************************
        /// A reduction of a sparse vector (vector<tuple(index,value)>) using a
        /// binary op or a monoid.
        template <typename D1, typename D3, typename BinaryOpT>
        bool reduction(
            D3                                                      &ans,
            std::vector<std::tuple<GraphBLAS::IndexType,D1> > const &vec,
            BinaryOpT                                                op)
        {
            if (vec.empty())
            {
                return false;
            }

            typedef typename BinaryOpT::result_type D3ScalarType;
            D3ScalarType tmp;

            if (vec.size() == 1)
            {
                tmp = static_cast<D3ScalarType>(std::get<1>(vec[0]));
            }
            else
            {
                /// @note Since op is associative and commutative left to right
                /// ordering is not strictly required.
                tmp = op(std::get<1>(vec[0]), std::get<1>(vec[1]));

                for (size_t idx = 2; idx < vec.size(); ++idx)
                {
                    tmp = op(tmp, std::get<1>(vec[idx]));
                }
            }

            ans = static_cast<D3>(tmp);
            return true;
        }

        //**********************************************************************
        /// Apply element-wise operation to union on sparse vectors.
        template <typename D1, typename D2, typename D3, typename BinaryOpT>
        void ewise_or(std::vector<std::tuple<GraphBLAS::IndexType,D3> >       &ans,
                      std::vector<std::tuple<GraphBLAS::IndexType,D1> > const &vec1,
                      std::vector<std::tuple<GraphBLAS::IndexType,D2> > const &vec2,
                      BinaryOpT                                                op)
        {
            ans.clear();
            auto v1_it = vec1.begin();
            auto v2_it = vec2.begin();

            D1 v1_val;
            D2 v2_val;
            GraphBLAS::IndexType v1_idx, v2_idx;

            // loop through both ordered sets to compute ewise_or
            while ((v1_it != vec1.end()) || (v2_it != vec2.end()))
            {
                if ((v1_it != vec1.end()) && (v2_it != vec2.end()))
                {
                    std::tie(v1_idx, v1_val) = *v1_it;
                    std::tie(v2_idx, v2_val) = *v2_it;

                    if (v2_idx == v1_idx)
                    {
                        //std::cerr << ans << " + " << v1_val << " * " << v2_val << " = ";
                        ans.push_back(std::make_tuple(v1_idx,
                                                      static_cast<D3>(op(v1_val, v2_val))));
                        //std::cerr << ans << std::endl;

                        ++v2_it;
                        ++v1_it;
                    }
                    else if (v2_idx > v1_idx)
                    {
                        //std::cerr << "Copying v1, Advancing v1_it" << std::endl;
                        ans.push_back(std::make_tuple(v1_idx,
                                                      static_cast<D3>(v1_val)));
                        ++v1_it;
                    }
                    else
                    {
                        //std::cerr << "Copying v2, Advancing v2_it" << std::endl;
                        ans.push_back(std::make_tuple(v2_idx,
                                                      static_cast<D3>(v2_val)));
                        ++v2_it;
                    }
                }
                else if (v1_it != vec1.end())
                {
                    std::tie(v1_idx, v1_val) = *v1_it;
                    ans.push_back(std::make_tuple(v1_idx, static_cast<D3>(v1_val)));
                    ++v1_it;
                }
                else // v2_it != vec2.end())
                {
                    std::tie(v2_idx, v2_val) = *v2_it;
                    ans.push_back(std::make_tuple(v2_idx, static_cast<D3>(v2_val)));
                    ++v2_it;
                }
            }
        }



       template < typename ZMatrixT,
                   typename CMatrixT,
                   typename TMatrixT,
                   typename BinaryOpT >
        void ewise_or_opt_accum(ZMatrixT         &Z,
                                CMatrixT const   &C,
                                TMatrixT const   &T,
                                BinaryOpT         accum)
        {
            typedef typename ZMatrixT::ScalarType ZScalarType;

            typedef std::vector<std::tuple<IndexType,ZScalarType> > ZRowType;

            ZRowType tmp_row;
            IndexType nRows(Z.nrows());
            for (IndexType row_idx = 0; row_idx < nRows; ++row_idx)
            {
                ewise_or(tmp_row, C.getRow(row_idx), T.getRow(row_idx), accum);
                Z.setRow(row_idx, tmp_row);
            }
        }

        //**********************************************************************

        // Specialized version that gets used when we don't have an accumulator
        template < typename ZMatrixT,
                   typename CMatrixT,
                   typename TMatrixT>
        void ewise_or_opt_accum(ZMatrixT                    &Z,
                                CMatrixT const              &C,
                                TMatrixT const              &T,
                                GraphBLAS::NoAccumulate )
        {
            T.printInfo(std::cout, "ewise_or_opt_accum, T (src)");
            Z.printInfo(std::cout, "ewise_or_opt_accum, Z (dst)");
            sparse_copy(Z, T);
        }

        template <typename ZScalarT,
                  typename WScalarT,
                  typename TScalarT,
                  typename BinaryOpT>
        void opt_accum_scalar(ZScalarT       &z,
                              WScalarT const &w,
                              TScalarT const &t,
                              BinaryOpT       accum)
        {
            z = static_cast<ZScalarT>(accum(w, t));
        }

        // Specialized version that gets used when we don't have an accumulator
        template <typename ZScalarT,
                  typename WScalarT,
                  typename TScalarT>
        void opt_accum_scalar(ZScalarT                &z,
                              WScalarT const          &w,
                              TScalarT const          &t,
                              GraphBLAS::NoAccumulate  accum)
        {
            z = static_cast<ZScalarT>(t);
        }

        //************************************************************************
        //**********************************************************************
        //**********************************************************************
        /**
         * This merges the values of C and Z into the result vector based on the
         * values of the mask.
         *
         * If replace:
         *
         * L(C) = \f[ {(i,j,Zij):(i,j) \in (ind(Z) \cap ind(M))} \f]
         *
         * If NOT replace:
         *
         * L(C) = {(i,j,Zij):(i,j) \in (ind(C) \cap int(\not M))} \cup
         *        {(i,j,Zij):(i,j) \in (ind(Z) \cap int(\not M))}
         *
         * @tparam CScalarT The scalar type of the C vector input AND result.
         * @tparam ZScalarT The scalar type of the Z vector input.
         * @tparam MScalarT The scalar type of the mask vector.
         *
         * @param result Result vector.  We clear this first.
         * @param c_vec The original c values that may be carried through.
         * @param z_vec The new values to insert/overlay.
         * @param mask_vec The mask which specifies which values to use.
         * @param replace If true, we should always clear the values specified
         *                by the mask regardless if they are overlayed.
         */
        template < typename CScalarT,
                   typename ZScalarT,
                   typename MScalarT>
        void apply_with_mask(std::vector<std::tuple<IndexType, CScalarT> >          &result,
                             std::vector<std::tuple<IndexType, CScalarT> > const    &c_vec,
                             std::vector<std::tuple<IndexType, ZScalarT> > const    &z_vec,
                             std::vector<std::tuple<IndexType, MScalarT> > const    &mask_vec,
                             bool                                                    replace)
        {
            auto c_it = c_vec.begin();
            auto z_it = z_vec.begin();
            auto mask_it = mask_vec.begin();

            CScalarT c_val;
            ZScalarT z_val;
            MScalarT mask_val;
            GraphBLAS::IndexType c_idx, z_idx, mask_idx;

            //std::cerr << "Executing apply_with_mask with mask and replace: " << replace << std::endl;
            //print_vec(std::cerr, "c_vec", c_vec);
            //print_vec(std::cerr, "z_vec", z_vec);
            //print_vec(std::cerr, "m_vec", mask_vec);

            result.clear();

            // Design: This approach is driven by the mask.
            while (mask_it != mask_vec.end())
            {
                // Make sure the mask is on a valid value
                increment_until_true(mask_it, mask_vec.end());

                // If we run out of mask, we are done!
                if (mask_it == mask_vec.end())
                {
                    //std::cerr << "Mask exhausted(1)." << std::endl;
                    break;
                }

                // Get the mask values
                std::tie(mask_idx, mask_val) = *mask_it;

                // If replace then we don't consider original C values.
                // If not replace, then we want to keep C values outide the mask
                // Any values in C "outside" the mask should now be applied
                // So, we catch "c" up to the mask.  This is the intersection
                // of C and !M.

                if (!replace)
                {
                    // This is the part of (ind(C) \cap int(\not M)
                    increment_and_add_while_below(c_it, c_vec.end(), mask_idx,
                                                  result);
                }

                // "Catch up" the input to the mask.
                increment_while_below(z_it, z_vec.end(), mask_idx);

                // Now, at the mask point add the value from Z if we have one.
                if (z_it != z_vec.end())
                {
                    std::tie(z_idx, z_val) = *z_it;
                    if (z_idx == mask_idx)
                    {
                        result.push_back(std::make_tuple(mask_idx, static_cast<CScalarT>(z_val)));
                        //std::cerr << "Copying v1. val: " << std::to_string(z_val) << std::endl;
                    }
                }

                // If there is a C here, skip it
                if (c_it != c_vec.end())
                {
                    std::tie(c_idx, c_val) = *c_it;
                    if (c_idx == mask_idx)
                        ++c_it;
                }

                // Now move to the next mask entry.
                ++mask_it;
            } // while mask_it != end

            // Now, we need to add the remaining C values beyond the mask
            if (!replace)
            {
                // This is the part of (ind(C) \cap int(\not M)
                while (c_it != c_vec.end())
                {
                    //std::tie(c_idx, c_val) = *c_it;
                    //std::cerr << "Catch up= " << c_idx << ":" << c_val << std::endl;
                    //result.push_back(std::make_tuple(c_idx, static_cast<CScalarT>(c_val)));
                    result.push_back(*c_it);
                    ++c_it;
                }
            }

            //print_vec(std::cerr, "result", result);

        } // apply_with_mask

        //**********************************************************************
        // Matrix version

        template < typename CMatrixT,
                   typename ZMatrixT,
                   typename MMatrixT>
        void write_with_opt_mask(CMatrixT           &C,
                                 ZMatrixT   const   &Z,
                                 MMatrixT   const   &mask,
                                 bool               replace)
        {
            typedef typename CMatrixT::ScalarType CScalarType;
            typedef typename ZMatrixT::ScalarType ZScalarType;
            typedef typename MMatrixT::ScalarType MScalarType;

            typedef std::vector<std::tuple<IndexType, CScalarType> > CRowType;
            typedef std::vector<std::tuple<IndexType, ZScalarType> > ZRowType;
            typedef std::vector<std::tuple<IndexType, MScalarType> > MRowType;

            CRowType tmp_row;
            IndexType nRows(C.nrows());
            for (IndexType row_idx = 0; row_idx < nRows; ++row_idx)
            {
                apply_with_mask(tmp_row, C.getRow(row_idx), Z.getRow(row_idx),
                                mask.getRow(row_idx), replace);

                // Now, set the new one.  Yes, we can optimize this later
                C.setRow(row_idx, tmp_row);
            }
        }

        //**********************************************************************
        // Matrix version specialized for no mask

        template < typename CMatrixT,
                   typename ZMatrixT >
        void write_with_opt_mask(CMatrixT                   &C,
                                 ZMatrixT           const   &Z,
                                 backend::NoMask    const   &foo,
                                 bool                       replace)
        {
            sparse_copy(C, Z);
        }

   } // backend
} // GraphBLAS

#endif
