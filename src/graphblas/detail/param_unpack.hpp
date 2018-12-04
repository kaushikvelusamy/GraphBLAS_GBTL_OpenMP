#pragma once
#include "matrix_tags.hpp"

#define GB_INCLUDE_BACKEND_MATRIX 1

#include <backend_include.hpp>

//this file contains the variadic template parameters unpacking utility.


namespace GraphBLAS
{
    namespace detail
    {

        // Substitute template to decide if a tag goes into a given slot
        template<typename TagCategory, typename Tag>
        struct substitute {
            using type = TagCategory;
        };


        template<>
        struct substitute<detail::SparsenessCategoryTag, DenseTag> {
            using type = DenseTag;
        };

        template<>
        struct substitute<detail::SparsenessCategoryTag, SparseTag> {
            using type = SparseTag;
        };

        template<>
        struct substitute<detail::DirectednessCategoryTag, UndirectedMatrixTag> {
            using type = UndirectedMatrixTag;
        };

        template<>
        struct substitute<detail::DirectednessCategoryTag, DirectedMatrixTag> {
            using type = DirectedMatrixTag;
        };

        template<>
        struct substitute<detail::DirectednessCategoryTag, detail::NullTag> {
            //default values
            using type = DirectedMatrixTag; // default directedness
        };

        template<>
        struct substitute<detail::SparsenessCategoryTag, detail::NullTag> {
            using type = SparseTag; // default sparseness
        };


        // hidden part in the frontend (detail namespace somewhere) to unroll
        // template parameter pack

        struct matrix_generator {
            // recursive call: shaves off one of the tags and puts it in the right
            // place (no error checking yet)
            template<typename ScalarT, typename Sparseness, typename Directedness,
                typename InputTag, typename... Tags>
            struct result {
                using type = typename result<ScalarT,
                      typename detail::substitute<Sparseness, InputTag >::type,
                      typename detail::substitute<Directedness, InputTag >::type,
                      Tags... >::type;
            };

            //null tag shortcut:
            template<typename ScalarT, typename Sparseness, typename Directedness>
            struct result<ScalarT, Sparseness, Directedness, detail::NullTag, detail::NullTag>
            {
                using type = typename backend::Matrix<ScalarT,
                      typename detail::substitute<Sparseness, detail::NullTag >::type,
                      typename detail::substitute<Directedness, detail::NullTag >::type >;
            };

            // base case returns the matrix from the backend
            template<typename ScalarT, typename Sparseness, typename Directedness, typename InputTag>
            struct result<ScalarT, Sparseness, Directedness, InputTag>
            {
                using type = typename backend::Matrix<ScalarT,
                      typename detail::substitute<Sparseness, InputTag >::type,
                      typename detail::substitute<Directedness, InputTag >::type > ;
            };
        };


    }//end detail
}
