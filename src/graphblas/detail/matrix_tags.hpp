#pragma once

namespace GraphBLAS
{
    // The default matrix is sparse and directed, and the default vector is sparse,
    // so we need tags that modify that
    struct DirectedMatrixTag {};
    struct UndirectedMatrixTag {};
    struct DenseTag {};
    struct SparseTag {};

    namespace detail
    {
        // add category tags in the detail namespace
        struct SparsenessCategoryTag {};
        struct DirectednessCategoryTag {};
        struct NullTag {};
    } //end detail

}//end GraphBLAS

