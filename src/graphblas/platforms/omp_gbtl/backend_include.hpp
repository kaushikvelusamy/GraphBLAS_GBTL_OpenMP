#if(GB_INCLUDE_BACKEND_ALL)
#include <graphblas/platforms/omp_gbtl/omp_gbtl.hpp>
#endif

#if(GB_INCLUDE_BACKEND_MATRIX)
#include <graphblas/platforms/omp_gbtl/Matrix.hpp>
#undef GB_INCLUDE_BACKEND_MATRIX
#endif

#if(GB_INCLUDE_BACKEND_TRANSPOSE_VIEW)
#include <graphblas/platforms/omp_gbtl/TransposeView.hpp>
#undef GB_INCLUDE_BACKEND_TRANSPOSE_VIEW
#endif

#if(GB_INCLUDE_BACKEND_OPERATIONS)
#include <graphblas/platforms/omp_gbtl/operations.hpp>
#undef GB_INCLUDE_BACKEND_OPERATIONS
#endif
