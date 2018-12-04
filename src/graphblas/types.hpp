#ifndef GB_TYPES_HPP
#define GB_TYPES_HPP

#include <cstdint>
#include <exception>
#include <vector>
#include <iostream>

namespace GraphBLAS
{
    typedef uint64_t IndexType;
    typedef std::vector<IndexType> IndexArrayType;

    //**************************************************************************
    // Error codes to replace some of the exceptions
    enum Info {
        SUCCESS,
        INVALID_VALUE,
        INVALID_INDEX,
        DOMAIN_MISMATCH,
        DIMENSION_MISMATCH,
        OUTPUT_NOT_EMPTY,
        NO_VALUE,
        INDEX_OUT_OF_BOUNDS,
        OUT_OF_MEMORY
    };


    //**************************************************************************
    struct NoAccumulate
    {
        // It doesn't really matter what the type is, it never gets executed.
        typedef bool result_type;
        inline bool operator()(bool lhs, bool rhs) { return true; }
    };

    struct matrix_tag {};

} 

#define CHECK_STATUS(x)     do { Info status(x); if (status != SUCCESS) return status; } while(0)

namespace std
{

// @TODO; It seems that unit tests can't find this!
    inline std::ostream &
    operator<<(std::ostream &os, const std::vector<long unsigned int> vec)
    {
        bool first = true;
        for (auto it = vec.begin(); it != vec.end(); ++it)
        {
            os << (first ? "" : ",") << *it;
            first = false;
        }
        return os;
    }

} // namespace std

#endif // GB_TYPES_HPP
