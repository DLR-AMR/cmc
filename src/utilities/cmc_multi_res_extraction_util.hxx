#ifndef CMC_MULTI_RES_EXTRACTION_UTIL_HXX
#define CMC_MULTI_RES_EXTRACTION_UTIL_HXX


#include <utility>
#include <cstddef>

namespace cmc::lossless::multi::res::util
{

constexpr uint32_t kMSBBit = 0x80000000;

enum ResidualOperation {IntegerAddition, IntegerSubtraction};

inline
uint32_t
GetSignumForEncoding(const bool is_approximation_greater)
{
    if (is_approximation_greater)
    {
        return 0;
    } else 
    {
        return kMSBBit;
    }
}

inline
std::pair<ResidualOperation, uint32_t>
DecodeLZC(const uint32_t encoded_lzc)
{
    if(kMSBBit > encoded_lzc)
    {
        return std::make_pair(ResidualOperation::IntegerSubtraction, encoded_lzc);
    } else
    {
        return std::make_pair(ResidualOperation::IntegerAddition, encoded_lzc - kMSBBit);
    }
}


}

#endif /* !CMC_MULTI_RES_EXTRACTION_UTIL_HXX */
