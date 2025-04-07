#ifndef CMC_MULTI_RES_EXTRACTION_UTIL_HXX
#define CMC_MULTI_RES_EXTRACTION_UTIL_HXX

#include "utilities/cmc_byte_compression_entropy_alphabet.hxx"

#include <utility>
#include <cstddef>

namespace cmc::lossless::multi_res::util
{

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
        return cmc::entropy_coding::arithmetic_coding::kByteCompressionSignumBit;
    }
}

inline
std::pair<ResidualOperation, uint32_t>
DecodeLZC(const uint32_t encoded_lzc)
{
    if(cmc::entropy_coding::arithmetic_coding::kByteCompressionSignumBit > encoded_lzc)
    {
        return std::make_pair(ResidualOperation::IntegerSubtraction, encoded_lzc);
    } else
    {
        return std::make_pair(ResidualOperation::IntegerAddition, encoded_lzc - cmc::entropy_coding::arithmetic_coding::kByteCompressionSignumBit);
    }
}


}

#endif /* !CMC_MULTI_RES_EXTRACTION_UTIL_HXX */
