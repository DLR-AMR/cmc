#ifndef CMC_MULTI_RES_PAR_EXTRACTION_UTIL_HXX
#define CMC_MULTI_RES_PAR_EXTRACTION_UTIL_HXX

#include "utilities/cmc_bit_map.hxx"

#include <cstdint>

namespace cmc::lossless::par::multi_res
{

template <typename T>
constexpr uint32_t kIndicateProcEnd = 2 * sizeof(T) * bit_map::kCharBit + 1;

template <typename T>
constexpr size_t
GetNumEntropySymbols()
{
    return 2 * sizeof(T) * bit_map::kCharBit + 1 + 1;
}

template <typename T>
inline uint32_t
ConvertToSymbolInFrequencyTable(const bool residual_indication, const uint32_t value)
{
    uint32_t symbol = value;
    
    if (residual_indication)
    {
        symbol += sizeof(T) * bit_map::kCharBit + 1;
    }

    cmc_assert(symbol <= 2 * sizeof(T) * bit_map::kCharBit + 1);

    return symbol;
}

template <typename T>
constexpr inline
uint32_t GetProcessEndSymbol()
{
    return 2 * sizeof(T) * bit_map::kCharBit + 1;
}

template <typename T>
constexpr inline
uint32_t GetFullyExtractedSymbol()
{
    return sizeof(T) * bit_map::kCharBit;
}


template <typename T>
constexpr inline std::pair<bool, uint32_t>
ConvertFrequencySymbolToValue(const uint32_t freq_symbol)
{
    cmc_static_assert(freq_symbol <= 2 * sizeof(T) * bit_map::kCharBit + 1);

    /* Check if it is the process end symbol */
    if (freq_symbol == GetProcessEndSymbol<T>())
    {
        return std::make_pair(false, GetProcessEndSymbol<T>());
    }

    /* Check the residual indication */
    if (freq_symbol < sizeof(T) * bit_map::kCharBit + 1)
    {
        return std::make_pair(false, freq_symbol);
    } else
    {
        return std::make_pair(true, freq_symbol - sizeof(T) * bit_map::kCharBit + 1);
    }
}



}

#endif /* !CMC_MULTI_RES_PAR_EXTRACTION_UTIL_HXX */
