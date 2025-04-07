#ifndef CMC_IFACE_ARITHMETIC_ENCODING_FREQ_MODEL_HXX
#define CMC_IFACE_ARITHMETIC_ENCODING_FREQ_MODEL_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_bit_vector.hxx"

#include <cstddef>

namespace cmc::entropy_coding
{

class IACModel
{
public:
    virtual uint32_t GetFrequencyCountOf(const uint32_t symbol) = 0;
    virtual uint32_t GetTotalSymbolFrequencyCount() = 0;
    virtual uint32_t GetCumulativeCountOfAllSymbolFrequenciesLowerThan(const uint32_t symbol) = 0;
    virtual uint32_t GetCumulativeCountOfAllSymbolFrequenciesIncluding(const uint32_t symbol) = 0;
    virtual uint32_t GetAlphabetSize() const = 0;
    virtual bit_vector::BitVector EncodeAlphabet() const = 0;
    virtual uint32_t GetSymbolFromCumulativeFrequency(const uint32_t value) = 0;
    virtual bool SymbolExistsInAlphabet(const uint32_t symbol) const = 0;
    virtual ~IACModel(){};
};

}

#endif /* !CMC_IFACE_ARITHMETIC_ENCODING_FREQ_MODEL_HXX */