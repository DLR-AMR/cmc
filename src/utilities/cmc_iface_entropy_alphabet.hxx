#ifndef CMC_IFACE_ENTROPY_ALPHABET_HXX
#define CMC_IFACE_ENTROPY_ALPHABET_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"

#include <unordered_map>

namespace cmc::entropy_coding
{

struct Letter
{
    uint32_t symbol;
    uint32_t frequency;
};

typedef std::unordered_map<uint32_t, uint32_t> Alphabet;

class IEntropyAlphabet
{
public:
    virtual void InitializeSymbols(const size_t type_size_in_bytes) = 0;
    virtual void UpdateSymbol(const uint32_t symbol) = 0;
    virtual ~IEntropyAlphabet(){};
};

class IByteCompressionEntropyAlphabet
{
public:
    virtual void InitializeSymbols(const size_t type_size_in_bytes) = 0;
    virtual void UpdateSymbol(const uint32_t symbol) = 0;
    virtual size_t GetAlphabetSize() const = 0;
    virtual std::vector<uint32_t> CommunicateSymbolFrequencies(const MPI_Comm comm) = 0;
    virtual const std::vector<uint32_t>& GetSymbolFrequencies() const = 0;
    virtual uint32_t RevertArrayIndexToSymbol(const uint32_t index) = 0;
    virtual uint32_t TransformSymbolToArrayIndex(const uint32_t symbol) = 0;
    virtual bool DoesSymbolExistInAlphabet(const uint32_t symbol) = 0;

    virtual ~IByteCompressionEntropyAlphabet(){};
};

namespace arithmetic_coding
{

constexpr uint32_t kByteCompressionSymbolJumpToNextByte = 0x80000000;

constexpr inline bool
CheckIfJumpToNextByteIndicatorBitIsSet(const uint32_t symbol)
{
    return (symbol == kByteCompressionSymbolJumpToNextByte);
}

}
}


#endif /*! CMC_IFACE_ENTROPY_ALPHABET_HXX */
