#ifndef CMC_ENTROPY_ALPHABET_HXX
#define CMC_ENTROPY_ALPHABET_HXX

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
    virtual Alphabet GetSymbolFrequencies() = 0;
    virtual ~IEntropyAlphabet(){};
};

}


#endif /*! CMC_ENTROPY_ALPHABET_HXX */
