#ifndef CMC_ENTROPY_ALPHABET_HXX
#define CMC_ENTROPY_ALPHABET_HXX

#include <vector>

namespace cmc
{

struct Letter
{
    uint32_t symbol;
    uint32_t frequency;
};

typedef std::vector<Letter> Alphabet;

class IEntropyAlphabet
{
    virtual void InitializeSymbols() = 0;
    virtual void UpdateSymbol(const uint32_t symbol) = 0;
    virtual Alphabet GetSymbolFrequencies() = 0;
    virtual ~IEntropyAlphabet(){};
}

}


#endif /*! CMC_ENTROPY_ALPHABET_HXX */
