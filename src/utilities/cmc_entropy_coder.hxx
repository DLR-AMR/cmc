#ifndef CMC_ENTROPY_CODER_HXX
#define CMC_ENTROPY_CODER_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"
#include "utilities/cmc_iface_arithmetic_encoding_freq_model.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#endif

#include <memory>

namespace cmc::entropy_coding
{

class IEntropyCoder
{
public: 
    virtual void InitializeAlphabet(const size_t type_size_in_bytes) = 0;
    virtual void UpdateSymbolFrequency(const uint32_t symbol) = 0;
    virtual bit_vector::BitVector EncodeAlphabet() const = 0;

    virtual void SetupEncoding() = 0;
#ifdef CMC_ENABLE_MPI
    virtual void SetupEncoding(const MPI_Comm comm) = 0;
#endif
    virtual void EncodeSymbol(const uint32_t symbol) = 0;
    virtual void FinishEncoding() = 0;

    virtual bit_map::BitMap GetEncodedBitStream() const = 0;
    virtual void ClearEncodedBitStream() = 0;

    virtual void Reset() = 0;
    
    virtual ~IEntropyCoder(){};
protected:
    std::unique_ptr<IEntropyAlphabet> alphabet_;
    std::unique_ptr<IACModel> frequency_model_;
};

class IByteCompressionEntropyCoder
{
public: 
    virtual void InitializeAlphabet(const size_t type_size_in_bytes) = 0;
    virtual void UpdateSymbolFrequency(const uint32_t symbol) = 0;
    virtual bit_vector::BitVector EncodeAlphabet() const = 0;

    virtual void SetupEncoding() = 0;
#ifdef CMC_ENABLE_MPI
    virtual void SetupEncoding(const MPI_Comm comm) = 0;
#endif
    virtual void EncodeSymbol(const uint32_t symbol) = 0;
    virtual void FinishEncoding() = 0;

    virtual bit_map::BitMap GetEncodedBitStream() const = 0;
    virtual void ClearEncodedBitStream() = 0;

    virtual void Reset(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet) = 0;
    
    virtual ~IByteCompressionEntropyCoder(){};
protected:
    std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet_;
    std::unique_ptr<IACModel> frequency_model_;
};

}

#endif /* !CMC_ENTROPY_CODER_HXX */
