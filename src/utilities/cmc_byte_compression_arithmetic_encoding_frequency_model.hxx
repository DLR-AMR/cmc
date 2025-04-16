#ifndef CMC_BYTE_COMPRESSION_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX
#define CMC_BYTE_COMPRESSION_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"
#include "utilities/cmc_iface_arithmetic_encoding_freq_model.hxx"

#include <vector>
#include <algorithm>
#include <array>
#include <numeric>
#include <memory>

namespace cmc::entropy_coding::arithmetic_coding
{

template <typename T>
class ByteCompressionStaticFrequencyModel;

template <typename T, typename Iter>
std::unique_ptr<ByteCompressionStaticFrequencyModel<T>>
DecodeByteCompressionStaticFrequencyAlphabet(Iter pos, std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet);

constexpr uint32_t kTargetFrequencyCount = 16350;

template <typename T>
class ByteCompressionStaticFrequencyModel : public IACModel
{
private:
    struct Decoding{};
public:

    ByteCompressionStaticFrequencyModel() = delete;
    /* Constructors for encoding */
    ByteCompressionStaticFrequencyModel(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet, const std::vector<uint32_t>& symbol_frequencies)
    : alphabet_{std::move(alphabet)} {
        ComputeCumulativeFrequencies(symbol_frequencies);
    };
    /* Constructors for decoding */
    ByteCompressionStaticFrequencyModel(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet, const std::vector<uint32_t>& cumulative_frequencies, [[maybe_unused]] Decoding)
    : alphabet_{std::move(alphabet)}, cumulative_frequencies_{cumulative_frequencies} {};
    ByteCompressionStaticFrequencyModel(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet, std::vector<uint32_t>&& cumulative_frequencies, [[maybe_unused]] Decoding)
    : alphabet_{std::move(alphabet)}, cumulative_frequencies_{std::move(cumulative_frequencies)} {};

    ~ByteCompressionStaticFrequencyModel(){};

    uint32_t GetFrequencyCountOf(const uint32_t symbol) override;
    uint32_t GetTotalSymbolFrequencyCount() override {return cumulative_frequencies_.back();}
    uint32_t GetCumulativeCountOfAllSymbolFrequenciesLowerThan(const uint32_t symbol) override;
    uint32_t GetCumulativeCountOfAllSymbolFrequenciesIncluding(const uint32_t symbol) override;
    uint32_t GetAlphabetSize() const override {cmc_assert(alphabet_ != nullptr); return static_cast<uint32_t>(alphabet_->GetAlphabetSize());}
    bit_vector::BitVector EncodeAlphabet() const override;

    uint32_t GetSymbolFromCumulativeFrequency(const uint32_t cumul_freq) override;
    
private:
    template <typename U, typename Iter> friend std::unique_ptr<ByteCompressionStaticFrequencyModel<U>> DecodeByteCompressionStaticFrequencyAlphabet(Iter pos, std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet);

    inline bool SymbolExistsInAlphabet(const uint32_t symbol) const override {cmc_assert(alphabet_ != nullptr); return alphabet_->DoesSymbolExistInAlphabet(symbol);}
    inline uint32_t GetSymbolIndex(const uint32_t symbol) {cmc_assert(SymbolExistsInAlphabet(symbol)); return alphabet_->TransformSymbolToArrayIndex(symbol);}
    inline uint32_t RevertToSymbol(const uint32_t index) {cmc_assert(alphabet_ != nullptr); return alphabet_->RevertArrayIndexToSymbol(index);}
    void ComputeCumulativeFrequencies(const std::vector<uint32_t>& symbol_frequencies);

    std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet_{nullptr};
    std::vector<uint32_t> cumulative_frequencies_;
};

template <typename T>
uint32_t
ByteCompressionStaticFrequencyModel<T>::GetSymbolFromCumulativeFrequency(const uint32_t cumul_freq) 
{
    auto sym_iter = std::upper_bound(cumulative_frequencies_.begin(), cumulative_frequencies_.end(), cumul_freq);
    cmc_assert(sym_iter != cumulative_frequencies_.end());
    const uint32_t idx = std::distance(cumulative_frequencies_.begin(), sym_iter);
    return RevertToSymbol(idx);
}

template <typename T>
bit_vector::BitVector
ByteCompressionStaticFrequencyModel<T>::EncodeAlphabet() const
{
    /* Allocate memory for the encoded alphabet */
    bit_vector::BitVector encoded_alphabet;
    encoded_alphabet.Reserve((alphabet_->GetAlphabetSize() + 1) * sizeof(uint32_t));
    /* Store the number of symbols wihtin the alphabet first */
    std::array<uint8_t, sizeof(uint32_t)> serialized_value = SerializeValue(static_cast<uint32_t>(alphabet_->GetAlphabetSize()), Endian::Big);
    encoded_alphabet.AppendBytes<sizeof(uint32_t)>(serialized_value);
    /* We store the cumulative frequencies in the given order afterwards */
    for (auto cumul_freq_iter = cumulative_frequencies_.begin(); cumul_freq_iter != cumulative_frequencies_.end(); ++cumul_freq_iter)
    {
        serialized_value = SerializeValue(*cumul_freq_iter, Endian::Big);
        encoded_alphabet.AppendBytes<sizeof(uint32_t)>(serialized_value);
    }
    /* Remove the (potential) leftover empty byte at the end of the encoding */
    encoded_alphabet.TrimToContent();
    return encoded_alphabet;
};

template <typename T>
uint32_t
ByteCompressionStaticFrequencyModel<T>::GetCumulativeCountOfAllSymbolFrequenciesIncluding(const uint32_t symbol)
{
    cmc_assert(SymbolExistsInAlphabet(symbol));
    const uint32_t symbol_idx = GetSymbolIndex(symbol);
    return cumulative_frequencies_[symbol_idx];
}

template <typename T>
uint32_t
ByteCompressionStaticFrequencyModel<T>::GetCumulativeCountOfAllSymbolFrequenciesLowerThan(const uint32_t symbol)
{
    cmc_assert(SymbolExistsInAlphabet(symbol));
    const uint32_t symbol_idx = GetSymbolIndex(symbol);
    if (symbol_idx != 0)
    {
        return cumulative_frequencies_[symbol_idx - 1];
    } else
    {
        return 0;
    }
}

template <typename T>
uint32_t
ByteCompressionStaticFrequencyModel<T>::GetFrequencyCountOf(const uint32_t symbol)
{
    cmc_assert(SymbolExistsInAlphabet(symbol));
    const uint32_t symbol_idx = GetSymbolIndex(symbol);
    if (symbol_idx != 0)
    {
        return cumulative_frequencies_[symbol_idx] - cumulative_frequencies_[symbol_idx - 1];
    } else
    {
        return cumulative_frequencies_[symbol_idx];
    }
}

template <typename T>
void
ByteCompressionStaticFrequencyModel<T>::ComputeCumulativeFrequencies(const std::vector<uint32_t>& symbol_frequencies)
{
    /* Allocate memory for the cumulative frequency counts */
    cumulative_frequencies_.reserve(symbol_frequencies.size());

    /* Get the total frequency count */
    const size_t total_count_frequencies = std::accumulate(symbol_frequencies.begin(), symbol_frequencies.end(), 0);

    /* We scale the the frequencies to a certain total */
    const double scaling = static_cast<double>(total_count_frequencies) / static_cast<double>(kTargetFrequencyCount);

    /* Scale and accumulate the frequencies */
    uint32_t previous_scaled_cumul_frequencies = 0;
    for (size_t idx = 0; idx < symbol_frequencies.size(); ++idx)
    {
        uint32_t current_scaled_freq = symbol_frequencies[idx];
        if (current_scaled_freq != 0)
        {
            /* Scale the current frequency */
            current_scaled_freq = static_cast<uint32_t>(static_cast<double>(symbol_frequencies[idx]) / scaling);
        
            /* Check if frequency became zero */
            if (current_scaled_freq == 0) {current_scaled_freq = 1;}
        }

        /* Accumulate the frequency with the previous one */
        cumulative_frequencies_.push_back(previous_scaled_cumul_frequencies + current_scaled_freq);

        /* Update the frequency count */
        previous_scaled_cumul_frequencies = cumulative_frequencies_.back();
    }
}

template <typename T, typename Iter>
std::unique_ptr<ByteCompressionStaticFrequencyModel<T>>
DecodeByteCompressionStaticFrequencyAlphabet(Iter pos, std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet)
{
    using frequency_type_t = uint32_t;
    const size_t frequency_type_length = sizeof(uint32_t);

    std::array<uint8_t, frequency_type_length> serialized_value;

    uint32_t offset = 0;

    /* Get the size of the alphabet */
    serialized_value = DeserializeValue<uint32_t>(pos, Endian::Big);
    uint32_t _num_symbol_frequencies;
    std::memcpy(static_cast<void*>(&_num_symbol_frequencies), static_cast<const void*>(serialized_value.data()), sizeof(uint32_t));
    const uint32_t num_symbol_frequencies = _num_symbol_frequencies;

    offset += sizeof(uint32_t);

    /* Allocate a vector for the cumulative frequencies */
    std::vector<frequency_type_t> cumulative_symbol_frequencies(alphabet->GetAlphabetSize());

    cmc_assert(static_cast<size_t>(num_symbol_frequencies) == cumulative_symbol_frequencies.size());

    /* Read in all sybol frequencies */
    for (size_t sym_idx = 0; sym_idx < cumulative_symbol_frequencies.size(); ++sym_idx)
    {
        /* Get the next frequency */
        serialized_value = DeserializeValue<frequency_type_t>(pos + offset, Endian::Big);
        std::memcpy(static_cast<void*>(&cumulative_symbol_frequencies[sym_idx]), static_cast<const void*>(serialized_value.data()), frequency_type_length);
        offset += frequency_type_length;
    }

    typename ByteCompressionStaticFrequencyModel<T>::Decoding decoding_flag;

    return std::make_unique<ByteCompressionStaticFrequencyModel<T>>(std::move(alphabet), cumulative_symbol_frequencies, decoding_flag);
}

}

#endif /* !CMC_BYTE_COMPRESSION_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX */
