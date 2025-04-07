#ifndef CMC_BYTE_COMPRESSION_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX
#define CMC_BYTE_COMPRESSION_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"
#include "utilities/cmc_byte_compression_entropy_alphabet.hxx"
#include "utilities/cmc_iface_arithmetic_encoding_freq_model.hxx"

#include <vector>
#include <algorithm>
#include <array>
#include <numeric>

namespace cmc::entropy_coding::arithmetic_coding
{

template <typename T>
class ByteCompressionStaticFrequencyModel;

template <typename T, typename Iter>
std::pair<ByteCompressionStaticFrequencyModel<T>, size_t>
DecodeByteCompressionStaticFrequencyAlphabet(Iter pos);

constexpr uint32_t kTargetFrequencyCount = 16350;


template <typename T>
void hello();

template <typename T>
class ByteCompressionStaticFrequencyModel : public IACModel
{
public:
    constexpr static size_t N = ByteCompressionAlphabet<T>::N;

    ByteCompressionStaticFrequencyModel() = delete;
    ByteCompressionStaticFrequencyModel(const std::array<uint32_t, N>& symbol_frequencies)
    {
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
            cumulative_frequencies_[idx] = previous_scaled_cumul_frequencies + current_scaled_freq;

            /* Update the frequency count */
            previous_scaled_cumul_frequencies = cumulative_frequencies_[idx];
        }
    };

    ~ByteCompressionStaticFrequencyModel(){};

    inline bool SymbolExistsInAlphabet(const uint32_t symbol) const override
    {
        return DoesSymbolExistsInAlphabet<T>(symbol);
    }
    uint32_t GetFrequencyCountOf(const uint32_t symbol) override
    {
        cmc_assert(SymbolExistsInAlphabet(symbol));
        const uint32_t symbol_idx = TransformSymbolToArrayIndex<T>(symbol);
        if (symbol_idx != 0)
        {
            return cumulative_frequencies_[symbol_idx] - cumulative_frequencies_[symbol_idx - 1];
        } else
        {
            return cumulative_frequencies_[symbol_idx];
        }
    }
    uint32_t GetTotalSymbolFrequencyCount() override
    {
        return cumulative_frequencies_.back();
    }
    uint32_t GetCumulativeCountOfAllSymbolFrequenciesLowerThan(const uint32_t symbol) override
    {
        cmc_assert(SymbolExistsInAlphabet(symbol));
        const uint32_t symbol_idx = GetIndexOfSymbol(symbol);
        if (symbol_idx != 0)
        {
            return cumulative_frequencies_[symbol_idx - 1];
        } else
        {
            return 0;
        }
    }
    uint32_t GetCumulativeCountOfAllSymbolFrequenciesIncluding(const uint32_t symbol) override
    {
        cmc_assert(SymbolExistsInAlphabet(symbol));
        const uint32_t symbol_idx = GetIndexOfSymbol(symbol);
        return cumulative_frequencies_[symbol_idx];
    }
    uint32_t GetAlphabetSize() const override
    {
        return static_cast<uint32_t>(N);
    }

    bit_vector::BitVector EncodeAlphabet() const override {
        /* Allocate memory for the encoded alphabet */
        bit_vector::BitVector encoded_alphabet;
        encoded_alphabet.Reserve((N + 1) * sizeof(uint32_t));

        /* Store the number of symbols wihtin the alphabet first */
        std::array<uint8_t, sizeof(uint32_t)> serialized_value = SerializeValue(static_cast<uint32_t>(N), Endian::Big);
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

    uint32_t GetSymbolFromCumulativeFrequency(const uint32_t cumul_freq) override
    {
        auto sym_iter = std::lower_bound(cumulative_frequencies_.begin(), cumulative_frequencies_.end(), cumul_freq);
        cmc_assert(sym_iter != cumulative_frequencies_.end());

        const uint32_t idx = std::distance(cumulative_frequencies_.begin(), sym_iter);

        return RevertArrayIndexToSymbol<T>(idx);
    }
    
private:
    struct Decoding{};
    ByteCompressionStaticFrequencyModel(const std::array<uint32_t, N>& cumulative_frequencies, [[maybe_unused]] Decoding&&)
    : cumulative_frequencies_{cumulative_frequencies} {};
    
    template <typename U, typename Iter> friend std::pair<ByteCompressionStaticFrequencyModel<U>, size_t> DecodeByteCompressionStaticFrequencyAlphabet(Iter pos);

    inline uint32_t GetIndexOfSymbol(const uint32_t symbol) const
    {
        cmc_assert(SymbolExistsInAlphabet(symbol));
        return TransformSymbolToArrayIndex<T>(symbol);
    }

    std::array<uint32_t, N> cumulative_frequencies_;
};

template <typename T, typename Iter>
std::pair<ByteCompressionStaticFrequencyModel<T>, size_t>
DecodeByteCompressionStaticFrequencyAlphabet(Iter pos)
{
    using frequency_type_t = uint32_t;
    const size_t frequency_type_length = sizeof(uint32_t);

    std::array<uint8_t, frequency_type_length> serialized_value;

    uint32_t decoded_byte_count = 0;
    uint32_t& offset = decoded_byte_count;

    /* Get the size of the alphabet */
    serialized_value = DeserializeValue<uint32_t>(pos, Endian::Big);
    uint32_t _num_symbol_frequencies;
    std::memcpy(static_cast<void*>(&_num_symbol_frequencies), static_cast<const void*>(serialized_value.data()), sizeof(uint32_t));
    const uint32_t num_symbol_frequencies = _num_symbol_frequencies;

    offset += sizeof(uint32_t);

    std::array<frequency_type_t, ByteCompressionAlphabet<T>::N> cumulative_symbol_frequencies;

    cmc_assert(static_cast<size_t>(num_symbol_frequencies) == cumulative_symbol_frequencies.size());

    /* Read in all sybol frequencies */
    for (size_t sym_idx = 0; sym_idx < cumulative_symbol_frequencies.size(); ++sym_idx)
    {
        /* Get the next frequency */
        serialized_value = DeserializeValue<frequency_type_t>(pos + offset, Endian::Big);
        std::memcpy(static_cast<void*>(&cumulative_symbol_frequencies[sym_idx]), static_cast<const void*>(serialized_value.data()), frequency_type_length);
        offset += frequency_type_length;
    }

    return std::make_pair(ByteCompressionStaticFrequencyModel<T>(cumulative_symbol_frequencies, typename ByteCompressionStaticFrequencyModel<T>::Decoding()), decoded_byte_count);
}

}

#endif /* !CMC_BYTE_COMPRESSION_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX */
