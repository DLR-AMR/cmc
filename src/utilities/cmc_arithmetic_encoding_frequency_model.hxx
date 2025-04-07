#ifndef CMC_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX
#define CMC_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"
#include "utilities/cmc_iface_arithmetic_encoding_freq_model.hxx"

#include <vector>
#include <algorithm>

namespace cmc::entropy_coding::arithmetic_coding
{

class StaticFrequencyModel;

template <typename Iter>
std::pair<StaticFrequencyModel, size_t>
DecodeStaticFrequencyAlphabet(Iter pos);

constexpr uint32_t kTargetFrequencyCount = 16350;

class StaticFrequencyModel : public IACModel
{
public:
    StaticFrequencyModel() = delete;
    StaticFrequencyModel(std::vector<Letter> alphabet)
    : alphabet_{alphabet}
    {
        /* The alphabet is sorted lexicographically by the symbols */
        std::sort(alphabet_.begin(), alphabet_.end(), [](const Letter& val1, const Letter& val2){return val1.symbol < val2.symbol;});

        size_t total_count_frequencies = 0;

        /* Calculate the cumulative frequency */
        for (auto letter_iter = alphabet_.begin(); letter_iter != alphabet_.end(); ++letter_iter)
        {
            total_count_frequencies += static_cast<size_t>(letter_iter->frequency);
        }

        /* We scale the the frequencies to a certain total */
        const double scaling = static_cast<double>(total_count_frequencies) / static_cast<double>(kTargetFrequencyCount);

        uint32_t scaled_total_frequency = 0;
        /* Scale the frequencies to match the target frequency count */
        for (size_t idx = 0; idx < alphabet_.size(); ++idx)
        {
            if (alphabet_[idx].frequency != 0)
            {
            alphabet_[idx].frequency = static_cast<uint32_t>(static_cast<double>(alphabet_[idx].frequency) / scaling);
            if (alphabet_[idx].frequency == 0) {++alphabet_[idx].frequency;}
            scaled_total_frequency += alphabet_[idx].frequency;
            }
        }

        /* Accumulate the frequencies */
        cumulative_frequencies_.reserve(alphabet_.size() + 1);
        cumulative_frequencies_.push_back(0);
        for (size_t idx = 0; idx < alphabet_.size(); ++idx)
        {
            cumulative_frequencies_.push_back(alphabet_[idx].frequency + cumulative_frequencies_.back());
        }

        /* Validate the cumulative frequencies */
        for (size_t idx = 1; idx < cumulative_frequencies_.size(); ++idx)
        {
            if (cumulative_frequencies_[idx - 1] >= cumulative_frequencies_[idx])
            {
                cmc_err_msg("At least one symbol from the alphabet has a zero frequency which therefore cannot be modeled.");
            }
        }
    };
    ~StaticFrequencyModel(){};

    inline bool SymbolExistsInAlphabet(const uint32_t symbol) const override
    {
        auto iter = GetIteratorToSymbol(symbol);
        return (iter != alphabet_.end() ? true : false);
    }
    uint32_t GetFrequencyCountOf(const uint32_t symbol) override
    {
        auto symbol_iter = GetIteratorToSymbol(symbol);
        return symbol_iter->frequency;
    }
    uint32_t GetTotalSymbolFrequencyCount() override
    {
        return cumulative_frequencies_.back();
    }
    uint32_t GetCumulativeCountOfAllSymbolFrequenciesLowerThan(const uint32_t symbol) override
    {
        cmc_assert(SymbolExistsInAlphabet(symbol));
        const int32_t symbol_idx = GetIndexOfSymbol(symbol);
        return cumulative_frequencies_[symbol_idx];
    }
    uint32_t GetCumulativeCountOfAllSymbolFrequenciesIncluding(const uint32_t symbol) override
    {
        cmc_assert(SymbolExistsInAlphabet(symbol));
        const int32_t symbol_idx = GetIndexOfSymbol(symbol);
        return cumulative_frequencies_[symbol_idx + 1];
    }
    uint32_t GetAlphabetSize() const override
    {
        return static_cast<uint32_t>(alphabet_.size());
    }

    bit_vector::BitVector EncodeAlphabet() const override {
        const uint32_t alphabet_size = GetAlphabetSize();

        /* Allocate memory for the encoded alphabet */
        bit_vector::BitVector encoded_alphabet;
        encoded_alphabet.Reserve(2 * alphabet_size * sizeof(uint32_t) + sizeof(uint32_t));

        /* Store the alphabet size first */
        std::array<uint8_t, sizeof(uint32_t)> serialized_value = SerializeValue(alphabet_size, Endian::Big);
        encoded_alphabet.AppendBytes<sizeof(uint32_t)>(serialized_value);

        /* Afterwards, we store the symbol and the frequency */
        for (size_t idx = 0; idx < alphabet_.size(); ++idx)
        {
            /* Store the symbol */
            serialized_value = SerializeValue(alphabet_[idx].symbol, Endian::Big);
            encoded_alphabet.AppendBytes<sizeof(uint32_t)>(serialized_value);

            /* Store the frequency of the symbol */
            serialized_value = SerializeValue(alphabet_[idx].frequency, Endian::Big);
            encoded_alphabet.AppendBytes<sizeof(uint32_t)>(serialized_value);
        }

        /* Remove the leftover empty byte at the end of the encoding */
        encoded_alphabet.TrimToContent();

        return encoded_alphabet;
    };

    uint32_t GetSymbolFromCumulativeFrequency(const uint32_t value) override
    {
        for (auto iter = std::next(cumulative_frequencies_.begin()); iter != cumulative_frequencies_.end(); ++iter)
        {
            if (*iter > value)
            {
                const int32_t idx = std::distance(cumulative_frequencies_.begin(), iter) - 1;
                return alphabet_[idx].symbol;
            }
        }
        return std::numeric_limits<uint32_t>::max();
    }
    
    template <typename Iter> friend std::pair<StaticFrequencyModel, size_t> DecodeStaticFrequencyAlphabet(Iter pos);
private:
    struct Decoding{};
    StaticFrequencyModel(std::vector<Letter>&& alphabet, [[maybe_unused]] Decoding&&)
    : alphabet_{std::move(alphabet)}
    {
        /* Accumulate the frequencies */
        cumulative_frequencies_.reserve(alphabet_.size() + 1);
        cumulative_frequencies_.push_back(0);
        for (size_t idx = 0; idx < alphabet_.size(); ++idx)
        {
            cumulative_frequencies_.push_back(alphabet_[idx].frequency + cumulative_frequencies_.back());
        }

        /* Validate the cumulative frequencies */
        for (size_t idx = 1; idx < cumulative_frequencies_.size(); ++idx)
        {
            if (cumulative_frequencies_[idx - 1] >= cumulative_frequencies_[idx])
            {
                cmc_err_msg("At least one symbol from the alphabet has a zero frequency which therefore cannot be modeled.");
            }
        }
    };

    inline std::vector<Letter>::const_iterator GetIteratorToSymbol(const uint32_t symbol) const 
    {
        std::vector<Letter>::const_iterator symbol_iter = std::lower_bound(alphabet_.begin(), alphabet_.end(), Letter{symbol, 0}, [](const Letter& val1, const Letter& val2){return val1.symbol < val2.symbol;});
        return symbol_iter;
    }
    
    inline int32_t GetIndexOfSymbol(const uint32_t symbol) const
    {
        auto symbol_iter = GetIteratorToSymbol(symbol);
        if (symbol_iter == alphabet_.end()) {return -1;}
        return std::distance(alphabet_.begin(), symbol_iter);
    }

    std::vector<Letter> alphabet_;
    std::vector<uint32_t> cumulative_frequencies_;
};

template <typename Iter>
std::pair<StaticFrequencyModel, size_t>
DecodeStaticFrequencyAlphabet(Iter pos)
{
    using symbol_type_t = uint32_t;
    using frequency_type_t = uint32_t;
    const uint32_t symbol_type_length = sizeof(uint32_t);
    const uint32_t frequency_type_length = sizeof(uint32_t);

    std::array<uint8_t, sizeof(uint32_t)> serialized_uint32_t;
    uint32_t decoded_byte_count = 0;
    uint32_t& offset = decoded_byte_count;

    /* Get the size of the alphabet */
    serialized_uint32_t = DeserializeValue<uint32_t>(pos, Endian::Big);
    uint32_t _num_letters;
    std::memcpy(static_cast<void*>(&_num_letters), static_cast<const void*>(serialized_uint32_t.data()), sizeof(uint32_t));
    const uint32_t count_letters = _num_letters;

    offset += sizeof(uint32_t);

    std::vector<Letter> alphabet;
    alphabet.reserve(count_letters);

    for (uint32_t idx = 0; idx < count_letters; ++idx)
    {
        /* Ceate a new letter */
        alphabet.emplace_back();

        /* Get the symbol from the encoded stream */
        serialized_uint32_t = DeserializeValue<symbol_type_t>(pos + offset, Endian::Big);
        std::memcpy(static_cast<void*>(&alphabet.back().symbol), static_cast<const void*>(serialized_uint32_t.data()), sizeof(symbol_type_length));
        offset += symbol_type_length;

        /* Get the frequency from the encoded stream */
        serialized_uint32_t = DeserializeValue<frequency_type_t>(pos + offset, Endian::Big);
        std::memcpy(static_cast<void*>(&alphabet.back().frequency), static_cast<const void*>(serialized_uint32_t.data()), sizeof(frequency_type_length));
        offset += frequency_type_length;
    }

    return std::make_pair(StaticFrequencyModel(std::move(alphabet), StaticFrequencyModel::Decoding()), decoded_byte_count);
}

}

#endif /* !CMC_ARITHMETIC_ENCODING_FREQUENCY_MODEL_HXX */
