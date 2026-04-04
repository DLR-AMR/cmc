#ifndef CMC_HUFFMAN_CODER_HXX
#define CMC_HUFFMAN_CODER_HXX

#include "cmc.hxx"
#include "utilities/cmc_bits.hxx"
#include "utilities/cmc_bits_vector.hxx"

#include <vector>
#include <cmath>
#include <unordered_map>
#include <queue>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <algorithm>

namespace cmc::entropy_coding::huffman
{

template <typename T>
class HuffmanTree;

using DefaultSymbolType = int32_t;

using FrequencyType = uint32_t;
using HuffmanCodeInfoType = uint32_t;

using HuffmanCodeWord = uint64_t;
using HuffmanCodeLength = uint64_t;

struct HuffmanCode
{
    HuffmanCodeWord code_word;
    HuffmanCodeLength code_length;
};

template <typename T>
using HuffmanCodeMap = std::unordered_map<T, HuffmanCode>;

template <typename T>
using HuffmanDecodeMap = std::unordered_map<HuffmanCodeWord, T>;

template <typename T>
struct EntropySymbol
{
    EntropySymbol(const T symbol_, const FrequencyType frequency_)
    : symbol{symbol_}, frequency{frequency_} {}

    T symbol;
    FrequencyType frequency;
};

//The codes are permitted to be of a maximum length of 56 bits, because one byte is needed for the code length during encoding
inline void 
AppendSetBit(HuffmanCode& huff_code)
{
    if (huff_code.code_length >= 56) [[unlikely]]
    {
        cmc_err_msg("The code word cannot be larger than 64 bit!");
    }

    huff_code.code_word <<= 1;
    huff_code.code_word |= HuffmanCodeWord{1};
    ++(huff_code.code_length);
}

//The codes are permitted to be of a maximum length of 56 bits, because one byte is needed for the code length during encoding
inline void 
AppendUnsetBit(HuffmanCode& huff_code)
{
    if (huff_code.code_length >= 56) [[unlikely]]
    {
        cmc_err_msg("The code word cannot be larger than 64 bit!");
    }

    huff_code.code_word <<= 1;
    ++(huff_code.code_length);
}

inline uint64_t
EncodeHuffmanCode(const HuffmanCode code)
{
    /* Shift the code word to the most significant bit */
    uint64_t serialized = code.code_word << (64 - code.code_length);
    /* Encode the codelength in the least significant byte */
    serialized |= code.code_length;
    return serialized;
}

inline uint64_t
CreateStartEncodedHuffmanCode(const bool bit)
{
    return (uint64_t{bit} << cmc::bits::kBitIndexStart) + uint64_t{1};
}

inline void
CreateNextEncodedHuffmanCode(uint64_t& code, const bool bit)
{
    /* Append the next bit */
    code |= (uint64_t{bit} << (cmc::bits::kBitIndexStart - (code & uint64_t{0x00000000000000FF})));
    /* Increment the code_length*/
    ++code;
}

template <typename T>
struct INode
{
public:
    virtual ~INode() {};

    const T frequency;

protected:
    INode(const T frequency)
    : frequency(frequency) {};
};

template <typename T>
struct InternalNode : public INode<FrequencyType>
{
public:
    InternalNode(INode<FrequencyType>* new_child1, INode<FrequencyType>* new_child2)
    : INode<FrequencyType>(new_child1->frequency + new_child2->frequency), left(new_child1), right(new_child2) {};

    ~InternalNode()
    {
        if (left != nullptr)
        {
            delete left;
        }
        if (right != nullptr)
        {
            delete right;
        }
    };

    INode<FrequencyType>* const left{nullptr};
    INode<FrequencyType>* const right{nullptr};
};


template <typename T>
class LeafNode : public INode<FrequencyType>
{
public:
    LeafNode(const FrequencyType frequency, const T symbol)
    : INode<FrequencyType>(frequency), symbol(symbol) {};

    const T symbol;
};

template <typename T>
struct NodeCompare
{
    bool operator()(const INode<T>* lhs, const INode<T>* rhs) const {return lhs->frequency > rhs->frequency;}
};

template<typename T>
class HuffmanCoder 
{
public:
    HuffmanCoder() = delete;
    HuffmanCoder(const std::vector<EntropySymbol<T>>& symbols_and_frequencies)
    : tree_(symbols_and_frequencies)
    {
        if (symbols_and_frequencies.empty()) [[unlikely]]
        {
            cmc_err_msg("The symbol frequency table for the HuffmanCoder is empty.");
        }

        /* Get the codes from the Huffman tree */
        codes_ = tree_.GetHuffmanCodes();
    }

    ~HuffmanCoder() = default;

    HuffmanCode EncodeSymbol(const T symbol) const;
    std::vector<uint8_t> SerializeHuffmanCodes() const;
    std::vector<uint64_t> SerializeHuffmanCodesBEPadded() const;

private:
    HuffmanTree<T> tree_;
    HuffmanCodeMap<T> codes_;
};

template <typename T>
inline 
HuffmanCode
HuffmanCoder<T>::EncodeSymbol(const T symbol) const
{
    /* Find the code for the given symbol in the generated Huffman codes */
    auto code = codes_.find(symbol);

    cmc_assert(code != codes_.end());

    if (code == codes_.end()) [[unlikely]]
    {
        cmc_err_msg("The symbol ", symbol, " is not in the symbol-frequency-table of the HuffmanCoder.");
    }

    /* Return the bits as well as the length of the code */
    return code->second;
}

template <typename T>
inline std::vector<uint8_t>
HuffmanCoder<T>::SerializeHuffmanCodes() const
{
    std::vector<uint8_t> serialized_codes;
    serialized_codes.reserve(codes_.size() * (sizeof(T) + sizeof(HuffmanCodeWord) + 2 * sizeof(HuffmanCodeInfoType)));

    /* Push back the number of symbols/codes */
    cmc::bits::SerializeBEToByteStream(serialized_codes, static_cast<HuffmanCodeInfoType>(codes_.size()));

    /* Push back the type of the symbol */
    cmc::bits::SerializeBEToByteStream(serialized_codes, static_cast<HuffmanCodeInfoType>(ConvertToCmcType<T>()));

    /* Iterate through the code book and serialize the values */
    for (const auto&[symbol, huffcode] : codes_)
    {
        /* Store the codeword first */
        cmc::bits::SerializeBEToByteStream<HuffmanCodeWord>(serialized_codes, EncodeHuffmanCode(huffcode));

        /* Store the symbol afetrwards */
        cmc::bits::SerializeBEToByteStream<T>(serialized_codes, symbol);
    }

    return serialized_codes;
}


template <typename T>
inline std::vector<uint64_t>
HuffmanCoder<T>::SerializeHuffmanCodesBEPadded() const
{
    static_assert(2 * sizeof(HuffmanCodeInfoType) == sizeof(uint64_t));

    cmc::bits::vector serialized_codes;
    serialized_codes.Reserve(((codes_.size() + 1) * 2 * sizeof(HuffmanCodeInfoType) * cmc::bits::kCharBit));

    /* Push back the number of symbols/codes */
    serialized_codes.AppendBits(static_cast<HuffmanCodeInfoType>(codes_.size()), 0, 0);

    /* Push back the type of the symbol */
    serialized_codes.AppendBits(static_cast<HuffmanCodeInfoType>(ConvertToCmcType<T>()), 0, 0);

    /* Iterate through the code book and serialize the values */
    for (const auto&[symbol, huffcode] : codes_)
    {
        /* Store the codeword first */
        serialized_codes.AppendBits(static_cast<HuffmanCodeWord>(EncodeHuffmanCode(huffcode)), 0, 0);

        /* Store the symbol afetrwards */
        serialized_codes.AppendBits(static_cast<T>(symbol), 0, 0);
    }

    return serialized_codes.GetSerializedByteStreamBE();
}

template<typename T>
class HuffmanTree
{
public:
    HuffmanTree() = delete;
    HuffmanTree(const std::vector<EntropySymbol<T>> symbols_and_frequencies)
    {
        if (symbols_and_frequencies.empty())
        {
            cmc_err_msg("The symbol frequency table for the HuffmanCoder is empty.");
        }
        /* Construct the Huffman tree */
        this->ConstructTree(symbols_and_frequencies);
    }

    ~HuffmanTree() {
        if (root_ != nullptr)
        {
            delete root_;
        }
    };

    HuffmanCodeMap<T> GetHuffmanCodes() const;
private:
    void ConstructTree(const std::vector<EntropySymbol<T>>& symbols_and_frequencies);
    void GenerateCodes(const INode<FrequencyType>* node, const HuffmanCode& prefix, HuffmanCodeMap<T>& codes) const;

    INode<FrequencyType>* root_{nullptr};
};


template<typename T>
void HuffmanTree<T>::ConstructTree(const std::vector<EntropySymbol<T>>& symbols_and_frequencies)
{
    cmc_assert(symbols_and_frequencies.size() > static_cast<size_t>(1));

    std::priority_queue<INode<FrequencyType>*, std::vector<INode<FrequencyType>*>, NodeCompare<FrequencyType>> nodes;

    /* Create the leaf nodes from the non-zero frequency symbols */
    for (auto sym_freq_iter = symbols_and_frequencies.begin(); sym_freq_iter != symbols_and_frequencies.end(); ++sym_freq_iter)
    {
        if(sym_freq_iter->frequency > 0)
        {
            nodes.push(new LeafNode<T>(sym_freq_iter->frequency, sym_freq_iter->symbol));
        }
    }

    /* Construct the structure and root element of a Huffman tree */
    while (nodes.size() > 1)
    {
        /* Get the leaf or internal node with the lowest frequency */
        INode<FrequencyType>* lowest_frequency_node = nodes.top();
        nodes.pop();

        /* Get the next leaf or internal node with the lowest frequency */
        INode<FrequencyType>* second_lowest_frequency_node = nodes.top();
        nodes.pop();

        /* Create a new internal node with the two lowest frequency nodes */
        INode<FrequencyType>* parent = new InternalNode<T>(lowest_frequency_node, second_lowest_frequency_node);
        nodes.push(parent);
    }

    /* Store the root of the tree */
    root_ = nodes.top();
}

template<typename T>
void HuffmanTree<T>::GenerateCodes(const INode<FrequencyType>* node, const HuffmanCode& prefix, HuffmanCodeMap<T>& codes) const
{
    if (const LeafNode<T>* leaf = dynamic_cast<const LeafNode<T>*>(node))
    {
        /* If it is a leaf node, we store the prefix code for the given symbol in the map */
        codes[leaf->symbol] = prefix;
    }
    else if (const InternalNode<T>* internal_node = dynamic_cast<const InternalNode<T>*>(node))
    {
        /* If it is an internal node, we append a bit to the prefix code and recurse into the left and right child */
        HuffmanCode left_prefix = prefix;
        AppendUnsetBit(left_prefix);
        GenerateCodes(internal_node->left, left_prefix, codes);

        HuffmanCode right_prefix = prefix;
        AppendSetBit(right_prefix);
        GenerateCodes(internal_node->right, right_prefix, codes);
    }
}

template<typename T>
HuffmanCodeMap<T>
HuffmanTree<T>::GetHuffmanCodes() const
{
    std::unordered_map<T, HuffmanCode> codes;

    if (root_ != nullptr)
    {
        this->GenerateCodes(root_, HuffmanCode{}, codes);
    }

    return codes;
}

}

#endif /* !CMC_HUFFMAN_CODER_HXX */
