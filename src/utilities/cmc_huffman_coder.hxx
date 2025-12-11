#ifndef CMC_HUFFMAN_CODER_HXX
#define CMC_HUFFMAN_CODER_HXX

#include "cmc.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_serialization.hxx"

#include <vector>
#include <cmath>
#include <map>
#include <queue>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <algorithm>

namespace cmc::entropy_coding::huffman
{

template <typename T>
class HuffmanTree;

using FrequencyType = uint32_t;
using CodeLengthType = uint8_t;
using SymbolInfoType = int32_t;
using HuffmanCodeInfoType = uint32_t;

using HuffmanCode = bit_vector::BitVector;
template <typename T>
using HuffmanCodeMap = std::map<T, HuffmanCode>;

constexpr bool kLeftBranch = false;
constexpr bool kRightBranch = true;

constexpr cmc::Endian kSerializationEndianness = Endian::Big;

template <typename T>
struct HuffmanSymbol
{
    HuffmanSymbol(const T symbol_, const FrequencyType frequency_)
    : symbol{symbol_}, frequency{frequency_} {}

    T symbol;
    FrequencyType frequency;
};

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
    bool operator()(const INode<T>* lhs, const INode<T>* rhs) const { return lhs->frequency > rhs->frequency; }
};

template<typename T>
class HuffmanCoder 
{
public:
    HuffmanCoder() = delete;
    HuffmanCoder(const std::vector<HuffmanSymbol<T>>& symbols_and_frequencies)
    : tree_(symbols_and_frequencies)
    {
        if (symbols_and_frequencies.empty())
        {
            cmc_err_msg("The symbol frequency table for the HuffmanCoder is empty.");
        }

        /* Get the codes from the Huffman tree */
        codes_ = tree_.GetHuffmanCodes();

        /* Serialize the Huffman Frequency Table */
        serialized_sym_freq_table_ = this->SerializeSymbolFrequencyTable(symbols_and_frequencies);
    }

    ~HuffmanCoder() = default;

    std::pair<std::vector<uint8_t>, size_t> EncodeSymbol(const T symbol) const;
    std::vector<uint8_t> GetSerializedSymbolFrequencyTable() const {return serialized_sym_freq_table_;}

private:
    std::vector<uint8_t> SerializeSymbolFrequencyTable(const std::vector<HuffmanSymbol<T>>& symbols_and_frequencies) const;

    HuffmanTree<T> tree_;
    HuffmanCodeMap<T> codes_;
    std::vector<uint8_t> serialized_sym_freq_table_;
};

template <typename T>
inline 
std::pair<std::vector<uint8_t>, size_t>
HuffmanCoder<T>::EncodeSymbol(const T symbol) const
{
    /* Find the code for the given symbol in the generated Huffman codes */
    auto code = codes_.find(symbol);

    cmc_assert(code != codes_.end());

    if (code == codes_.end())
    {
        cmc_err_msg("The symbol ", symbol, " is not in the symbol-frequency-table of the HuffmanCoder.");
    }

    /* Return the bits as well as the length of the code */
    return code->second.GetBits();
}

template <typename T>
std::vector<uint8_t>
HuffmanCoder<T>::SerializeSymbolFrequencyTable(const std::vector<HuffmanSymbol<T>>& symbols_and_frequencies) const
{
    std::vector<uint8_t> serialized_symbol_frequency_table;
    serialized_symbol_frequency_table.reserve(symbols_and_frequencies.size() * sizeof(T) + symbols_and_frequencies.size() * sizeof(FrequencyType) + sizeof(HuffmanCodeInfoType) + sizeof(SymbolInfoType));

    HuffmanCodeInfoType num_symbols = static_cast<HuffmanCodeInfoType>(symbols_and_frequencies.size());

    /* Count out zero frequencies */
    for (const auto&[_, frequency] : symbols_and_frequencies)
    {
        if (frequency <= 0)
        {
            --num_symbols;
        }
    }

    /* Push back the count of symbols */
    PushBackValueToByteStream<HuffmanCodeInfoType>(serialized_symbol_frequency_table, num_symbols, kSerializationEndianness);
    
    /* Push back the type of the symbol */
    PushBackValueToByteStream<SymbolInfoType>(serialized_symbol_frequency_table, static_cast<SymbolInfoType>(ConvertToCmcType<T>()), kSerializationEndianness);
    
    /* Iterate over the symbol frequency table */
    for (const auto&[symbol, frequency] : symbols_and_frequencies)
    {
        if (frequency > 0)
        {
            /* Push back the symbol */
            PushBackValueToByteStream<T>(serialized_symbol_frequency_table, static_cast<T>(symbol), kSerializationEndianness);
            
            /* Push back the frequency count */
            PushBackValueToByteStream<FrequencyType>(serialized_symbol_frequency_table, static_cast<FrequencyType>(frequency), kSerializationEndianness);
        }
    }

    return serialized_symbol_frequency_table;
}

template<typename T>
class HuffmanTree
{
public:
    HuffmanTree() = delete;
    HuffmanTree(const std::vector<HuffmanSymbol<T>> symbols_and_frequencies)
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

    T GetNextSymbol(bit_vector::BitVectorView& view) const;

private:
    void ConstructTree(const std::vector<HuffmanSymbol<T>>& symbols_and_frequencies);
    void GenerateCodes(const INode<FrequencyType>* node, const HuffmanCode& prefix, HuffmanCodeMap<T>& codes) const;

    INode<FrequencyType>* root_{nullptr};
};

template<typename T>
T
HuffmanTree<T>::GetNextSymbol(bit_vector::BitVectorView& view) const
{
    /* Start at the root element */
    INode<FrequencyType>* node = root_;

    /* Iterate until we will find a leaf element */
    while (const InternalNode<T>* current_node = dynamic_cast<const InternalNode<T>*>(node))
    {
        const bool flag = view.IsCurrentBitSet();
        view.MoveToNextBit();

        if (flag == kLeftBranch)
        {
            node = current_node->left;
        } else
        {
            cmc_assert(flag == kRightBranch);
            node = current_node->right;
        }
    }

    /* If a leaf element is reached, we will get the symbol from it and return it */
    const LeafNode<T>* leaf = dynamic_cast<const LeafNode<T>*>(node);

    cmc_assert(leaf != nullptr);
    
    return leaf->symbol;
}

template<typename T>
void HuffmanTree<T>::ConstructTree(const std::vector<HuffmanSymbol<T>>& symbols_and_frequencies)
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
        left_prefix.AppendBit(kLeftBranch);
        GenerateCodes(internal_node->left, left_prefix, codes);
        HuffmanCode right_prefix = prefix;
        right_prefix.AppendBit(kRightBranch);
        GenerateCodes(internal_node->right, right_prefix, codes);
    }
}

template<typename T>
HuffmanCodeMap<T>
HuffmanTree<T>::GetHuffmanCodes() const
{
    std::map<T, HuffmanCode> codes;

    if (root_ != nullptr)
    {
        this->GenerateCodes(root_, bit_vector::BitVector(), codes);
    }

    return codes;
}

template<typename T>
class HuffmanDecoder 
{
public:
    HuffmanDecoder() = delete;
    HuffmanDecoder(const uint8_t* start_encoding_pos);

    ~HuffmanDecoder()
    {
        if (tree_ != nullptr)
        {
            delete tree_;
        }
    };

    void StartDecoding(const bit_vector::BitVectorView& encoding);
    T DecodeNextSymbol() {cmc_assert(tree_ != nullptr); return tree_->GetNextSymbol(encoded_stream_view_);}
    size_t GetNumberOfProcessedBytesForSymbolFrequencyTable() const {return num_processed_bytes_sym_freq_table_;}
private:
    std::pair<std::vector<HuffmanSymbol<T>>, size_t> ReconstructHuffmanSymbolFrequencyTable(const uint8_t* start_encoding_pos);

    HuffmanTree<T>* tree_{nullptr};
    size_t num_processed_bytes_sym_freq_table_{0};
    bit_vector::BitVectorView encoded_stream_view_;
};

template <typename T>
HuffmanDecoder<T>::HuffmanDecoder(const uint8_t* start_encoding_pos)
{
    /* Reconstruct the symbol frequency table */
    const auto [symbol_frequency_table, num_processed_bytes] = this->ReconstructHuffmanSymbolFrequencyTable(start_encoding_pos);
    
    /* Construct the tree */
    tree_ = new HuffmanTree<T>(symbol_frequency_table);

    /* Store the processed bytes for the decoding of the symbol frequency table */
    num_processed_bytes_sym_freq_table_ = num_processed_bytes;
}

template <typename T>
void
HuffmanDecoder<T>::StartDecoding(const bit_vector::BitVectorView& encoding)
{
    encoded_stream_view_ = encoding;
}

template <typename T>
std::pair<std::vector<HuffmanSymbol<T>>, size_t>
HuffmanDecoder<T>::ReconstructHuffmanSymbolFrequencyTable(const uint8_t* start_encoding_pos)
{
    size_t offset{0};

    const HuffmanCodeInfoType num_symbols = GetValueFromByteStream<HuffmanCodeInfoType>(start_encoding_pos, kSerializationEndianness);
    offset += sizeof(HuffmanCodeInfoType);

    const SymbolInfoType data_type = GetValueFromByteStream<SymbolInfoType>(start_encoding_pos + offset, kSerializationEndianness);
    offset += sizeof(SymbolInfoType);

    if (static_cast<SymbolInfoType>(ConvertToCmcType<T>()) != data_type)
    {
        cmc_err_msg("The template parameter does not coincide with the symbol type of the Huffman symbol frequency table.");
    }

    std::vector<HuffmanSymbol<T>> symbol_frequency_table;
    symbol_frequency_table.reserve(num_symbols);

    /* Iterate until the symbol frequency table has been re-created */
    for (HuffmanCodeInfoType iter{0}; iter < num_symbols; ++iter)
    {
        const T deserialized_symbol = GetValueFromByteStream<T>(start_encoding_pos + offset, kSerializationEndianness);
        offset += sizeof(T);
    
        const FrequencyType deserialized_freq = GetValueFromByteStream<FrequencyType>(start_encoding_pos + offset, kSerializationEndianness);
        offset += sizeof(FrequencyType);
        
        symbol_frequency_table.emplace_back(deserialized_symbol, deserialized_freq);
    }

    return std::make_pair(symbol_frequency_table, offset);
}


}

#endif /* !CMC_HUFFMAN_CODER_HXX */
