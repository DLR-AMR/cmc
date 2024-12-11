#ifndef CMC_HUFFMAN_HXX
#define CMC_HUFFMAN_HXX

#include "utilities/cmc_bit_map.hxx"

#include <iostream>
#include <queue>
#include <map>
#include <climits>
#include <iterator>
#include <algorithm>
#include <limits>

namespace cmc
{

namespace huffman
{

typedef bit_vector::BitVector HuffmanCode;
typedef std::map<uint8_t, HuffmanCode> HuffmanCodeMap;


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
struct InternalNode : public INode<T>
{
public:
    InternalNode(INode<T>* new_child1, INode<T>* new_child2)
    : INode<T>(new_child1->frequency + new_child2->frequency), left(new_child1), right(new_child2) {};

    ~InternalNode()
    {
        delete left;
        delete right;
    };

    INode<T>* const left;
    INode<T>* const right;
};

template <typename T>
class LeafNode : public INode<T>
{
public:
    LeafNode(const T frequency, const uint8_t symbol)
    : INode<T>(frequency), symbol(symbol) {};

    const uint8_t symbol;
};

template <typename T>
struct NodeCompare
{
    bool operator()(const INode<T>* lhs, const INode<T>* rhs) const { return lhs->frequency > rhs->frequency; }
};

template <typename FrequenceType>
class HuffmanTree
{
public:
    HuffmanTree() = delete;
    HuffmanTree(const std::vector<FrequenceType>& value_frequency_table)
    {
        this->ConstructTree(value_frequency_table);
        this->GenerateCodes(root_, bit_vector::BitVector(), this->codes_);
    };

    ~HuffmanTree() {
        delete root_;
    };

    std::pair<std::vector<uint8_t>, size_t> EncodeSymbol(const uint8_t symbol) const;
    void EncodeSymbols(const std::vector<uint8_t>& symbols) const;

    const HuffmanCodeMap& GetHuffmanCodes() const {return codes_;}

    bit_vector::BitVector SerializeHuffmanTree(const size_t symbol_bit_size) const;

private:
    void ConstructTree(const std::vector<FrequenceType>& value_frequency_table);
    void GenerateCodes(const INode<FrequenceType>* node, const HuffmanCode& initial_prefix, HuffmanCodeMap& codes) const;

    std::map<uint8_t, HuffmanCode> codes_;
    INode<FrequenceType>* root_;
};

template <typename FrequenceType>
inline void
HuffmanTree<FrequenceType>::ConstructTree(const std::vector<FrequenceType>& value_frequency_table)
{
    std::priority_queue<INode<FrequenceType>*, std::vector<INode<FrequenceType>*>, NodeCompare<FrequenceType>> nodes;

    /* Create the leaf nodes from the non-zero frequency symbols */
    for (auto freq_iter = value_frequency_table.begin(); freq_iter != value_frequency_table.end(); ++freq_iter)
    {
        if(*freq_iter != 0)
        {
            cmc_assert(std::distance(value_frequency_table.begin(), freq_iter) >= 0 && std::distance(value_frequency_table.begin(), freq_iter) <= std::numeric_limits<uint8_t>::max());
            nodes.push(new LeafNode<FrequenceType>(*freq_iter, static_cast<uint8_t>(std::distance(value_frequency_table.begin(), freq_iter))));
        }
    }

    /* Construct the structure and root element of a Huffman tree */
    while (nodes.size() > 1)
    {
        /* Get the leaf or internal node with the lowest frequency */
        INode<FrequenceType>* lowest_frequency_node = nodes.top();
        nodes.pop();

        /* Get the next leaf or internal node with the lowest frequency */
        INode<FrequenceType>* second_lowest_frequency_node = nodes.top();
        nodes.pop();

        /* Create a new internal node with the two lowest frequency nodes */
        INode<FrequenceType>* parent = new InternalNode(lowest_frequency_node, second_lowest_frequency_node);
        nodes.push(parent);
    }

    /* Store the root of the tree */
    root_ = nodes.top();
}

template <typename FrequenceType>
inline void
HuffmanTree<FrequenceType>::GenerateCodes(const INode<FrequenceType>* node, const HuffmanCode& prefix, HuffmanCodeMap& codes) const
{
    if (const LeafNode<FrequenceType>* leaf = dynamic_cast<const LeafNode<FrequenceType>*>(node))
    {
        /* If it is a leaf node, we store the prefix code for the given symbol in the map */
        codes[leaf->symbol] = prefix;
    }
    else if (const InternalNode<FrequenceType>* internal_node = dynamic_cast<const InternalNode<FrequenceType>*>(node))
    {
        /* If it is an internal node, we append a bit to the prefix code and recurse into the left and right child */
        HuffmanCode left_prefix = prefix;
        left_prefix.AppendBit(false);
        GenerateCodes(internal_node->left, left_prefix, codes);
        HuffmanCode right_prefix = prefix;
        right_prefix.AppendBit(true);
        GenerateCodes(internal_node->right, right_prefix, codes);
    }
}

template <typename FrequenceType>
inline bit_vector::BitVector
HuffmanTree<FrequenceType>::SerializeHuffmanTree(const size_t symbol_bit_size) const
{
    cmc_assert(symbol_bit_size <= bit_vector::kCharBit && symbol_bit_size != 0);

    bit_vector::BitVector serialized_map;
    serialized_map.Reserve(codes_.size() * symbol_bit_size);

    const int symbol_size = static_cast<int>(symbol_bit_size);
    //const int bit_shift = bit_vector::kCharBit - symbol_bit_size;

    /* In case the codebook consists of only a single code, the corresponding code would be of length zero.
     * Therefore, we only encode the symbol */
    if (codes_.size() == 1)
    {
        const uint8_t& symbol = codes_.begin()->first;
        serialized_map.AppendBits(symbol, symbol_size);
        return serialized_map;
    }

    cmc_assert(codes_.size() >= 2);

    /* Iterate over all codes and serialize them in the order [symbol; num_bits_code; code] */
    for (auto code_iter = codes_.begin(); code_iter != codes_.end(); ++code_iter)
    {
        //cmc_debug_msg("Huffman Code Bit size: ", code_iter->second.size_bits());
        /* Get the symbol */
        const uint8_t& symbol = code_iter->first;

        /* Store the trimmed (bit-wise) symbol in the BitVector */
        //serialized_map.AppendBits(symbol << bit_shift, symbol_size);//Does not need to be shifted
        //cmc_debug_msg("The symbol will be appended, num bits: ", symbol_size);
        serialized_map.AppendBits(symbol, symbol_size);

        /* Get the code as well as the length of the code */
        const auto [bits, num_bits] = code_iter->second.GetBits();

        //cmc_debug_msg("Serialize huffman tree: num_bits: ", num_bits);

        cmc_assert(num_bits < (1 << symbol_bit_size));

        /* Set the length of the code */
        //serialized_map.AppendBits(static_cast<uint8_t>(num_bits) << bit_shift, symbol_size);//Same here: no shifting required
        //cmc_debug_msg("The code length will be appended: bit_length: ", symbol_size);
        serialized_map.AppendBits(static_cast<uint8_t>(num_bits), symbol_size);

        //cmc_debug_msg("Huffman Code will be appened: num_bits: ", num_bits);
        /* Set the code itself */
        serialized_map.AppendBits(bits, static_cast<int>(num_bits));

        //cmc_debug_msg("Huffman code has been serialized");
    }

    return serialized_map;
}

#if 0
inline void
HuffmanTree::EncodeSymbols(const std::vector<uint8_t>& symbols) const
{

    //continue...

    for (auto it = codes_.begin(); it != codes_.end(); ++it)
    {
        std::cout << static_cast<int>(it->first) << " encoded with num bits: " << it->second.size_bits() << std::endl;

        auto [bits, sizebits] = it->second.GetBits();
        //cmc_debug_msg("Code hat ", sizebits, " bits und sieht aus: ", std::bitset<8>(bits.front()));

        //std::cout << static_cast<int>(it->first) << " encoded as: ";
        //std::copy(it->second.begin(), it->second.end(),
        //          std::ostream_iterator<bool>(std::cout));
        //std::cout << std::endl;
    }


    //Fuer symbol size 5
    // 10 * (5 + 5) + 2+3+3+2+3+4+6+5+7+7 = 100 + 42 = 142
    // In the canonical variation it should be 10 * (5+5) = 100
    //5 is sufficient for four byte types; for double 6 would be needed
    const size_t num_bits_needed_to_store_symbol = 5;
    bit_vector::BitVector serialized_huffman = SerializeHuffmanTree(num_bits_needed_to_store_symbol);

    cmc_debug_msg("Huffman Tree has been serialized in num bits: ", serialized_huffman.size_bits());



}
#endif


template <typename FrequenceType>
inline 
std::pair<std::vector<uint8_t>, size_t>
HuffmanTree<FrequenceType>::EncodeSymbol(const uint8_t symbol) const
{
    /* Find the code for the given symbol in the generated Huffman codes */
    auto code = codes_.find(symbol);

    cmc_assert(code != codes_.end());

    /* Return the bits as well as the length of the code */
    return code->second.GetBits();
}

}
}

#endif /* !CMC_HUFFMAN_HXX */
