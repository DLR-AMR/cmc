#ifndef CMC_PATCH_AMR_COMPRESSION_HXX
#define CMC_PATCH_AMR_COMPRESSION_HXX

#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_compression_settings.hxx"

#include <vector>
#include <variant>

namespace cmc
{
    
namespace patch_amr
{

template <typename T>
class CompressionBlock
{
public:
    CompressionBlock() = default;
    CompressionBlock(const int block_id, const Hyperslab& domain, const DataLayout layout)
    : block_id_{block_id}, domain_{domain}, layout_{layout} {};
    CompressionBlock(const int block_id, Hyperslab&& domain, const DataLayout layout)
    : block_id_{block_id}, domain_{std::move(domain)}, layout_{layout} {};

    void OrderMortonCurveCompliant();
    int GetDimensionality() const {return GetDimensionalityOfDataLayout(layout_);}
    const std::vector<T>& GetInitialData() const {return data_;};
    bool IsMortonOrdered() const {return is_morton_ordered_;}

    void SetData(std::vector<T>&& data) {data_ = std::move(data);}
    int GetBlockID() const {return block_id_;};

    size_t size() const {return data_.size();};
    const Hyperslab& GetHyperslab() const {return domain_;};

    void ClearData() {data_.clear(); invalidated_ = true;};
private:
    int block_id_;
    Hyperslab domain_;
    DataLayout layout_;
    std::vector<T> data_;
    bool is_morton_ordered_{false};
    bool invalidated_{false};
};

template <typename T>
void
CompressionBlock<T>::OrderMortonCurveCompliant()
{
    if (is_morton_ordered_) {return;}

    //cmc_debug_msg("Vor reordering sieht data so aus:");
    //for (auto diter = data_.begin(); diter != data_.end(); ++diter)
    //{
    //    cmc_debug_msg(*diter);
    //}
    //cmc_debug_msg("\n\n\n");

    /* Create a hyperslab with the same count but without a start offset */
    Hyperslab resetted_hs = domain_;
    resetted_hs.NullifyStartIndices();

    /* Get the corresponding hyperslab iteration function */
    const UpdateHyperslabCoordinateFn HsCoordUpdateFn = GetHyperslabCoordinatesIterationFunction(layout_);

    const int dimensionality = GetDimensionalityOfDataLayout(layout_);

    /* Create a vector of the given dimensionality to hold the refernce coordinates which can then be transformed to Morton indices */
    std::vector<HyperslabIndex> reference_coordinates(dimensionality);

    /* Allocate a vector holding the ordered data */
    std::vector<T> ordered_data(data_.size());

    /* Get the number of coordinates the hyperslab covers */
    const HyperslabIndex num_coords = resetted_hs.GetNumberCoordinates();

    /* Iterate over all Morton indices and sort the data accordingly */
    for (HyperslabIndex coord_iter = 0; coord_iter < num_coords; ++coord_iter)
    {
        /* Update the refernce coordinates for the current hyperlsab index */
        HsCoordUpdateFn(resetted_hs, reference_coordinates, coord_iter);

        /* Get the corresponding Morton index */
        const MortonIndex pos = GetMortonIndex(reference_coordinates, dimensionality);
        //cmc_debug_msg("Für coord_iter: ", coord_iter, " kommt morton index: ", pos);
        /* Set the data at the coorect position */
        ordered_data[pos] = data_[coord_iter];
    }

    /* Exchange the ordered with the unordered data */
    std::swap(ordered_data, data_);

    //cmc_debug_msg("Nach reordering sieht data so aus:");
    //for (auto diter = data_.begin(); diter != data_.end(); ++diter)
    //{
    //    cmc_debug_msg(*diter);
    //}
    /* Set the flag accordingly */
    is_morton_ordered_ = true;

    return;
}

//using GeneralCompressionBlock = std::variant<CompressionBlock<int8_t>, CompressionBlock<char>, CompressionBlock<int16_t>, CompressionBlock<int32_t>,
//                                             CompressionBlock<float>, CompressionBlock<double>, CompressionBlock<uint8_t>, CompressionBlock<uint16_t>,
//                                             CompressionBlock<uint32_t>, CompressionBlock<int64_t>, CompressionBlock<uint64_t>>;
//
//class Block
//{
//public:
//    GeneralCompressionBlock compression_block_;
//};


class Compressor
{
public:
    Compressor() = default;

    Compressor(const std::vector<InputVar>& variables_to_compress, const CompressionSettings& settings)
    : variables_{variables_to_compress}, compression_settings_{settings} {};
    Compressor(std::vector<InputVar>&& variables_to_compress, const CompressionSettings& settings)
    : variables_{std::move(variables_to_compress)}, compression_settings_{settings} {};
    Compressor(const std::vector<InputVar>& variables_to_compress, CompressionSettings&& settings)
    : variables_{variables_to_compress}, compression_settings_(std::move(settings)) {};
    Compressor(std::vector<InputVar>&& variables_to_compress, CompressionSettings&& settings)
    : variables_{std::move(variables_to_compress)}, compression_settings_(std::move(settings))  {};

    Compressor(const Compressor& other) = default;
    Compressor& operator=(const Compressor& other) = default;
    Compressor(Compressor&& other) = default;
    Compressor& operator=(Compressor&& other) = default;

    ~Compressor() = default;

    
    void Compress();
    //void WriteCompressedData(const std::string& file_name, const int time_step) const;

private:
    //New members
    //only input variables which are then blocked and comrpessed that way
    std::vector<InputVar> variables_;
    CompressionSettings compression_settings_;
};

}

}

#endif /* !CMC_PATCH_AMR_COMPRESSION_HXX */
