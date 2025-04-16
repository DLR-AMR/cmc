#ifndef CMC_IFACE_MESH_DECODER_HXX
#define CMC_IFACE_MESH_DECODER_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_serialization.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <vector>
#include <utility>

namespace cmc::mesh_compression
{

class IMeshDecoder
{
public:

    void IntializeDecompressionIteration();
    void FinalizeDecompressionIteration();
    bool WillNextElementBeRefined();

    bool IsDecompressionProgressing();

    int GetNumberOfDecompressionIterations() const {return decompression_iteration_count_;};
    t8_gloidx_t GetCurrentNumberOfGlobalElements() const {return current_num_global_elements_;};

    t8_forest_t DecodeRootMesh();

    virtual ~IMeshDecoder(){};

protected:
    IMeshDecoder() = delete;
    IMeshDecoder(const std::vector<uint8_t>& encoded_mesh_byte_stream)
    : encoded_mesh_(encoded_mesh_byte_stream) {};

    /* Return the decoded root level mesh as well as the processes bytes needed for decoding the mesh */
    virtual std::pair<t8_forest_t, size_t> DecodeRootLevelMesh(const std::vector<uint8_t>& encoded_mesh_stream_) = 0;

private:
    const std::vector<uint8_t>& encoded_mesh_;
    bit_map::BitMapView level_refinement_indications_;
    t8_gloidx_t current_num_global_elements_{0};
    size_t processed_bytes_{0};
    int decompression_iteration_count_{0};
    
};

inline bool
IMeshDecoder::IsDecompressionProgressing()
{
    return (processed_bytes_ < encoded_mesh_.size());
}

inline t8_forest_t
IMeshDecoder::DecodeRootMesh()
{
    /* Decode the root level mesh */
    auto [root_mesh, num_processed_bytes] = DecodeRootLevelMesh(encoded_mesh_);

    /* Store the offset for the already processed bytes of the root mesh */
    processed_bytes_ = num_processed_bytes;

    /* Return the decoded root mesh */
    return root_mesh;
}

inline void
IMeshDecoder::IntializeDecompressionIteration()
{
    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = encoded_mesh_.data();

    /* Get the number of elements on this level which is equal to the number of refinement indications */
    current_num_global_elements_ = static_cast<t8_gloidx_t>(GetValueFromByteStream<size_t>(data_start_ptr + processed_bytes_));
    
    cmc_debug_msg("The number of refinement indication (which is equal to the number of elements) is ", current_num_global_elements_, " during this decompression iteration.");

    /* Update the processes offset */
    processed_bytes_ += sizeof(size_t);

    /* Set the view on the current level refinement bytes */
    level_refinement_indications_ = bit_map::BitMapView(data_start_ptr + processed_bytes_, current_num_global_elements_);

    /* Update the processed bytes to the next level */
    processed_bytes_ += (current_num_global_elements_ / bit_map::kCharBit) + (current_num_global_elements_ % bit_map::kCharBit != 0 ? 1 : 0);
}

inline void
IMeshDecoder::FinalizeDecompressionIteration()
{
    ++decompression_iteration_count_;
}

inline bool
IMeshDecoder::WillNextElementBeRefined()
{
    return level_refinement_indications_.GetNextBit();
}

}


#endif /* !CMC_IFACE_MESH_DECODER_HXX */
