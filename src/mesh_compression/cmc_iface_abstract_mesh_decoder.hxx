#ifndef CMC_IFACE_ABSTRACT_MESH_DECODER_HXX
#define CMC_IFACE_ABSTRACT_MESH_DECODER_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_serialization.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <vector>
#include <utility>

namespace cmc::mesh_compression
{

class IAbstractMeshDecoder
{
public:

    void IntializeDecompressionIteration();
    void FinalizeDecompressionIteration();
    bool WillNextElementBeRefined();

    bool IsDecompressionProgressing();

    int GetNumberOfDecompressionIterations() const {return decompression_iteration_count_;};
    t8_gloidx_t GetCurrentNumberOfGlobalElements() const {return current_num_global_elements_;};

    virtual ~IAbstractMeshDecoder(){};

protected:
    IAbstractMeshDecoder() = delete;
    IAbstractMeshDecoder(const std::vector<uint8_t>& encoded_mesh_byte_stream)
    : encoded_mesh_(encoded_mesh_byte_stream) {};

    void AddToProcessedBytes(const size_t num_bytes_to_add) {processed_bytes_ += num_bytes_to_add;};
    const std::vector<uint8_t>& GetEncodedMeshStreamData() const {return encoded_mesh_;};
    
private:
    const std::vector<uint8_t>& encoded_mesh_;
    bit_map::BitMapView level_refinement_indications_;
    t8_gloidx_t current_num_global_elements_{0};
    size_t processed_bytes_{0};
    int decompression_iteration_count_{0};
    
};

inline bool
IAbstractMeshDecoder::IsDecompressionProgressing()
{
    return (processed_bytes_ < encoded_mesh_.size());
}

inline void
IAbstractMeshDecoder::IntializeDecompressionIteration()
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
IAbstractMeshDecoder::FinalizeDecompressionIteration()
{
    ++decompression_iteration_count_;
}

inline bool
IAbstractMeshDecoder::WillNextElementBeRefined()
{
    return level_refinement_indications_.GetNextBit();
}

}


#endif /* !CMC_IFACE_ABSTRACT_MESH_DECODER_HXX */
