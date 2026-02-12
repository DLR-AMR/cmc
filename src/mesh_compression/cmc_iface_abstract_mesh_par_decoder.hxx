#ifndef CMC_IFACE_ABSTRACT_MESH_PAR_DECODER_HXX
#define CMC_IFACE_ABSTRACT_MESH_PAR_DECODER_HXX


#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_serialization.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <vector>
#include <tuple>

namespace cmc::mesh_compression
{

class IAbstractMeshParDecoder
{
public:

    void IntializeDecompressionIteration(const t8_gloidx_t global_elem_offset);
    void FinalizeDecompressionIteration();
    bool WillNextElementBeRefined();

    int GetNumberOfDecompressionIterations() const {return decompression_iteration_count_;};

    bit_map::BitMapView GetGlobalMeshEncodingStep(const int decompression_step) const;

    virtual ~IAbstractMeshParDecoder(){};

protected:
    IAbstractMeshParDecoder() = delete;
    IAbstractMeshParDecoder(const uint8_t* encoded_mesh_lvl_byte_stream)
    : global_encoded_mesh_(encoded_mesh_lvl_byte_stream) {}
    
    /* Return the decoded root level mesh as well as the processes bytes needed for decoding the mesh */
    std::pair<t8_forest_t, int> DecodeRootLevelMesh(const t8_cmesh_t cmesh, const t8_scheme *scheme);

private:
    /* Return the decoded root level mesh as well as the processes bytes needed for decoding the mesh */
    virtual std::tuple<t8_forest_t, size_t, int> DecodeRootLevel(const uint8_t* encoded_mesh_stream_, const t8_cmesh_t cmesh, const t8_scheme *scheme) = 0;

    const uint8_t* global_encoded_mesh_{nullptr};
    bit_map::BitMapView level_refinement_indications_;
    size_t offset_refinement_structure_encoding_{0};
    int decompression_iteration_count_{0};
};

inline std::pair<t8_forest_t, int>
IAbstractMeshParDecoder::DecodeRootLevelMesh(const t8_cmesh_t cmesh, const t8_scheme *scheme)
{
    /* Recreate the base mesh and get the processed bytes from the stream */
    auto [base_mesh, processed_bytes, dim] = this->DecodeRootLevel(global_encoded_mesh_, cmesh, scheme);

    /* Store the offset to the first encoded refinement structure */
    offset_refinement_structure_encoding_ = processed_bytes;

    ++decompression_iteration_count_;

    return std::make_pair(base_mesh, dim);
}

bit_map::BitMapView
IAbstractMeshParDecoder::GetGlobalMeshEncodingStep(const int decompression_step) const
{
    size_t offset = offset_refinement_structure_encoding_;

    /* Iterate until we find the correct offset for the specified level */
    for (int step_iter{0}; step_iter < decompression_step; ++step_iter)
    {
        /* Getthe number of global elements in this level */
        const uint64_t current_num_global_elements = GetValueFromByteStream<uint64_t>(global_encoded_mesh_ + offset);

        /* Add the global level bit count */
        offset += sizeof(uint64_t);

        /* Move to the first byte of the succeeding level */
        offset += (current_num_global_elements / bit_map::kCharBit) + (current_num_global_elements % bit_map::kCharBit != 0 ? 1 : 0);
    }

    /* Apply the offset of the current level */
    const uint64_t num_level_elements = GetValueFromByteStream<uint64_t>(global_encoded_mesh_ + offset);
    offset += sizeof(uint64_t);

    /* Create the BitMapView */
    return bit_map::BitMapView(global_encoded_mesh_ + offset, num_level_elements);
}

inline void
IAbstractMeshParDecoder::IntializeDecompressionIteration(const t8_gloidx_t global_elem_offset)
{
    /* Get the number of elements on this level which is equal to the number of refinement indications */
    const uint64_t current_num_global_elements = GetValueFromByteStream<uint64_t>(global_encoded_mesh_);
    
    cmc_debug_msg("The number of refinement indications (which is equal to the number of elements) is ", current_num_global_elements, " during this decompression iteration.");

    /* Set the view on the current level refinement bytes */
    level_refinement_indications_ = bit_map::BitMapView(global_encoded_mesh_ + offset_refinement_structure_encoding_ + sizeof(uint64_t), current_num_global_elements);

    cmc_assert(global_elem_offset <= current_num_global_elements);

    /* Move to the correct process-local corresponding offset */
    level_refinement_indications_.MoveToStartBit(global_elem_offset);
}

inline void
IAbstractMeshParDecoder::FinalizeDecompressionIteration()
{
    ++decompression_iteration_count_;
}

inline bool
IAbstractMeshParDecoder::WillNextElementBeRefined()
{
    return level_refinement_indications_.GetNextBit();
}

}

#endif /* !CMC_IFACE_ABSTRACT_MESH_PAR_DECODER_HXX */
