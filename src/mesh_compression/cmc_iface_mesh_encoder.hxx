#ifndef CMC_IFACE_MESH_ENCODER_HXX
#define CMC_IFACE_MESH_ENCODER_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_serialization.hxx"
#include "t8code/cmc_t8_mesh.hxx"

namespace cmc::mesh_compression
{

class IMeshEncoder
{
public:

    void IntializeCompressionIteration(const t8_locidx_t num_local_elems);
    void FinalizeCompressionIteration();
    void IndicateCoarsening();
    void IndicateElementStaysUnchanged();
    int GetNumberOfCompressionIterations() const {return compression_iteration_count_;};

    std::vector<uint8_t> GetEncodedLevelData();

    virtual std::vector<uint8_t> EncodeRootLevelMesh(t8_forest_t root_level) = 0;

    virtual ~IMeshEncoder(){};

protected:
    IMeshEncoder() = default;
        
private:
    bit_map::BitMap encoded_level_data_;
    int compression_iteration_count_{0};
};

inline void
IMeshEncoder::IntializeCompressionIteration(const t8_locidx_t num_local_elems)
{
    encoded_level_data_ = bit_map::BitMap();
    encoded_level_data_.Reserve(num_local_elems);
}

inline void
IMeshEncoder::IndicateCoarsening()
{
    encoded_level_data_.AppendSetBit();
}

inline void
IMeshEncoder::IndicateElementStaysUnchanged()
{
    encoded_level_data_.AppendUnsetBit();
}

inline void
IMeshEncoder::FinalizeCompressionIteration()
{
    ++compression_iteration_count_;
}

inline std::vector<uint8_t>
IMeshEncoder::GetEncodedLevelData()
{
    /* Get the encoded data from this level and return it */
    std::vector<uint8_t> encoded_stream;
    encoded_stream.reserve(encoded_level_data_.size_bytes() + sizeof(size_t));

    const size_t encoded_bits = encoded_level_data_.size();
    PushBackValueToByteStream(encoded_stream, encoded_bits);

    std::copy_n(encoded_level_data_.begin_bytes(), encoded_level_data_.size_bytes(), std::back_insert_iterator(encoded_stream));

    return encoded_stream;
}


}


#endif /* !CMC_IFACE_MESH_ENCODER_HXX */
