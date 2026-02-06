#ifndef CMC_IFACE_MESH_ENCODER_HXX
#define CMC_IFACE_MESH_ENCODER_HXX

#include "mesh_compression/cmc_iface_abstract_mesh_encoder.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#endif

namespace cmc::mesh_compression
{

class IMeshEncoder : public IAbstractMeshEncoder
{
public:

    virtual std::vector<uint8_t> EncodeRootLevelMesh(t8_forest_t root_level) = 0;
#ifdef CMC_ENABLE_MPI
    virtual std::vector<uint8_t> GatherRootLevelMeshEncodingOnTheRootRank(const t8_forest_t root_level, const MPI_Comm comm) = 0;
#endif

    virtual ~IMeshEncoder(){};

protected:
    IMeshEncoder() = default;
        
private:

};

}


#endif /* !CMC_IFACE_MESH_ENCODER_HXX */
