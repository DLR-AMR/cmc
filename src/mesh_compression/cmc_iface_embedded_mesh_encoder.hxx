#ifndef CMC_IFACE_EMBEDDED_MESH_ENCODER_HXX
#define CMC_IFACE_EMBEDDED_MESH_ENCODER_HXX

#include "mesh_compression/cmc_iface_abstract_mesh_encoder.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <vector>

namespace cmc::mesh_compression
{

class IEmbeddedMeshEncoder : public IAbstractMeshEncoder
{
public:

    virtual std::vector<uint8_t> EncodeRootLevelMesh(AmrMesh root_level_mesh, const GeoDomain& global_domain) = 0;    

    virtual ~IEmbeddedMeshEncoder(){};

protected:
    IEmbeddedMeshEncoder() = default;
        
private:

};

}


#endif /* !CMC_IFACE_EMBEDDED_MESH_ENCODER_HXX */
