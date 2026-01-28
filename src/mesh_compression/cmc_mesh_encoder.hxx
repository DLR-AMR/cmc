#ifndef CMC_MESH_ENCODER_HXX
#define CMC_MESH_ENCODER_HXX

#include "mesh_compression/cmc_iface_mesh_encoder.hxx"
#include "utilities/cmc_serialization.hxx"

#include <vector>

namespace cmc::mesh_compression
{

class MeshEncoder : public IMeshEncoder
{
public: 
    MeshEncoder() = default;

    std::vector<uint8_t> EncodeRootLevelMesh(t8_forest_t root_level_mesh) override;    

#ifdef CMC_ENABLE_MPI
    std::vector<uint8_t> EncodeRootLevelMeshPar(t8_forest_t root_level, const MPI_Comm comm) override;
#endif
};


inline std::vector<uint8_t>
MeshEncoder::EncodeRootLevelMesh([[maybe_unused]] t8_forest_t root_level_mesh)
{
    size_t num_bytes_root_mesh_encoding_ = 1;

    std::vector<uint8_t> encoded_root_level;
    encoded_root_level.reserve(sizeof(size_t) + 1);

    /* Push the number of root level bytes to the stream */
    PushBackValueToByteStream(encoded_root_level, num_bytes_root_mesh_encoding_);
    encoded_root_level.push_back(0);

    return encoded_root_level;
}


#ifdef CMC_ENABLE_MPI
inline std::vector<uint8_t>
MeshEncoder::EncodeRootLevelMeshPar([[maybe_unused]] t8_forest_t root_level, const MPI_Comm comm)
{
    #ifdef CMC_ENABLE_DEBUG
    MPI_Comm mesh_comm = t8_forest_get_mpicomm(root_level);
    int are_mpi_comms_equal{0};
    const int ret_val_comm_compare = MPI_Comm_compare(comm, mesh_comm, &are_mpi_comms_equal);
    MPICheckError(ret_val_comm_compare);
    if (are_mpi_comms_equal != MPI_IDENT)
    {
        cmc_err_msg("The MPI-Communicator of the mesh and the variable differs. Therefore, no compression can be applied.");
    }
    #endif

    int rank{0};
    const int rv_rank = MPI_Comm_rank(comm, &rank);
    MPICheckError(rv_rank);

    /* Define a root rank */
    const int root_rank{0};

    /* Declare an encoded output vector (only to be filled by the root rank) */
    std::vector<uint8_t> encoded_root_level;

    if (rank == root_rank)
    {
        encoded_root_level.reserve(sizeof(uint64_t));

        /* Only the root rank stores the number of root elements */
        const uint64_t num_root_elements = static_cast<uint64_t>(t8_forest_get_global_num_leaf_elements(root_level));
        PushBackValueToByteStream(encoded_root_level, num_root_elements);
    }

    return encoded_root_level;
}
#endif

}


#endif /* !CMC_MESH_ENCODER_HXX */