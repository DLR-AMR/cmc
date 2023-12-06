#include "t8code/cmc_t8_mesh.hxx"

namespace cmc {

inline t8_forest_t
AmrMesh::GetMesh() const
{
    cmc_assert(mesh_ != nullptr);
    return mesh_;    
}

inline void
AmrMesh::SetMesh(t8_forest_t mesh)
{
    cmc_assert(mesh != nullptr);
    mesh_ = mesh;    
}

bool
AmrMesh::IsValid()
{
    return (mesh_ != nullptr ? true : false);
}

inline t8_gloidx_t
AmrMesh::GetNumberGlobalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_global_num_elements(mesh_);
}

inline t8_locidx_t
AmrMesh::GetNumberLocalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_local_num_elements(mesh_);
}

inline int
AmrMesh::GetInitialRefinementLevel() const
{
    return initial_refinement_level_;
}

inline int
AmrMesh::SetInitialRefinementLevel(const int initial_refinement_level)
{
    initial_refinement_level_ = initial_refinement_level;
}

}
