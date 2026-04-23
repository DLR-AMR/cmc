#ifndef CMC_T8_MESH_HXX
#define CMC_T8_MESH_HXX
/**
 * @file cmc_t8_mesh.hxx
 */

#include "utilities/cmc_utilities.hxx"

#include <climits>

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx> 
#include <t8_forest/t8_forest_iterate.h>
#include <t8_eclass/t8_eclass.h>
#include <t8_element/t8_element.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#endif

namespace cmc
{

constexpr int kInitialRefinementLevelIsUnknown = INT_MIN;
constexpr int kMeshCorrespondsToNoneVariables = INT_MIN;
constexpr int kMeshCorrespondsToAllVariables = INT_MIN + 1;
constexpr int kDimensionalityIsUnknown = INT_MIN;

class AmrMesh
{
public:
    AmrMesh() = default;
    AmrMesh(t8_forest_t mesh)
    : mesh_{mesh} {};
    AmrMesh(t8_forest_t mesh, const int initial_refinement_level)
    : mesh_{mesh}, initial_refinement_level_{initial_refinement_level} {};
    AmrMesh(t8_forest_t mesh, const int initial_refinement_level, const int dimensionality)
    : mesh_{mesh}, initial_refinement_level_{initial_refinement_level}, dimensionality_{dimensionality} {};
    ~AmrMesh();

    AmrMesh(const AmrMesh& other);
    AmrMesh& operator=(const AmrMesh& other);
    AmrMesh(AmrMesh&& other);
    AmrMesh& operator=(AmrMesh&& other);

    bool IsValid() const;
    int GetInitialRefinementLevel() const;
    void SetInitialRefinementLevel(const int initial_refinement_level);
    int GetDimensionality() const;
    void SetDimensionality(const int dimensionality);
    t8_gloidx_t GetNumberGlobalElements() const;
    t8_locidx_t GetNumberLocalElements() const;
    t8_gloidx_t GetNumberGlobalTrees() const; 
    t8_forest_t GetMesh() const;
    void SetMesh(t8_forest_t mesh);
    void SetNullMesh();
private:
    t8_forest_t mesh_{nullptr};
    int initial_refinement_level_{kInitialRefinementLevelIsUnknown};
    int dimensionality_{kDimensionalityIsUnknown};
};

inline
AmrMesh::~AmrMesh()
{
    if (mesh_ != nullptr)
    {
        /* Deallocate the mesh (if there is one) */
        t8_forest_unref(&mesh_);
    }
}

inline
AmrMesh::AmrMesh(const AmrMesh& other)
: mesh_{other.mesh_},
  initial_refinement_level_{other.initial_refinement_level_},
  dimensionality_{other.dimensionality_}
{
    if (other.mesh_ != nullptr)
    {
        t8_forest_ref(other.mesh_);
    }
}

inline AmrMesh&
AmrMesh::operator=(const AmrMesh& other)
{
    if (mesh_ != nullptr)
    {
        t8_forest_unref(&mesh_);
    }
    std::cout << std::endl;
    return *this = AmrMesh(other);
}

inline
AmrMesh::AmrMesh(AmrMesh&& other)
: mesh_{std::move(other.mesh_)}, initial_refinement_level_{other.initial_refinement_level_},
  dimensionality_{other.dimensionality_}
{
    other.mesh_ = nullptr;
}

inline AmrMesh&
AmrMesh::operator=(AmrMesh&& other)
{
    this->mesh_ = std::move(other.mesh_);
    other.mesh_ = nullptr;
    this->initial_refinement_level_ = other.initial_refinement_level_;
    this->dimensionality_ = other.dimensionality_;
    return *this;
}
    
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

inline void
AmrMesh::SetNullMesh()
{
    mesh_ = nullptr; 
}

inline int
AmrMesh::GetDimensionality() const
{
    return dimensionality_;
}

inline void
AmrMesh::SetDimensionality(const int dimensionality)
{
    cmc_assert(dimensionality >= 2 && dimensionality <= 3);
    dimensionality_ = dimensionality;
}

inline bool
AmrMesh::IsValid() const
{
    return (mesh_ != nullptr ? true : false);
}

inline t8_gloidx_t
AmrMesh::GetNumberGlobalTrees() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_num_global_trees(mesh_);
} 

inline t8_gloidx_t
AmrMesh::GetNumberGlobalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_global_num_leaf_elements(mesh_);
}

inline t8_locidx_t
AmrMesh::GetNumberLocalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_local_num_leaf_elements(mesh_);
}

inline int
AmrMesh::GetInitialRefinementLevel() const
{
    return initial_refinement_level_;
}

inline void
AmrMesh::SetInitialRefinementLevel(const int initial_refinement_level)
{
    initial_refinement_level_ = initial_refinement_level;
}


}

#endif /* !CMC_T8_MESH_HXX */
