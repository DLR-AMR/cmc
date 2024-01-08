#ifndef CMC_T8_MESH_HXX
#define CMC_T8_MESH_HXX
/**
 * @file cmc_t8_mesh.hxx
 */


#ifdef CMC_WITH_T8CODE
#include "forest/t8_forest_general.h"
#endif

namespace cmc {

constexpr int kInitialRefinementLevelIsUnknown = INT_MIN;
constexpr int kMeshCorrespondsToNoneVariables = INT_MIN;
constexpr int kMeshCorrespondsToAllVariables = INT_MIN + 1;

class AmrMesh
{
public:
    AmrMesh(){};
    ~AmrMesh(){
        if (mesh_ != nullptr)
        {
            /* Deallocate the mesh (if there is one) */
            t8_forest_unref(&mesh_);
        }
    };

    AmrMesh(const AmrMesh& other)
    : mesh{other.mesh_},
      initial_refinement_level_{other.initial_refinement_level_},
      corresponding_variable_id{other.corresponding_variable_id_}
    {
        if (other.mesh_ != nullptr)
        {
            t8_forest_ref(other.mesh_);
        }
    }

    AmrMesh& operator=(const AmrMesh& other)
    {
        if (mesh_ != nullptr)
        {
            t8_forest_unref(mesh_);
        }
        return *this = AmrMesh(other);
    }

    AmrMesh(AmrMesh&& other) = default;
    AmrMesh& operator=(AmrMesh&& other) = default;

    bool IsValid() const;
    int GetInitialRefinementLevel() const;
    int SetInitialRefinementLevel(const int initial_refinement_level);
    t8_gloidx_t GetNumberGlobalElements() const;
    t8_locidx_t GetNumberLocalElements() const;

    t8_forest_t GetMesh() const;
    void SetMesh(t8_forest_t mesh);
    
private:
    t8_forest_t mesh_{nullptr};
    int initial_refinement_level_{kInitialRefinementLevelIsUnknown};
    //maybe an enum indicating whether it is a congruent or embedded mesh
};

struct CoarseningSample
{
    CoarseningSample(const int variable_id)
    : corresponding_variable_id{variable_id}{};
    const int corresponding_variable_id{kMeshCorrespondsToNoneVariables};
};

}

bool
IsMeshElementWithinGeoDomain(const GeoDomain& reference_domain, const AmrMesh& mesh, const t8_element_t* element, t8_eclass_scheme_c* ts);

bool
IsAnyElementWithinGeoDomain(const GeoDomain& reference_domain, const AmrMesh& mesh, const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c* ts);

#endif /* !CMC_T8_MESH_HXX */
