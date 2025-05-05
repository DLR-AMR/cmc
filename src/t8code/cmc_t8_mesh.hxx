#ifndef CMC_T8_MESH_HXX
#define CMC_T8_MESH_HXX
/**
 * @file cmc_t8_mesh.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_hyperslab.hxx"

#include <climits>

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx> 
#include <t8_forest/t8_forest_iterate.h>
#include <t8_eclass.h>
#include <t8_element_cxx.hxx>
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

private:
    t8_forest_t mesh_{nullptr};
    int initial_refinement_level_{kInitialRefinementLevelIsUnknown};
    int dimensionality_{kDimensionalityIsUnknown};
};

struct CoarseningSample
{
    CoarseningSample(const int variable_id)
    : corresponding_variable_id{variable_id}{};
    const int corresponding_variable_id{kMeshCorrespondsToNoneVariables};
};

t8_eclass_t
DimensionToElementClass(const int dimensionality);

bool
IsMeshElementWithinGlobalDomain(const t8_element_t* element, t8_eclass_scheme_c* ts, const GeoDomain& global_domain, const int initial_refinement_level, const DataLayout initial_layout);

bool
IsAnyElementWithinGlobalDomain(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c* ts, const GeoDomain& global_domain, const int initial_refinement_level, const DataLayout initial_layout);

MortonIndex GetMortonIndexOnLevel(const t8_element_t* elem, t8_eclass_scheme_c* ts, const int dimensioanlity, const int level);

}

#endif /* !CMC_T8_MESH_HXX */
