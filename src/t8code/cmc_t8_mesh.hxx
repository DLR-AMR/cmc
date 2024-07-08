#ifndef CMC_T8_MESH_HXX
#define CMC_T8_MESH_HXX
/**
 * @file cmc_t8_mesh.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_span.hxx"

#include <climits>

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#endif

namespace cmc {

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
    ~AmrMesh(){
        if (mesh_ != nullptr)
        {
            /* Deallocate the mesh (if there is one) */
            t8_forest_unref(&mesh_);
        }
    };

    AmrMesh(const AmrMesh& other)
    : mesh_{other.mesh_},
      initial_refinement_level_{other.initial_refinement_level_},
      dimensionality_{other.dimensionality_},
      are_dummy_elements_present_{other.are_dummy_elements_present_}
    {
        if (other.mesh_ != nullptr)
        {
            t8_forest_ref(other.mesh_);
        }
    };

    AmrMesh& operator=(const AmrMesh& other)
    {
        if (mesh_ != nullptr)
        {
            t8_forest_unref(&mesh_);
        }
        std::cout << std::endl;
        return *this = AmrMesh(other);
    };

    AmrMesh(AmrMesh&& other)
    : mesh_{std::move(other.mesh_)}, initial_refinement_level_{other.initial_refinement_level_},
      dimensionality_{other.dimensionality_}, are_dummy_elements_present_{other.are_dummy_elements_present_}
    {
        other.mesh_ = nullptr;
    };
    AmrMesh& operator=(AmrMesh&& other)
    {
        this->mesh_ = std::move(other.mesh_);
        other.mesh_ = nullptr;
        this->initial_refinement_level_ = other.initial_refinement_level_;
        this->dimensionality_ = other.dimensionality_;
        this->are_dummy_elements_present_ = other.are_dummy_elements_present_;
        return *this;
    };

    bool IsValid() const;
    int GetInitialRefinementLevel() const;
    void SetInitialRefinementLevel(const int initial_refinement_level);
    int GetDimensionality() const;
    void SetDimensionality(const int dimensionality);
    t8_gloidx_t GetNumberGlobalElements() const;
    t8_locidx_t GetNumberLocalElements() const;
    void IndicateWhetherDummyElementsArePresent(const bool are_dummy_elements_present);
    bool AreDummyElementsPresent() const;
    t8_forest_t GetMesh() const;
    void SetMesh(t8_forest_t mesh);
    
private:
    t8_forest_t mesh_{nullptr};
    int initial_refinement_level_{kInitialRefinementLevelIsUnknown};
    int dimensionality_{kDimensionalityIsUnknown};
    bool are_dummy_elements_present_{true};
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
IsMeshElementWithinGeoDomain(const t8_element_t* element, t8_eclass_scheme_c* ts, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout);

bool
IsAnyElementWithinGeoDomain(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c* ts, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout);

int
DetermineNumberOfElementsOnReferenceLevel(const t8_element_t* element, t8_eclass_scheme_c* ts, const int num_children, const int reference_level);

int
DetermineNumberOfDecompressedElements(const bool restrict_to_domain, const t8_element_t* element, t8_eclass_scheme_c* ts, const int num_children, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout);

Hyperslab
DetermineHyperslabOfDecompressedElements(const t8_element_t* element, t8_eclass_scheme_c* ts, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout);

MortonIndex GetMortonIndexOnLevel(const t8_element_t* elem, t8_eclass_scheme_c* ts, const int dimensioanlity, const int level);

t8_forest_t
ReconstructMeshFromRefinementBits(const VectorView<uint8_t>& refinement_bits, const int dimensionality, const MPI_Comm comm);

t8_forest_t
ReconstructBaseMesh(const int dimensionality, const MPI_Comm comm);

std::tuple<t8_forest_t, int, int>
BuildInitialMesh(const GeoDomain& domain, const DataLayout initial_layout, MPI_Comm comm);

}

#endif /* !CMC_T8_MESH_HXX */
