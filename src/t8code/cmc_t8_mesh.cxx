#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_morton.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8_schemes/t8_default/t8_default.hxx>
#include <p4est.h>
#include <p8est.h>
#endif

#if 0
#ifdef CMC_WITH_T8CODE
#include <t8_eclass.h>
#include <t8_element.h>
#include <t8_schemes/t8_scheme.hxx>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.hxx>
#include "t8_schemes/t8_default/t8_default_c_interface.h"
#include "t8_element_c_interface.h"
#include <p4est.h>
#include <p8est.h>
#endif
#endif

namespace cmc
{

AmrMesh::~AmrMesh()
{
    if (mesh_ != nullptr)
    {
        /* Deallocate the mesh (if there is one) */
        t8_forest_unref(&mesh_);
    }
}

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

AmrMesh& AmrMesh::operator=(const AmrMesh& other)
{
    if (mesh_ != nullptr)
    {
        t8_forest_unref(&mesh_);
    }
    std::cout << std::endl;
    return *this = AmrMesh(other);
}

AmrMesh::AmrMesh(AmrMesh&& other)
: mesh_{std::move(other.mesh_)}, initial_refinement_level_{other.initial_refinement_level_},
  dimensionality_{other.dimensionality_}
{
    other.mesh_ = nullptr;
}

AmrMesh& AmrMesh::operator=(AmrMesh&& other)
{
    this->mesh_ = std::move(other.mesh_);
    other.mesh_ = nullptr;
    this->initial_refinement_level_ = other.initial_refinement_level_;
    this->dimensionality_ = other.dimensionality_;
    return *this;
}
    
t8_forest_t
AmrMesh::GetMesh() const
{
    cmc_assert(mesh_ != nullptr);
    return mesh_;    
}

void
AmrMesh::SetMesh(t8_forest_t mesh)
{
    cmc_assert(mesh != nullptr);
    mesh_ = mesh;    
}

int
AmrMesh::GetDimensionality() const
{
    return dimensionality_;
}

void
AmrMesh::SetDimensionality(const int dimensionality)
{
    cmc_assert(dimensionality >= 2 && dimensionality <= 3);
    dimensionality_ = dimensionality;
}

bool
AmrMesh::IsValid() const
{
    return (mesh_ != nullptr ? true : false);
}

t8_gloidx_t
AmrMesh::GetNumberGlobalTrees() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_num_global_trees(mesh_);
} 

t8_gloidx_t
AmrMesh::GetNumberGlobalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_global_num_leaf_elements(mesh_);
}

t8_locidx_t
AmrMesh::GetNumberLocalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_local_num_leaf_elements(mesh_);
}

int
AmrMesh::GetInitialRefinementLevel() const
{
    return initial_refinement_level_;
}

void
AmrMesh::SetInitialRefinementLevel(const int initial_refinement_level)
{
    initial_refinement_level_ = initial_refinement_level;
}


t8_eclass_t
DimensionToElementClass(const int dimensionality)
{
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    return (dimensionality == 3 ? T8_ECLASS_HEX : T8_ECLASS_QUAD);
}

static
std::array<int, 3>
GetElementAnchorOfElement(const t8_eclass_t tree_class, const t8_element_t* element, const t8_scheme_c* scheme)
{
    cmc_assert(t8_eclass_scheme_is_default(scheme, tree_class) != 0);

    std::array<int, 3> element_anchor;
    element_anchor.fill(0);
    //int* array_ptr = element_anchor.data();

    /* Receive the integer anchor coordinates of the element */
    scheme->element_get_anchor (tree_class, element, element_anchor.data());

    return element_anchor;
}

bool
IsMeshElementWithinGlobalDomain(const t8_eclass_t tree_class, const t8_element_t* element, const t8_scheme_c* ts, const GeoDomain& global_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(global_domain.GetDimensionality() == 2 || global_domain.GetDimensionality() == 3);

    const int dimensionality = global_domain.GetDimensionality();

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    /* Get the anchor coordinates of the element */
    std::array<int, 3> element_anchor = GetElementAnchorOfElement(tree_class, element, ts);

    /* Transform coordinates into the range of the initial-refinement-level coordinate values */
    for (int index{0}; index < dimensionality; ++index)
    {
        element_anchor[index] >>= (element_anchor_max_lvl - initial_refinement_level);
    }

    /* Check if the anchor coordinates of the element lie within the "lat x lon x lev" mesh */
    /** \note: Different geo-spatial coordinate combinations are corresponding to certain x- and y- coordinates. This labeling is consistent. 
     *         Data only concerning latitude and longitude will be ordered, s.t. longitude equals x, latitude equals y
     *         Data only concerning latitude and elevation will be ordered, s.t. latitude equals x, elevation equals y
     *         Data only concerning longitude and elevation will be ordered, s.t. longitude equals x, elevation equals y
    **/
    if (dimensionality == 2)
    {
        /* 2D case */
        switch (initial_layout)
        {
            case DataLayout::Lat_Lon:
                [[fallthrough]];
            case DataLayout::Lon_Lat:
                if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lon) &&
                    element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lat))
                {
                    /* The 2D element is inside the "lon x lat" mesh */
                    return true;
                }
            break;
            case DataLayout::Lat_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lat:
                if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lat) &&
                    element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lev))
                {
                    /* The 2D element is inside the "lev x lat" mesh */
                    return true;
                }
            break;
            case DataLayout::Lon_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lon:
                if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lon) &&
                    element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lev))
                {
                    /* The 2D element is inside the "lev x lon" mesh */
                    return true;
                }
            break;
            default:
                cmc_err_msg("There was no valid 2D data layout supplied to determine whether or not the element is within the domain.");
        }

    } else {
        /* 3D case */
        if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lon) &&
            element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lat) &&
            element_anchor[2] >= 0 && element_anchor[2] < global_domain.GetDimensionLength(Dimension::Lev))
        {
            /* The 3D is inside the "lat x lon x lev" mesh */
            return true;
        }
    }

    /* The element is not inside the ("lat x lon x lev") mesh */
    return false;

    #else
    return static_cast<bool>(CMC_ERR);
    #endif
}

bool
IsMeshElementWithinGlobalDomain(const std::vector<DomainIndex>& element_anchor, const GeoDomain& global_domain, const DataLayout initial_layout)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(global_domain.GetDimensionality() == 2 || global_domain.GetDimensionality() == 3);

    const int dimensionality = global_domain.GetDimensionality();

    /* Check if the anchor coordinates of the element lie within the "lat x lon x lev" mesh */
    /** \note: Different geo-spatial coordinate combinations are corresponding to certain x- and y- coordinates. This labeling is consistent. 
     *         Data only concerning latitude and longitude will be ordered, s.t. longitude equals x, latitude equals y
     *         Data only concerning latitude and elevation will be ordered, s.t. latitude equals x, elevation equals y
     *         Data only concerning longitude and elevation will be ordered, s.t. longitude equals x, elevation equals y
    **/
    if (dimensionality == 2)
    {
        /* 2D case */
        switch (initial_layout)
        {
            case DataLayout::Lat_Lon:
                [[fallthrough]];
            case DataLayout::Lon_Lat:
                if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lon) &&
                    element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lat))
                {
                    /* The 2D element is inside the "lon x lat" mesh */
                    return true;
                }
            break;
            case DataLayout::Lat_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lat:
                if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lat) &&
                    element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lev))
                {
                    /* The 2D element is inside the "lev x lat" mesh */
                    return true;
                }
            break;
            case DataLayout::Lon_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lon:
                if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lon) &&
                    element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lev))
                {
                    /* The 2D element is inside the "lev x lon" mesh */
                    return true;
                }
            break;
            default:
                cmc_err_msg("There was no valid 2D data layout supplied to determine whether or not the element is within the domain.");
        }

    } else {
        /* 3D case */
        if (element_anchor[0] >= 0 && element_anchor[0] < global_domain.GetDimensionLength(Dimension::Lon) &&
            element_anchor[1] >= 0 && element_anchor[1] < global_domain.GetDimensionLength(Dimension::Lat) &&
            element_anchor[2] >= 0 && element_anchor[2] < global_domain.GetDimensionLength(Dimension::Lev))
        {
            /* The 3D is inside the "lat x lon x lev" mesh */
            return true;
        }
    }

    /* The element is not inside the ("lat x lon x lev") mesh */
    return false;

    #else
    return static_cast<bool>(CMC_ERR);
    #endif
}

bool
IsAnyElementWithinGlobalDomain(const t8_eclass_t tree_class, const int num_elements, const t8_element_t* elements[], const t8_scheme_c* ts, const GeoDomain& global_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    cmc_assert(global_domain.GetDimensionality() == 2 || global_domain.GetDimensionality() == 3);

    for (int elem_id = 0; elem_id < num_elements; ++elem_id)
    {
        if (IsMeshElementWithinGlobalDomain(tree_class, elements[elem_id], ts, global_domain, initial_refinement_level, initial_layout))
        {
            return true;
        }
    }

    /* The elements are not inside the global domain */
    return false;
}

MortonIndex
GetMortonIndexOnLevel(const t8_eclass_t tree_class, const t8_element_t* elem, const t8_scheme_c* ts, const int dimensionality, const int level)
{
    cmc_assert(t8_eclass_scheme_is_default(ts, tree_class) != 0);

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    int element_anchor[3];

    /* Receive the integer anchor coordinates of the element */
    ts->element_get_anchor (tree_class, elem, element_anchor);

    std::vector<DomainIndex> coords;

    /* Transform coordinates into the range of the initial-refinement-level coordinate values */
    for (int index{0}; index < dimensionality; ++index)
    {
        coords.emplace_back(element_anchor[index] >> (element_anchor_max_lvl - level));
    }

    return GetMortonIndex(coords, dimensionality);
}

}
