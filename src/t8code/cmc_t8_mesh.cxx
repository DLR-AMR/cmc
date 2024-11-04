#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "t8code/cmc_t8_morton.hxx"
#include "t8code/cmc_t8_adapt_callbacks.hxx"

#if CMC_WITH_T8CODE
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.hxx>
#include "t8_schemes/t8_default/t8_default_c_interface.h"
#include "t8_element_c_interface.h"
#include <p4est.h>
#include <p8est.h>
#endif

namespace cmc
{

t8_eclass_t
DimensionToElementClass(const int dimensionality)
{
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    return (dimensionality == 3 ? T8_ECLASS_HEX : T8_ECLASS_QUAD);
}

static
std::array<int, 3>
GetElementAnchorOfElement(const t8_element_t* element, t8_eclass_scheme_c* ts)
{
    std::array<int, 3> element_anchor;

    /* Get the element scheme */
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    /* Receive the integer anchor coordinates of the element */
    ts_c->t8_element_anchor (element, element_anchor.data());

    return element_anchor;
}
bool
IsMeshElementWithinGlobalDomain(const t8_element_t* element, t8_eclass_scheme_c* ts, const GeoDomain& global_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(global_domain.GetDimensionality() == 2 || global_domain.GetDimensionality() == 3);

    const int dimensionality = global_domain.GetDimensionality();

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    /* Get the anchor coordinates of the element */
    std::array<int, 3> element_anchor = GetElementAnchorOfElement(element, ts);

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
IsMeshElementWithinGeoDomain(const t8_element_t* element, t8_eclass_scheme_c* ts, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(reference_domain.GetDimensionality() == 2 || reference_domain.GetDimensionality() == 3);

    const int dimensionality = reference_domain.GetDimensionality();

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    /* Get the anchor coordinates of the element */
    std::array<int, 3> element_anchor = GetElementAnchorOfElement(element, ts);

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
                if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lon) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lon) &&
                    element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lat) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lat))
                {
                    /* The 2D element is inside the "lon x lat" mesh */
                    return true;
                }
            break;
            case DataLayout::Lat_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lat:
                if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lat) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lat) &&
                    element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lev) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lev))
                {
                    /* The 2D element is inside the "lev x lat" mesh */
                    return true;
                }
            break;
            case DataLayout::Lon_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lon:
                if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lon) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lon) &&
                    element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lev) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lev))
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
        if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lon) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lon) &&
            element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lat) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lat) &&
            element_anchor[2] >= reference_domain.GetDimensionStartIndex(Dimension::Lev) && element_anchor[2] < reference_domain.GetDimensionEndIndex(Dimension::Lev))
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
IsAnyElementWithinGlobalDomain(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c* ts, const GeoDomain& global_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    cmc_assert(global_domain.GetDimensionality() == 2 || global_domain.GetDimensionality() == 3);

    for (int elem_id = 0; elem_id < num_elements; ++elem_id)
    {
        if (IsMeshElementWithinGlobalDomain(elements[elem_id], ts, global_domain, initial_refinement_level, initial_layout))
        {
            return true;
        }
    }

    /* The element is not inside the ("lat x lon x lev") mesh */
    return false;
}


bool
IsAnyElementWithinGeoDomain(const int num_elements, const t8_element_t* elements[], t8_eclass_scheme_c* ts, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    cmc_assert(reference_domain.GetDimensionality() == 2 || reference_domain.GetDimensionality() == 3);

    //TODO:remove an define ar more abstract class of geodomains for the different element types
    //return true;//remove, just for tests now

    const int dimensionality = reference_domain.GetDimensionality();

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    for (int elem_id = 0; elem_id < num_elements; ++elem_id)
    {
        /* Get the anchor coordinates of the element */
        std::array<int, 3> element_anchor = GetElementAnchorOfElement(elements[elem_id], ts);

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
                    if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lon) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lon) &&
                        element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lat) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lat))
                    {
                        /* The 2D element is inside the "lon x lat" mesh */
                        return true;
                    }
                break;
                case DataLayout::Lat_Lev:
                    [[fallthrough]];
                case DataLayout::Lev_Lat:
                    if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lat) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lat) &&
                        element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lev) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lev))
                    {
                        /* The 2D element is inside the "lev x lat" mesh */
                        return true;
                    }
                break;
                case DataLayout::Lon_Lev:
                    [[fallthrough]];
                case DataLayout::Lev_Lon:
                    if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lon) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lon) &&
                        element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lev) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lev))
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
            if (element_anchor[0] >= reference_domain.GetDimensionStartIndex(Dimension::Lon) && element_anchor[0] < reference_domain.GetDimensionEndIndex(Dimension::Lon) &&
                element_anchor[1] >= reference_domain.GetDimensionStartIndex(Dimension::Lat) && element_anchor[1] < reference_domain.GetDimensionEndIndex(Dimension::Lat) &&
                element_anchor[2] >= reference_domain.GetDimensionStartIndex(Dimension::Lev) && element_anchor[2] < reference_domain.GetDimensionEndIndex(Dimension::Lev))
            {
                /* The 3D is inside the "lat x lon x lev" mesh */
                return true;
            }
        }
    }

    /* The element is not inside the ("lat x lon x lev") mesh */
    return false;
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
AmrMesh::GetNumberGlobalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_global_num_elements(mesh_);
}

t8_locidx_t
AmrMesh::GetNumberLocalElements() const
{
    cmc_assert(mesh_ != nullptr);
    return t8_forest_get_local_num_elements(mesh_);
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

void
AmrMesh::IndicateWhetherDummyElementsArePresent(const bool are_dummy_elements_present)
{
    are_dummy_elements_present_ = are_dummy_elements_present;
}

bool
AmrMesh::AreDummyElementsPresent() const
{
    return are_dummy_elements_present_;
}

MortonIndex
GetMortonIndexOnLevel(const t8_element_t* elem, t8_eclass_scheme_c* ts, const int dimensionality, const int level)
{
    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    int element_anchor[3];

    /* Get the element scheme */
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    /* Receive the integer anchor coordinates of the element */
    ts_c->t8_element_anchor (elem, element_anchor);

    std::vector<DomainIndex> coords;

    /* Transform coordinates into the range of the initial-refinement-level coordinate values */
    for (int index{0}; index < dimensionality; ++index)
    {
        coords.emplace_back(element_anchor[index] >> (element_anchor_max_lvl - level));
    }

    return GetMortonIndex(coords, dimensionality);
}

int
GetFirstElementIDOnReferenceLevel(const t8_element_t* element, t8_eclass_scheme_c* ts, const int reference_level)
{
    cmc_assert(t8_element_level(ts, element) <= reference_level);
    //const int element_level = t8_element_level(ts, element);
    //Get Linear ID oder get anchor coords -> shift in range -> GetMorton
    cmc_err_msg("Not implmeneted yet");
    return -1;
}

int
DetermineNumberOfElementsOnReferenceLevel(const t8_element_t* element, t8_eclass_scheme_c* ts, const int num_children, const int reference_level)
{
    cmc_assert(t8_element_level(ts, element) <= reference_level);
    return std::pow(num_children, reference_level - t8_element_level(ts, element));
}

int 
DetermineNumberOfDecompressedElements(const bool restrict_to_domain, const t8_element_t* element, t8_eclass_scheme_c* ts, const int num_children, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    if (!restrict_to_domain)
    {
        //return std::pow(num_children, initial_refinement_level - t8_element_level(ts, element));
        return DetermineNumberOfElementsOnReferenceLevel(element, ts, num_children, initial_refinement_level);
    }

    const int dimensionality = reference_domain.GetDimensionality();

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    const int element_level = t8_element_level(ts, element);

    const int max_decompressed_elems_per_dimension = std::pow(2, initial_refinement_level - element_level);

    /* Get the anchor coordinates of the element */
    std::array<int, 3> element_anchor = GetElementAnchorOfElement(element, ts);

    /* Transform coordinates into the range of the initial-refinement-level coordinate values */
    for (int index{0}; index < dimensionality; ++index)
    {
        element_anchor[index] >>= (element_anchor_max_lvl - initial_refinement_level);
    }
        
    if (dimensionality == 2)
    {
        /* 2D case */
        switch (initial_layout)
        {
            case DataLayout::Lat_Lon:
                [[fallthrough]];
            case DataLayout::Lon_Lat:
            {
                const int lon_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lon) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lon) - element_anchor[0]);
                const int lat_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lat) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lat) - element_anchor[1]);
                const int lon_max = std::min(lon_spread, max_decompressed_elems_per_dimension);
                const int lat_max = std::min(lat_spread, max_decompressed_elems_per_dimension);
                return lon_max * lat_max;
            }
            break;
            case DataLayout::Lat_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lat:
            {
                const int lat_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lat) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lat) - element_anchor[0]);
                const int lev_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lev) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lev) - element_anchor[1]);
                const int lat_max = std::min(lat_spread, max_decompressed_elems_per_dimension);
                const int lev_max = std::min(lev_spread, max_decompressed_elems_per_dimension);
                return lat_max * lev_max;
            }
            break;
            case DataLayout::Lon_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lon:
            {
                const int lon_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lon) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lon) - element_anchor[0]);
                const int lev_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lev) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lev) - element_anchor[1]);
                const int lon_max = std::min(lon_spread, max_decompressed_elems_per_dimension);
                const int lev_max = std::min(lev_spread, max_decompressed_elems_per_dimension);
                return lon_max * lev_max;
            }
            break;
            default:
                cmc_err_msg("There was no valid 2D data layout supplied to determine whether or not the element is within the domain.");
        }

    } else {
        /* 3D case */
        const int lon_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lon) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lon) - element_anchor[0]);
        const int lat_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lat) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lat) - element_anchor[1]);
        const int lev_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lev) <= element_anchor[2] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lev) - element_anchor[2]);
        const int lon_max = std::min(lon_spread, max_decompressed_elems_per_dimension);
        const int lat_max = std::min(lat_spread, max_decompressed_elems_per_dimension);
        const int lev_max = std::min(lev_spread, max_decompressed_elems_per_dimension);
        return lon_max * lat_max * lev_max;
    }

    return 0;
}




Hyperslab
DetermineHyperslabOfDecompressedElements(const t8_element_t* element, t8_eclass_scheme_c* ts, const GeoDomain& reference_domain, const int initial_refinement_level, const DataLayout initial_layout)
{
    const int dimensionality = reference_domain.GetDimensionality();

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    const int element_level = t8_element_level(ts, element);

    const int max_decompressed_elems_per_dimension = std::pow(2, initial_refinement_level - element_level);

    /* Get the anchor coordinates of the element */
    std::array<int, 3> element_anchor = GetElementAnchorOfElement(element, ts);

    /* Transform coordinates into the range of the initial-refinement-level coordinate values */
    for (int index{0}; index < dimensionality; ++index)
    {
        element_anchor[index] >>= (element_anchor_max_lvl - initial_refinement_level);
    }
    
    if (dimensionality == 2)
    {
        /* 2D case */
        switch (initial_layout)
        {
            case DataLayout::Lat_Lon:
                [[fallthrough]];
            case DataLayout::Lon_Lat:
            {    
                const int lon_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lon) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lon) - element_anchor[0]);
                const int lat_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lat) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lat) - element_anchor[1]);
                const int lon_max = std::min(lon_spread, max_decompressed_elems_per_dimension);
                const int lat_max = std::min(lat_spread, max_decompressed_elems_per_dimension);
                return Hyperslab(DimensionInterval(Dimension::Lon, element_anchor[0], element_anchor[0] + lon_max),
                                 DimensionInterval(Dimension::Lat, element_anchor[1], element_anchor[1] + lat_max));
            }
            break;
            case DataLayout::Lat_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lat:
            {
                const int lat_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lat) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lat) - element_anchor[0]);
                const int lev_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lev) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lev) - element_anchor[1]);
                const int lat_max = std::min(lat_spread, max_decompressed_elems_per_dimension);
                const int lev_max = std::min(lev_spread, max_decompressed_elems_per_dimension);
                return Hyperslab(DimensionInterval(Dimension::Lat, element_anchor[0], element_anchor[0] + lat_max),
                                 DimensionInterval(Dimension::Lev, element_anchor[1], element_anchor[1] + lev_max));
            }
            break;
            case DataLayout::Lon_Lev:
                [[fallthrough]];
            case DataLayout::Lev_Lon:
            {
                const int lon_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lon) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lon) - element_anchor[0]);
                const int lev_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lev) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lev) - element_anchor[1]);
                const int lon_max = std::min(lon_spread, max_decompressed_elems_per_dimension);
                const int lev_max = std::min(lev_spread, max_decompressed_elems_per_dimension);
                return Hyperslab(DimensionInterval(Dimension::Lon, element_anchor[0], element_anchor[0] + lon_max),
                                 DimensionInterval(Dimension::Lev, element_anchor[1], element_anchor[1] + lev_max));
            }
            break;
            default:
                cmc_err_msg("There was no valid 2D data layout supplied to determine whether or not the element is within the domain.");
        }

    } else {
        /* 3D case */
        const int lon_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lon) <= element_anchor[0] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lon) - element_anchor[0]);
        const int lat_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lat) <= element_anchor[1] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lat) - element_anchor[1]);
        const int lev_spread = (reference_domain.GetDimensionEndIndex(Dimension::Lev) <= element_anchor[2] ? 0 : reference_domain.GetDimensionEndIndex(Dimension::Lev) - element_anchor[2]);
        const int lon_max = std::min(lon_spread, max_decompressed_elems_per_dimension);
        const int lat_max = std::min(lat_spread, max_decompressed_elems_per_dimension);
        const int lev_max = std::min(lev_spread, max_decompressed_elems_per_dimension);
        return Hyperslab(DimensionInterval(Dimension::Lon, element_anchor[0], element_anchor[0] + lon_max),
                         DimensionInterval(Dimension::Lat, element_anchor[1], element_anchor[1] + lat_max),
                         DimensionInterval(Dimension::Lev, element_anchor[2], element_anchor[2] + lev_max));
    }

    return Hyperslab();
}

t8_forest_t
ReconstructBaseMesh(const int dimensionality, const MPI_Comm comm)
{
    /* Get a base cmesh corresponding to the dimensionality */
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(DimensionToElementClass(dimensionality), comm, 0, 0, 1);

    /* Construct a forest from the cmesh */
    t8_forest_t base_mesh;
    t8_forest_init(&base_mesh);
    t8_forest_set_cmesh(base_mesh, cmesh, comm);
    t8_forest_set_scheme(base_mesh, t8_scheme_new_default_cxx());
    t8_forest_set_level(base_mesh, 0);
    t8_forest_commit(base_mesh);

    return base_mesh;
}

t8_forest_t
ReconstructMeshFromRefinementBits(const VectorView<uint8_t>& refinement_bits, const int dimensionality, const MPI_Comm comm)
{
    /* Get a base cmesh corresponding to the dimensionality */
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(DimensionToElementClass(dimensionality), comm, 0, 0, 1);

    /* Construct a forest from the cmesh */
    t8_forest_t decompressed_forest;
    t8_forest_init(&decompressed_forest);
    t8_forest_set_cmesh(decompressed_forest, cmesh, comm);
    t8_forest_set_scheme(decompressed_forest, t8_scheme_new_default_cxx());
    t8_forest_set_level(decompressed_forest, 0);
    t8_forest_commit(decompressed_forest);

    const size_t num_global_bytes = refinement_bits.size();
    size_t used_bytes = 0;

    /* Get a pointer to the start of the underlying array */
    const uint8_t* initial_data_ptr = refinement_bits.data();

    while (used_bytes < num_global_bytes)
    {
        const t8_gloidx_t num_global_elements = t8_forest_get_global_num_elements(decompressed_forest);

        /* Get the number of bytes needed for this decompression/refinement step */
        const t8_gloidx_t num_encoded_bytes = num_global_elements / 8 + (num_global_elements % 8 != 0);

        DecompressionRefinementBits adapt_data(VectorView<uint8_t>(initial_data_ptr + used_bytes, num_encoded_bytes));

        decompressed_forest = t8_forest_new_adapt(decompressed_forest, ApplyRefinementBits, 0, 0, static_cast<void*>(&adapt_data)); 

        used_bytes += num_encoded_bytes;
    }

    return decompressed_forest;
}

static
int
DetermineDimensionalityOfTheData(const GeoDomain& global_domain)
{
    return global_domain.GetDimensionality();
}

static
int
CalculateInitialRefinementLevel(const GeoDomain& global_domain)
{
    const size_t max_elem_per_direction = global_domain.GetLargestDimensionLength();
    /* Calculate the induced initial refinement level needed in order to build an enclosing mesh */
    return static_cast<int>(std::ceil(std::log2(max_elem_per_direction) + std::numeric_limits<double>::epsilon()));
}

static
void ValidateInitialRefinementLevelForDimensionality(const int initial_refinement_level, const int dimensionality)
{
    #ifdef CMC_ENABLE_DEBUG
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    
    const int maximum_possible_refinement_level = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);
    
    if (initial_refinement_level < 0 ||
        initial_refinement_level > maximum_possible_refinement_level)
    {
        cmc_err_msg("The corresponding refinement level is not within the range of an computationally posiible refinement level.");
    }

    #endif
}

std::tuple<t8_forest_t, int, int>
BuildInitialMesh(const GeoDomain& domain, const DataLayout initial_layout, MPI_Comm comm)
{
    const int dimensionality = DetermineDimensionalityOfTheData(domain);

    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(DimensionToElementClass(dimensionality), comm, 0, 0, 1);

    const int initial_refinement_level = CalculateInitialRefinementLevel(domain);
    cmc_debug_msg("initial ref level is: ", initial_refinement_level);
    ValidateInitialRefinementLevelForDimensionality(initial_refinement_level, dimensionality);

    /* Construct a forest from the cmesh */
    t8_forest_t initial_forest;
    t8_forest_init(&initial_forest);
    t8_forest_set_cmesh(initial_forest, cmesh, comm);
    t8_forest_set_scheme(initial_forest, t8_scheme_new_default_cxx());
    t8_forest_set_level(initial_forest, 0);
    t8_forest_commit(initial_forest);
    cmc_debug_msg("Constructed initial forest with: ", t8_forest_get_local_num_elements(initial_forest));
    AdaptDataInitialMesh adapt_data(domain, initial_refinement_level, initial_layout);

    for (int i = 0; i < 4; ++i)
    {
        cmc_debug_msg("Domain length at ", i, " is ", domain.GetDimensionLength(static_cast<Dimension>(i)));
        cmc_debug_msg("Start index: ", domain.GetDimensionStartIndex(static_cast<Dimension>(i)), " and End index: ", domain.GetDimensionEndIndex(static_cast<Dimension>(i)));
    }
    //TODO: Recursive construction without repartitioning may not be optimal in parallel
    /* Build the initial mesh via recursive refinement */
    //initial_forest = t8_forest_new_adapt(initial_forest, RefineToInitialMesh, 1, 0, static_cast<void*>(&adapt_data));

    t8_forest_t adapted_forest;
    for (int adaptation_steps = 0; adaptation_steps <= initial_refinement_level; ++adaptation_steps)
    {
        cmc_debug_msg("Iteration: ", adaptation_steps, "\n\n\n\n");
        t8_forest_init(&adapted_forest);
        t8_forest_set_adapt(adapted_forest, initial_forest, RefineToInitialMesh, 0);
        const int set_partition_for_coarsening = 0; //TODO change to one later
        t8_forest_set_partition(adapted_forest, NULL, set_partition_for_coarsening);
        t8_forest_set_user_data(adapted_forest, static_cast<void*>(&adapt_data));
        t8_forest_commit(adapted_forest);

        initial_forest = adapted_forest;
    }

    ///* Repartition the forest */
    //t8_forest_t partitioned_forest;
    //t8_forest_init(&partitioned_forest);
    //const int partition_for_coarsening = 0; //TODO: change to one in the future
    //t8_forest_set_partition(partitioned_forest, initial_forest, partition_for_coarsening);
    //t8_forest_commit(partitioned_forest);

    cmc_debug_msg("Num local elements: ", t8_forest_get_local_num_elements(initial_forest));
    return std::make_tuple(initial_forest, initial_refinement_level, dimensionality); 
}

}
