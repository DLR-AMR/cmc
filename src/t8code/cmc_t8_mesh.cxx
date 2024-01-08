#include "t8code/cmc_t8_mesh.hxx"

namespace cmc {

static
std::array<int, 3>
GetElementAnchorOfElement(const t8_element_t* element, const t8_eclass_scheme_c* ts, const AmrMesh& mesh)
{
    //Depenedent on the initial mesh the Element anchor has to be adapted (embedded or congrunt mesh)
    // In case of an embedded mesh, we have only one tree and nothing else has to be done except the below
    std::array<int, 3> element_anchor;

    /* Get the element scheme */
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    /* Receive the integer anchor coordinates of the element */
    ts_c->t8_element_anchor (element, element_anchor.data());

    return elemet_anchor;
}

bool
IsMeshElementWithinGeoDomain(const AmrMesh& mesh, const t8_element_t* element, t8_eclass_scheme_c* ts, const GeoDomain& reference_domain)
{
    #ifdef CMC_WITH_T8CODE
    cmc_assert(reference_domain->GetDimensionality() == 2 || reference_domain->GetDimensionality() == 3);

    const int dimensionality = reference_domain->GetDimensionality();

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);

    /* Get the anchor coordinates of the element */
    std::array<int, 3> element_anchor = GetElementAnchorOfElement(element, element_anchor, mesh);

    const int initial_refinement_level = mesh.GetInitialRefinementLevel();

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
        switch (reference_domain.GetDataLayout())
        {
            case CMC_2D_LAT_LON:
                [[fallthrough]];
            case CMC_2D_LON_LAT:
                if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max &&
                    element_anchor[1] >= lat_min && element_anchor[1] < lat_max)
                {
                    /* The 2D element is inside the "lon x lat" mesh */
                    return true;
                }
            break;
            case CMC_2D_LAT_LEV:
                [[fallthrough]];
            case CMC_2D_LEV_LAT:
                if (element_anchor[0] >= lat_min && element_anchor[0] < lat_max &&
                    element_anchor[1] >= lev_min && element_anchor[1] < lev_max)
                {
                    /* The 2D element is inside the "lev x lat" mesh */
                    return true;
                }
            break;
            case CMC_2D_LON_LEV:
                [[fallthrough]];
            case CMC_2D_LEV_LON:
                if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max &&
                    element_anchor[1] >= lev_min && element_anchor[1] < lev_max)
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
        if (element_anchor[0] >= lon_min && element_anchor[0] < lon_max &&
            element_anchor[1] >= lat_min && element_anchor[1] < lat_max &&
            element_anchor[2] >= lev_min && element_anchor[2] < lev_max)
        {
            /* The 3D is inside the "lat x lon x lev" mesh */
            return true;
        }
    }
    /* The element is not inside the "lat x lon x lev" mesh */
    return false;

    #else
    return static_cast<bool>(CMC_ERR);
    #endif
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
