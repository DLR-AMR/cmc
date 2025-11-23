#include "utilities/cmc_morton.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_geo_utilities.hxx"

namespace cmc
{

static inline
MortonIndex
PrepareValueFor2DInterleaving(MortonIndex value)
{
    /** Prepare the value for interleaving by fragmenting the bits, like:
     * (...b4 b3 b2 b1 b0) -> (...0 b4 0 b3 0 b2 0 b1 0 b0) */
    value = (value | (value << 16)) & 0x0000FFFF0000FFFF;
    value = (value | (value << 8))  & 0x00FF00FF00FF00FF;
    value = (value | (value << 4))  & 0x0F0F0F0F0F0F0F0F;
    value = (value | (value << 2))  & 0x3333333333333333;
    value = (value | (value << 1))  & 0x5555555555555555;
    return value;
}

static inline
MortonIndex
PrepareValueFor3DInterleaving(MortonIndex value)
{
    /** Prepare the value for interleaving by fragmenting the bits, like:
     * (...b4 b3 b2 b1 b0) -> (...0 0 b4 0 0 b3 0 0 b2 0 0 b1 0 0 b0) */
    value = (value | (value << 14)  | (value << 28)) & 0x0001FC000FE0007F;
    value = (value | (value << 10)) & 0x06007C3003E1801F;
    value = (value | ((value << 6)  & ~(0x00000C0000600000))) & 0x06181C30C0E18607;
    value = (value | (value << 2))  & 0x1248649243249219;
    value = (value | (value << 2))  & 0x1249249249249249;
    return value;
}

static
MortonIndex
GetMortonIndex2D(const std::vector<DomainIndex>& coordinates)
{
    cmc_assert(coordinates.size() == 2);

    cmc_assert(static_cast<MortonIndex>(coordinates[0]) < (static_cast<MortonIndex>(1) << 32) &&
               static_cast<MortonIndex>(coordinates[1]) < (static_cast<MortonIndex>(1) << 32));

    const MortonIndex x = PrepareValueFor2DInterleaving(static_cast<MortonIndex>(coordinates[0]));
    const MortonIndex y = PrepareValueFor2DInterleaving(static_cast<MortonIndex>(coordinates[1]));

    /* Interleave both coordinates */
    return x | (y << 1);
}

static
MortonIndex
GetMortonIndex3D(const std::vector<DomainIndex>& coordinates)
{
    cmc_assert(coordinates.size() == 3);

    cmc_assert(static_cast<MortonIndex>(coordinates[0]) < (static_cast<MortonIndex>(1) << 21) &&
               static_cast<MortonIndex>(coordinates[1]) < (static_cast<MortonIndex>(1) << 21) &&
               static_cast<MortonIndex>(coordinates[2]) < (static_cast<MortonIndex>(1) << 21));

    const MortonIndex x = PrepareValueFor3DInterleaving(static_cast<MortonIndex>(coordinates[0]));
    const MortonIndex y = PrepareValueFor3DInterleaving(static_cast<MortonIndex>(coordinates[1]));
    const MortonIndex z = PrepareValueFor3DInterleaving(static_cast<MortonIndex>(coordinates[2]));

    /* Interleave all three coordinates */
    return x | (y << 1) | (z << 2);
}


MortonIndex
GetMortonIndex(const std::vector<DomainIndex>& coordinates, const int dimensionality)
{
    cmc_assert(coordinates.size() >= 2 && coordinates.size() <=3);

    switch(dimensionality)
    {
        case 2:
            return GetMortonIndex2D(coordinates);
        break;
        case 3:
            return GetMortonIndex3D(coordinates);
        break;
        default:
            cmc_err_msg("The dimensionality (", dimensionality, "D) is not considered for the computation of Morton indices.");
            return static_cast<MortonIndex>(CMC_ERR);
    }
}


std::vector<MortonIndex>
TransformHyperslabCoordinatesToMortonIndices(const std::vector<Hyperslab>& hyperslabs, const DataLayout layout, const GeoDomain& global_domain)
{
    /* Obtain the overall (local) number of coordinates */
    DomainIndex num_coordinates = 0;
    for (auto iter = hyperslabs.begin(); iter != hyperslabs.end(); ++iter)
    {
        num_coordinates += iter->GetNumberCoordinates();
    }

    std::vector<MortonIndex> morton_indices;
    morton_indices.reserve(num_coordinates);

    const int dimensionality = GetDimensionalityOfDataLayout(layout);

    std::vector<DomainIndex> linear_dimension_indices(dimensionality, 0);

    const UpdateHyperslabCoordinateFn HsUpdateFn = GetHyperslabCoordinatesIterationFunction(layout);

    /* Transform each hyperslab coordinate to a Morton index */
    for (auto iter = hyperslabs.begin(); iter != hyperslabs.end(); ++iter)
    {
        /* The global domain may not start with zeros as start indices for the domain. Therefore, we need to create a normalized hyperslab
        starting (globally) at (0,0,0,0), otherwise the Morton index calculation is shifted */
        const Hyperslab normalized_hyperslab = SubtractGeoDomainOffset(*iter, global_domain);

        const DomainIndex num_hyperslab_coordinates = normalized_hyperslab.GetNumberCoordinates();

        for (DomainIndex hs_coord_index = 0; hs_coord_index < num_hyperslab_coordinates; ++hs_coord_index)
        {
            HsUpdateFn(normalized_hyperslab, linear_dimension_indices, hs_coord_index);

            morton_indices.push_back(GetMortonIndex(linear_dimension_indices, dimensionality));
        }
    }

    return morton_indices;
}


}
