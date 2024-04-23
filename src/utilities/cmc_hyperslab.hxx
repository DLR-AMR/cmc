#ifndef CMC_HYPERSLAB_HXX
#define CMC_HYPERSLAB_HXX

#include "utilities/cmc_dimension_interval.hxx"

#include <array>
#include <vector>
#include <functional>

namespace cmc
{

typedef int64_t HyperslabIndex;

class Hyperslab
{
public:
    using iterator = std::array<HyperslabIndex, Dimension::NumCoordinates>::iterator;
    using const_iterator = std::array<HyperslabIndex, Dimension::NumCoordinates>::const_iterator;

    Hyperslab() = default;
    /* 1D Constructor */
    Hyperslab(const DimensionInterval& dimension1);
    /* 2D Constructor */
    Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2);
    /* 3D Constructor */
    Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
              const DimensionInterval& dimension3);
    /* 4D constructor */
    Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
              const DimensionInterval& dimension3, const DimensionInterval& dimension4);
    Hyperslab(const HyperslabIndex lon_start, const HyperslabIndex lon_count,
              const HyperslabIndex lat_start, const HyperslabIndex lat_count,
              const HyperslabIndex lev_start, const HyperslabIndex lev_count,
              const HyperslabIndex time_start, const HyperslabIndex time_count);
    ~Hyperslab(){};

    Hyperslab(const Hyperslab& other) = default;
    Hyperslab& operator=(const Hyperslab& other) = default;
    Hyperslab(Hyperslab&& other) = default;
    Hyperslab& operator=(Hyperslab&& other) = default;

    HyperslabIndex GetNumberCoordinates() const;

    HyperslabIndex GetNumberCoordinatesWithoutCertainDimension(const Dimension excluded_dimension) const;

    HyperslabIndex GetDimensionStart(const Dimension dimension) const;

    HyperslabIndex GetDimensionLength(const Dimension dimension) const;

    int GetDimensionality() const;

    iterator StartIndicesBegin() { return start_indices_.begin(); };
    iterator StartIndicesEnd() { return start_indices_.end(); };
    const_iterator StartIndicesBegin() const { return start_indices_.begin(); };
    const_iterator StartIndicesEnd() const { return start_indices_.end(); };
    const_iterator StartIndicesCBegin() const { return start_indices_.cbegin(); };
    const_iterator StartIndicesCEnd() const { return start_indices_.cend(); };

    iterator CountIndicesBegin() { return count_indices_.begin(); };
    iterator CountIndicesEnd() { return count_indices_.end(); };
    const_iterator CountIndicesBegin() const { return count_indices_.begin(); };
    const_iterator CountIndicesEnd() const { return count_indices_.end(); };
    const_iterator CountIndicesCBegin() const { return count_indices_.cbegin(); };
    const_iterator CountIndicesCEnd() const { return count_indices_.cend(); };

private:
    void SetupDimension(const DimensionInterval& dimension);
    
    std::array<HyperslabIndex, Dimension::NumCoordinates> start_indices_{0,0,0,0};
    std::array<HyperslabIndex, Dimension::NumCoordinates> count_indices_{1,1,1,1};
};

using UpdateHyperslabCoordinateFn = std::function<void(const Hyperslab&, std::vector<HyperslabIndex>&, const HyperslabIndex)>;

UpdateHyperslabCoordinateFn
GetHyperslabCoordinatesIterationFunction(const DataLayout layout);

void UpdateHyperslabIndexLonLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLatLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLatLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLevLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLonLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLevLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLonLatLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLevLonLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLonLevLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLevLatLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLatLevLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexLatLonLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);

void UpdateHyperslabIndexError(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index);


using DimensionValueExtractionFn = std::function<HyperslabIndex(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension)>;

DimensionValueExtractionFn
GetDimensionValueFunctionForReferenceCoords(const DataLayout layout);

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLonLat(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension);

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLatLev(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension);

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLonLev(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension);

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLonLatLev(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension);

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsError(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension);

using TrimRefCoordVectorFn = std::function<std::vector<HyperslabIndex>(const std::vector<HyperslabIndex>& reference_coordinates)>;

TrimRefCoordVectorFn
GetReferenceCoordTrimmingFunction(const DataLayout trimmed_data_layout);

std::vector<HyperslabIndex> TrimRefCoordsToLonLat(const std::vector<HyperslabIndex>& reference_coordinates);

std::vector<HyperslabIndex> TrimRefCoordsToLatLev(const std::vector<HyperslabIndex>& reference_coordinates);

std::vector<HyperslabIndex> TrimRefCoordsToLonLev(const std::vector<HyperslabIndex>& reference_coordinates);

std::vector<HyperslabIndex> TrimRefCoordsError(const std::vector<HyperslabIndex>& reference_coordinates);

}

#endif /* !CMC_HYPERSLAB_HXX */
