#ifndef CMC_GEO_DOMAIN_HXX
#define CMC_GEO_DOMAIN_HXX

#include "utilities/cmc_dimension_interval.hxx"

#include <array>

namespace cmc
{

class GeoDomain
{
public:
    GeoDomain() = default;
    /* 1D Constructor */
    GeoDomain(const DimensionInterval& dimension1);
    /* 2D Constructor */
    GeoDomain(const DimensionInterval& dimension1, const DimensionInterval& dimension2);
    /* 3D Constructor */
    GeoDomain(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
              const DimensionInterval& dimension3);
    /* 4D constructor */
    GeoDomain(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
              const DimensionInterval& dimension3, const DimensionInterval& dimension4);
    GeoDomain(const DomainIndex lon_dimension_start, const DomainIndex lon_dimension_end,
              const DomainIndex lat_dimension_start, const DomainIndex lat_dimension_end,
              const DomainIndex lev_dimension_start, const DomainIndex lev_dimension_end,
              const DomainIndex time_dimension_start, const DomainIndex time_dimension_end);

    GeoDomain(const GeoDomain& other) = default;
    GeoDomain& operator=(const GeoDomain& other) = default;
    GeoDomain(GeoDomain&& other) = default;
    GeoDomain& operator=(GeoDomain&& other) = default;

    ~GeoDomain() = default;

    DomainIndex GetDimensionStartIndex(const Dimension dimension) const;
    DomainIndex GetDimensionEndIndex(const Dimension dimension) const;

    DomainIndex GetDimensionLength(const Dimension dimension) const;
    DomainIndex GetLargestDimensionLength() const;
    DomainIndex GetNumberReferenceCoordsCovered() const;
    
    int GetDimensionality() const;

    void UpdateDimension(const DimensionInterval& dimension);

    void ClearDimension(const Dimension removed_dimension);

    friend bool CompareGeoDomains(const GeoDomain& domain1, const GeoDomain& domain2);
    friend GeoDomain ExtendGeoDomain(const GeoDomain& domain, const DimensionInterval& add_dimension);
    
    GeoDomain GetZeroOffsetDomain() const;
    
    bool IsValid() const;
private:
    void SetupDimension(const DimensionInterval& dimension);

    std::array<DomainIndex, Dimension::NumCoordinates> start_indices_{0,0,0,0};
    std::array<DomainIndex, Dimension::NumCoordinates> end_indices_{0,0,0,0};
};

}

#endif /* !CMC_GEO_DOMAIN_HXX */
