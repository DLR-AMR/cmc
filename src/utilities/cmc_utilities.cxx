#include "utilities/cmc_utilities.hxx"


inline int
GeoDomain::GetDimensionality() const
{
    return static_cast<int>(domain_.size());
}

inline void
GeoDomain::AddDimension(const Dimension dimension, const int start_idx, const int end_idx)
{
    if (GetDimensionality() < Dimension::NumCoordinates)
    {
        domain_.emplace_back(DimensionInterval(dimension, start_idx, end_idx));
    } else
    {
        cmc_err_msg("The GeoDomain is already fully defined in every dimension.");
    }
}

constexpr GeoDomain::GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1)
: domain_{DimensionInterval(dimension1, start_idx_dim1, end_idx_dim1)}
{};
constexpr GeoDomain::GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1,
                      const Dimension dimension2, const int start_idx_dim2, const int end_idx_dim2)
: domain_{DimensionInterval(dimension1, start_idx_dim1, end_idx_dim1),
          DimensionInterval(dimension2, start_idx_dim2, end_idx_dim2)}
{};
constexpr GeoDomain::GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1,
                     const Dimension dimension2, const int start_idx_dim2, const int end_idx_dim2,
                     const Dimension dimension3, const int start_idx_dim3, const int end_idx_dim3)
: domain_{DimensionInterval(dimension1, start_idx_dim1, end_idx_dim1),
          DimensionInterval(dimension2, start_idx_dim2, end_idx_dim2),
          DimensionInterval(dimension3, start_idx_dim3, end_idx_dim3)}
{};
constexpr GeoDomain::GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1,
                     const Dimension dimension2, const int start_idx_dim2, const int end_idx_dim2,
                     const Dimension dimension3, const int start_idx_dim3, const int end_idx_dim3,
                     const Dimension dimension4, const int start_idx_dim4, const int end_idx_dim4)
: domain_{DimensionInterval(dimension1, start_idx_dim1, end_idx_dim1),
          DimensionInterval(dimension2, start_idx_dim2, end_idx_dim2),
          DimensionInterval(dimension3, start_idx_dim3, end_idx_dim3),
          DimensionInterval(dimension4, start_idx_dim4, end_idx_dim4)}
{};
