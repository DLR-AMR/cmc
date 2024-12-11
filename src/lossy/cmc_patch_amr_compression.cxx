#include "lossy/cmc_patch_amr_compression.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_prefix.hxx"
#include "utilities/cmc_span.hxx"
#include "utilities/cmc_huffman.hxx"

#include <cstdint>
#include <variant>
#include <utility>
#include <type_traits>
#include <array>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <cmath>
#include <bitset>
/* Meshless AMR Compression in a block-based fashion */

/* Runtime optimzation may be to precompute the Morton indices reordering as a look up table for the fixed refinement levels
 *  as well as a reduction of single copy iterations in the hyperslab data extraction/indices computation or to re-formulate this all as a single function */
namespace cmc
{
    
namespace patch_amr
{

using LevelInteger = int_fast8_t;

constexpr LevelInteger kInvalidLevel = -1;

constexpr int kPatchAMRError = -1;

//Sollte auch einigermaßen einfach auf 1D und 4D erweitert werden können
constexpr int kInitialCompressionBlockLevel2D = 4;
constexpr int kInitialCompressionBlockLevel3D = 3;

static inline
int GetInitialRefinementLevel(const int dimensionality)
{
    if (dimensionality == 2)
    {
        return kInitialCompressionBlockLevel2D;
    } else if (dimensionality == 3)
    {
        return kInitialCompressionBlockLevel3D;
    } else
    {
        cmc_err_msg("PatchError in Level init");
        return kPatchAMRError;
    }
}

constexpr
int GetCompressionBlockSize(const int dimensionality)
{
    int exp = 0;
    if (dimensionality == 2)
    {
        exp = 2 * kInitialCompressionBlockLevel2D;
    } else if (dimensionality == 3)
    {
        exp = 3 * kInitialCompressionBlockLevel3D;
    }
    
    int size = 1;
    for (int exp_iter = 0; exp_iter < exp; ++exp_iter)
    {
        size *= 2;
    }

    return size;
}

constexpr
int GetCompressionBlockSizePerDimension(const int dimensionality)
{
    int exp = 0;
    if (dimensionality == 2)
    {
        exp = kInitialCompressionBlockLevel2D;
    } else if (dimensionality == 3)
    {
        exp = kInitialCompressionBlockLevel3D;
    }
    
    int size = 1;
    for (int exp_iter = 0; exp_iter < exp; ++exp_iter)
    {
        size *= 2;
    }

    return size;
}

constexpr int kCompressionBlockSize2D = GetCompressionBlockSize(2);
constexpr int kCompressionBlockSize3D = GetCompressionBlockSize(3);
constexpr int kMaxCompressionBlockSize = (kCompressionBlockSize2D >= kCompressionBlockSize3D ? kCompressionBlockSize2D : kCompressionBlockSize3D);

inline bool
IsCoarseningProgressing(const size_t previous_num_elems, const size_t num_elems)
{
    return (num_elems < previous_num_elems) && (num_elems > 1);
}

constexpr inline int
GetNumSiblings(const int dimensionality)
{
    switch (dimensionality)
    {
        case 1:
            return 2;
            break;
        case 2:
            return 4;
            break;
        case 3:
            return 8;
            break;
        case 4:
            return 16;
            break;
        default:
            cmc_err_msg("The supplied dimension is not supported (only 1D to 4D).");
            return 0;
    }
}

template<typename T>
void UpdateVectors(std::vector<T>& vec_old, std::vector<T>& vec_new)
{
    std::swap(vec_old, vec_new);
    vec_new.clear();
}

struct PermittedError
{
    PermittedError() = delete;
    PermittedError(const CompressionCriterion etype, const double permitted_error)
    : criterion{etype}, error{permitted_error}{};

    const CompressionCriterion criterion{CompressionCriterion::CriterionUndefined};
    const double error{0.0};
};


PermittedError
GetPermittedError()
{
    //TODO: Error is currently set here
    //return PermittedError(CompressionCriterion::AbsoluteErrorThreshold, 0.0625);
    return PermittedError(CompressionCriterion::RelativeErrorThreshold, 0.01);
}

template <typename T>
T
InterpolateValues(const std::vector<T>& values, const T missing_value)
{
    #if 1
    /* Calculate the mid-range */
    T max{std::numeric_limits<T>::min()};
    T min{std::numeric_limits<T>::max()};

    int num_elems = 0;
    for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
    {
        if (not ApproxCompare(*val_iter, missing_value))
        {
            if (max < *val_iter){max = *val_iter;}
            if (min > *val_iter){min = *val_iter;}
            ++num_elems;
        }
    }

    if (num_elems > 0)
    {
        return (max + min) / 2;
    } else
    {
        return missing_value;
    }
    #else
    /* Arithmetic mean */
    double sum{0.0};

    double num_elems = 0;
    for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
    {
        if (not ApproxCompare(*val_iter, missing_value))
        {
            sum += static_cast<double>(*val_iter);
            ++num_elems;
        }
    }

    if (num_elems > 0)
    {
        return static_cast<T>(sum / num_elems);
    } else
    {
        return missing_value;
    }
    #endif
}

static size_t
GetInitialValueCoverage(const size_t num_children, const size_t initial_refinement_level, const size_t level)
{
    cmc_assert(num_children != 0);
    cmc_assert(initial_refinement_level >= level);

    const size_t level_diff = initial_refinement_level - level;
    size_t coverage{1};

    for (size_t exp = 0; exp < level_diff; ++exp)
    {
        coverage *= num_children;
    }

    return coverage;
}


template<typename T>
inline auto
ComputeAbsoluteError(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value));
}

template<typename T>
inline auto
ComputeAbsoluteError(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value));
}


template<typename T>
inline auto
ComputeRelativeError(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the relative deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value)) / static_cast<double>(std::abs(initial_value));
}

template<typename T>
inline auto
ComputeRelativeError(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{

    /* Calculate the relative deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value)) / static_cast<double>(initial_value);
}

template<typename T>
auto
IsAbsoluteErrorCompliant(const T initial_value, const T nominal_value, const double error)
 -> std::enable_if_t<std::is_signed_v<T>, bool>
{
    /* Calculate the absolute deviation and check if it complies with the given error */
    return (ComputeAbsoluteError(initial_value, nominal_value) <= error);
}

template<typename T>
auto
IsAbsoluteErrorCompliant(const T initial_value, const T nominal_value, const double error)
 -> std::enable_if_t<std::is_unsigned_v<T>, bool>
{
    /* Calculate the absolute deviation and check if it complies with the given error */
    return (ComputeAbsoluteError(initial_value, nominal_value) <= error);
}


template<typename T>
auto
IsRelativeErrorCompliant(const T initial_value, const T nominal_value, const double error)
 -> std::enable_if_t<std::is_signed_v<T>, bool>
{
    /* Calculate the relative deviation and check if it complies with the given error */
    return (ComputeRelativeError(initial_value, nominal_value) <= error);
}

template<typename T>
auto
IsRelativeErrorCompliant(const T initial_value, const T nominal_value, const double error)
 -> std::enable_if_t<std::is_unsigned_v<T>, bool>
{

    /* Calculate the relative deviation and check if it complies with the given error */
    return (ComputeRelativeError(initial_value, nominal_value) <= error);
}


template <typename T>
bool
IsErrorCompliant(const VectorView<T>& initial_values, const T interpolated_value, const PermittedError& error_threshold, const T missing_value)
{
    //int i = 0;
    /* Iterate over all initial values angd check whether the interpolated value is error compliant */
    for (auto initial_val_iter = initial_values.begin(); initial_val_iter != initial_values.end(); ++initial_val_iter)
    {
        //cmc_debug_msg("Is error compliant step: ", i);
        //++i;

        if (ApproxCompare(*initial_val_iter, missing_value))
        {
            continue;
        }

        if (error_threshold.criterion == CompressionCriterion::AbsoluteErrorThreshold)
        {
            if (not IsAbsoluteErrorCompliant(*initial_val_iter, interpolated_value, error_threshold.error))
            {
                return false;
            }
        } else if (error_threshold.criterion == CompressionCriterion::RelativeErrorThreshold)
        {
            if (not IsRelativeErrorCompliant(*initial_val_iter, interpolated_value, error_threshold.error))
            {
                return false;
            }
        } else
        {
            cmc_err_msg("An undefined error criterion has been supplied.");
            return false;
        }
    }

    /* If all values comply, we can indicate that the coarsening is error compliant */
    return true;
}

template <typename T>
static void
LeaveValueUnchanged(const std::vector<T>& values, const std::vector<LevelInteger>& levels, std::vector<T>& data_new, std::vector<LevelInteger>& data_levels_new)
{
    data_new.push_back(values.front());
    data_levels_new.push_back(levels.front());
}

template <typename T>
static void
EvaluateAndApplyCoarsening(const std::vector<T>& values, const std::vector<LevelInteger>& levels, size_t& index, const VectorView<T>& initial_values,
                           std::vector<T>& data_new, std::vector<LevelInteger>& data_levels_new, const PermittedError& error_threshold,
                           const T missing_value)
{
    //cmc_assert(not values.empty());
    //cmc_assert(not levels.empty());
    //cmc_debug_msg("in evaluate and apply coarsening");
    /* Interpolate the values */
    const T interpolation_result = InterpolateValues(values, missing_value);
     //cmc_debug_msg("after interpolation");
    const bool is_error_compliant = IsErrorCompliant(initial_values, interpolation_result, error_threshold, missing_value);
     //cmc_debug_msg("After error compl check");

    if (is_error_compliant)
    {
        /* Store the coarsening results */
        data_new.push_back(interpolation_result);
        data_levels_new.push_back(levels.front() - 1);
    } else
    {
        /* The coarsening cannot be applied */
        (void) std::copy_n(values.begin(), values.size(), std::back_inserter(data_new));
        (void) std::copy_n(levels.begin(), levels.size(), std::back_inserter(data_levels_new));
    }
}

//TODO: Levels would not need to be a vector, a scalar would be sufficient
template <typename T>
static bool
GetNextValuesToEvaluate(const LevelInteger current_reference_level, std::vector<T>& values, std::vector<LevelInteger>& levels, const size_t index, const std::vector<T>& data, const std::vector<LevelInteger>& data_levels)
{
    cmc_assert(not values.empty());
    cmc_assert(not levels.empty());
    cmc_assert(index < data.size());

    const size_t end = std::min(index + values.size(), data.size());

    /* We only get values that are considered on the current reference level, All other values cannot be considered, since they
     * may have been already evaluated negative in previous iterations */
    for (size_t iter = index; iter < end; ++iter)
    {
        if (current_reference_level != data_levels[iter])
        {
            values[0] = data[index];
            levels[0] = data_levels[index];
            return false;
        }
    }

    const bool is_incomplete_family = (index + values.size() > data.size());

    if (not is_incomplete_family)
    {
        /* If we do have a family of elements */
        (void) std::copy_n(&data[index], values.size(), values.begin());
        (void) std::fill_n(levels.begin(), levels.size(), current_reference_level);
        return true;
    } else
    {
        /* If we do not have a family of elements */
        values[0] = data[index];
        levels[0] = data_levels[index];
        return false;
    }
}

template<typename T>
inline
auto int_pow(const T& base, const T& exp)
 -> std::enable_if_t<std::is_integral_v<T>, T>
{
    T result{static_cast<T>(1)};

    for (T i = 0; i < exp; ++i)
    {
        result *= base;
    }

    return result;
}

//TODO: Levels would not need to be a vector, a scalar would be sufficient
template <typename T>
static bool
GetNextValuesToEvaluateForCoarsening(const size_t initial_ref_lvl_index, const LevelInteger initial_reference_level, const LevelInteger current_reference_level, std::vector<T>& values, std::vector<LevelInteger>& levels, const size_t index, const std::vector<T>& data, const std::vector<LevelInteger>& data_levels)
{
    cmc_assert(initial_reference_level >= current_reference_level);
    cmc_assert(not values.empty());
    cmc_assert(not levels.empty());
    cmc_assert(index < data.size());

    const size_t num_siblings = values.size();
    const size_t end = std::min(index + num_siblings, data.size());

    /* We check whether the index corresponds to the start of a famiy, if not, we need to skip */
    //cmc_debug_msg("Initial RefLevel Index: ", initial_ref_lvl_index, " und int_pow: ", (int_pow(num_siblings, static_cast<size_t>(initial_reference_level - current_reference_level) + 1)), " mit resultat: ", (initial_ref_lvl_index % (int_pow(num_siblings, static_cast<size_t>(initial_reference_level - current_reference_level) + 1))));
    if (not (initial_ref_lvl_index % (int_pow(num_siblings, static_cast<size_t>(initial_reference_level - current_reference_level) + 1)) == 0))
    {
        values[0] = data[index];
        levels[0] = data_levels[index];
        return false;
    }

    /* We only get values that are considered on the current reference level, All other values cannot be considered, since they
     * may have been already evaluated negative in previous iterations */
    for (size_t iter = index; iter < end; ++iter)
    {
        if (current_reference_level != data_levels[iter])
        {
            values[0] = data[index];
            levels[0] = data_levels[index];
            return false;
        }
    }

    const bool is_incomplete_family = (index + values.size() > data.size());

    if (not is_incomplete_family)
    {
        /* If we do have a family of elements */
        (void) std::copy_n(&data[index], values.size(), values.begin());
        (void) std::fill_n(levels.begin(), levels.size(), current_reference_level);
        return true;
    } else
    {
        /* If we do not have a family of elements */
        values[0] = data[index];
        levels[0] = data_levels[index];
        return false;
    }
}

template <typename T>
std::pair<std::vector<LevelInteger>, std::vector<T>>
CoarsenValuesAdaptively(const InputVariable<T>& var, CompressionBlock<T>& compression_block, const PermittedError& error_thresholds)
{
    //cmc_debug_msg("In coarsen values adaptively");
    if (not compression_block.IsMortonOrdered())
    {
        compression_block.OrderMortonCurveCompliant();
    }

    const int dimensionality = compression_block.GetDimensionality();
    const size_t num_siblings = GetNumSiblings(dimensionality);
    const size_t initial_refinement_level = static_cast<size_t>(GetInitialRefinementLevel(dimensionality));
    
    const size_t num_elems = compression_block.size();

    //cmc_debug_msg("Num elems = ", num_elems);
    std::vector<LevelInteger> data_levels(num_elems, initial_refinement_level);
    std::vector<LevelInteger> data_levels_new;
    data_levels_new.reserve(data_levels.size());

    //cmc_debug_msg("Num initial values: ", compression_block.GetInitialData().size());
    const std::vector<T>& initial_values = compression_block.GetInitialData();
    std::vector<T> data = initial_values;
    std::vector<T> data_new;
    data_new.reserve(data.size());

    /* The level of data which considered for coarsening */
    LevelInteger current_reference_level = static_cast<LevelInteger>(initial_refinement_level);

    std::vector<T> values(num_siblings, static_cast<T>(0));
    std::vector<LevelInteger> levels(num_siblings, current_reference_level);

    //PermittedError error_thresholds = GetPermittedError(); //TODO: implement

    size_t previous_num_elems = data.size() + 1;
    //cmc_debug_msg("Size of data: ", data.size());
    while (IsCoarseningProgressing(previous_num_elems, data.size()))
    {
        cmc_assert(data_new.empty());
        cmc_assert(data_levels_new.empty());

        previous_num_elems = data.size();
        //cmc_debug_msg("Pt 1");
        size_t initial_values_index = 0;
        //size_t uniform_reference_level_id = 0;

        /* Iterate over all values and try to coarsen all element families */
        for (size_t index = 0; index < data.size();)
        {
            //cmc_debug_msg("Pt 2, index: ", index);
            const bool is_family = GetNextValuesToEvaluateForCoarsening(initial_values_index, initial_refinement_level, current_reference_level, values, levels, index, data, data_levels);
            
            if (is_family)
            {
                //cmc_debug_msg("Pt 3 is family");
                const size_t size_intial_value_coverage = GetInitialValueCoverage(num_siblings, initial_refinement_level, static_cast<size_t>(levels.front() - 1));
                //cmc_debug_msg("Pt 3 fam b");
                const VectorView<T> initial_values_view(initial_values.data() + initial_values_index, size_intial_value_coverage);
                //cmc_debug_msg("Pt 3 fam c");
                EvaluateAndApplyCoarsening(values, levels, index, initial_values_view, data_new, data_levels_new, error_thresholds, var.GetMissingValue());
                //cmc_debug_msg("Pt 3 fam d");
                initial_values_index += size_intial_value_coverage;
                //uniform_reference_level_id += size_intial_value_coverage;
                index += num_siblings;
            } else
            {
                //cmc_debug_msg("Pt 3 is no family ");
                const size_t size_intial_value_coverage = GetInitialValueCoverage(num_siblings, initial_refinement_level, static_cast<size_t>(levels.front()));
                LeaveValueUnchanged(values, levels, data_new, data_levels_new);
                initial_values_index += size_intial_value_coverage;
                //uniform_reference_level_id += size_intial_value_coverage;
                ++index;
            }

        }
        //cmc_debug_msg("Pt 4");
        /* Update the data and level vectors */
        UpdateVectors(data, data_new);
        UpdateVectors(data_levels, data_levels_new);

        //cmc_debug_msg("on ref level: ", current_reference_level, " hat data size: ", data.size());

        /* Update the reference level which will be considered in the next iteration */
        --current_reference_level;

        #if 0
        cmc_debug_msg("\n\nNach coarsening:");
        int ii = 0;
        for (auto viter = data_levels.begin(); viter != data_levels.end(); ++viter, ++ii)
        {
            cmc_debug_msg("Index: ", ii, ", level: ", static_cast<int>(*viter));
        }
        #endif
    }
    
    cmc_debug_msg("The block with ID ", compression_block.GetBlockID(), " has been coarsened from ", initial_values.size(), " values to ", data.size(), " values.");

    return std::make_pair(std::move(data_levels), std::move(data));
}

template<typename T>
std::vector<double>
GetErrorsAfterCoarsening(const CompressionBlock<T>& compression_block, const std::vector<LevelInteger>& levels, const std::vector<T>& coarsened_values, const T missing_value, const PermittedError& error_threshold)
{
    const int dimensionality = compression_block.GetDimensionality();
    const size_t num_siblings = GetNumSiblings(dimensionality);
    const size_t initial_reference_level = static_cast<size_t>(GetInitialRefinementLevel(dimensionality));

    const std::vector<T>& initial_values = compression_block.GetInitialData();
    
    std::vector<double> errors;
    errors.reserve(initial_values.size());
    
    size_t initial_val_index{0};
    size_t data_accessor{0};

    //cmc_debug_msg("\n\n Size of levels: ", levels.size(), "\n\n");
    for (auto level_iter = levels.begin(); level_iter != levels.end(); ++level_iter, ++data_accessor)
    {
        /* Get the end index of initial values that need to be checked */
        const size_t end_initial_index = initial_val_index + GetInitialValueCoverage(num_siblings, initial_reference_level, *level_iter);
        //cmc_debug_msg("Level: ", static_cast<int>(*level_iter), ", initial val index: ", initial_val_index, " und end index ", end_initial_index);
        /* Iterate over all values and check whether they comply with the error */
        for (; initial_val_index < end_initial_index; ++initial_val_index)
        {

            //if (initial_val_index >= initial_values.size()) {cmc_debug_msg("Initial val size: ", initial_values.size(), " und initial val index: ", initial_val_index);}
            if (ApproxCompare(initial_values[initial_val_index], missing_value))
            {
                /* If a missing value is present, we indicate that the maximum error has been reached in order to not compress it with supplementary methods */
                errors.push_back(std::numeric_limits<double>::max());
            }

            if (error_threshold.criterion == CompressionCriterion::AbsoluteErrorThreshold)
            {
                errors.push_back(ComputeAbsoluteError(initial_values[initial_val_index], coarsened_values[data_accessor]));
            } else if (error_threshold.criterion == CompressionCriterion::RelativeErrorThreshold)
            {
                errors.push_back(ComputeRelativeError(initial_values[initial_val_index], coarsened_values[data_accessor]));
            } else
            {
                cmc_err_msg("An undefined error criterion has been supplied.");
                errors.emplace_back();
            }

            //cmc_debug_msg("idx: ", initial_val_index, ": Inital val: ", initial_values[initial_val_index], ", coarsend val: ", coarsened_values[data_accessor], ", Error: ", errors.back());
            //if (std::isnan(errors.back()))
            //{
            //    cmc_warn_msg("Nan values in error at index: ", initial_val_index, ". Initial val was: ", initial_values[initial_val_index], " and coarsend val was: ", coarsened_values[data_accessor]);
            //}
        }
    }

    return errors;
}


template<typename T>
std::vector<double>
GetErrorsForCoarsenedValues(const CompressionBlock<T>& compression_block, const std::vector<LevelInteger>& levels, const std::vector<T>& coarsened_values, const T missing_value, const PermittedError& error_threshold)
{
    const int dimensionality = compression_block.GetDimensionality();
    const size_t num_siblings = GetNumSiblings(dimensionality);
    const size_t initial_reference_level = static_cast<size_t>(GetInitialRefinementLevel(dimensionality));

    const std::vector<T>& initial_values = compression_block.GetInitialData();
    
    std::vector<double> errors;
    errors.reserve(coarsened_values.size());
    
    size_t initial_val_index{0};
    size_t data_accessor{0};

    //cmc_debug_msg("\n\n Size of levels: ", levels.size(), "\n\n");
    for (auto level_iter = levels.begin(); level_iter != levels.end(); ++level_iter, ++data_accessor)
    {
        double max_err = 0.0;
        /* Get the end index of initial values that need to be checked */
        const size_t end_initial_index = initial_val_index + GetInitialValueCoverage(num_siblings, initial_reference_level, *level_iter);
        //cmc_debug_msg("Level: ", static_cast<int>(*level_iter), ", initial val index: ", initial_val_index, " und end index ", end_initial_index);
        /* Iterate over all values and check whether they comply with the error */
        for (; initial_val_index < end_initial_index; ++initial_val_index)
        {

            //if (initial_val_index >= initial_values.size()) {cmc_debug_msg("Initial val size: ", initial_values.size(), " und initial val index: ", initial_val_index);}
            if (ApproxCompare(initial_values[initial_val_index], missing_value))
            {
                continue;
            }

            if (error_threshold.criterion == CompressionCriterion::AbsoluteErrorThreshold)
            {
                const double abs_err = ComputeAbsoluteError(initial_values[initial_val_index], coarsened_values[data_accessor]);
                if (abs_err > max_err){max_err = abs_err;}

            } else if (error_threshold.criterion == CompressionCriterion::RelativeErrorThreshold)
            {
                const double rel_err = ComputeRelativeError(initial_values[initial_val_index], coarsened_values[data_accessor]);
                if (rel_err > max_err){max_err = rel_err;}
            } else
            {
                cmc_err_msg("An undefined error criterion has been supplied.");
            }

            //cmc_debug_msg("idx: ", initial_val_index, ": Inital val: ", initial_values[initial_val_index], ", coarsend val: ", coarsened_values[data_accessor], ", Error: ", errors.back());
            //if (std::isnan(errors.back()))
            //{
            //    cmc_warn_msg("Nan values in error at index: ", initial_val_index, ". Initial val was: ", initial_values[initial_val_index], " and coarsend val was: ", coarsened_values[data_accessor]);
            //}
        }

        errors.push_back(max_err);
    }

    return errors;
}



static
std::array<DomainIndex, Dimension::NumCoordinates>
GetNumBlocksPerDimension(const GeoDomain& global_domain, const DomainIndex block_size_per_dim)
{
    std::array<DomainIndex, Dimension::NumCoordinates> num_blocks_per_dim;
    num_blocks_per_dim.fill(1);

    for (int i = 0; i < Dimension::NumCoordinates; ++i)
    {
        const DomainIndex dim_length = global_domain.GetDimensionLength(static_cast<Dimension>(i));
        if (dim_length >= 1)
        {
            /* Determine the number of blocks needed in this direction */
            num_blocks_per_dim[i] = (dim_length / block_size_per_dim) + (dim_length % block_size_per_dim != 0 ? 1 : 0);
        }
    }

    return num_blocks_per_dim;
}

template<typename T>
std::vector<CompressionBlock<T>>
GetCompressionBlocks(const InputVariable<T>& var)
{
    const DataLayout layout = var.GetInitialDataLayout();
    const int dimensionality = GetDimensionalityOfDataLayout(layout);
    const GeoDomain& global_domain = var.GetGlobalDomain();

    const DomainIndex largest_dim = global_domain.GetLargestDimensionLength();
    const int initial_reference_level = GetInitialRefinementLevel(dimensionality);
    const int block_size = GetCompressionBlockSize(dimensionality);
    const int block_size_per_dim = GetCompressionBlockSizePerDimension(dimensionality);

    std::array<DomainIndex, Dimension::NumCoordinates> num_blocks_per_dim = GetNumBlocksPerDimension(global_domain, static_cast<DomainIndex>(block_size_per_dim));

    /* Gather the number of total blocks that need to be generated */
    const DomainIndex num_blocks = std::accumulate(num_blocks_per_dim.begin(), num_blocks_per_dim.end(), 1, std::multiplies<DomainIndex>());

    std::vector<CompressionBlock<T>> blocks;
    blocks.reserve(num_blocks);

    bool is_time_considered = (global_domain.GetDimensionLength(Dimension::Time) > 1 ? true : false);
    bool is_lev_considered = (global_domain.GetDimensionLength(Dimension::Lev) > 1 ? true : false);
    bool is_lat_considered = (global_domain.GetDimensionLength(Dimension::Lat) > 1 ? true : false);
    bool is_lon_considered = (global_domain.GetDimensionLength(Dimension::Lon) > 1 ? true : false);

    const DomainIndex time_offset = global_domain.GetDimensionStartIndex(Dimension::Time);
    const DomainIndex lev_offset = global_domain.GetDimensionStartIndex(Dimension::Lev);
    const DomainIndex lat_offset = global_domain.GetDimensionStartIndex(Dimension::Lat);
    const DomainIndex lon_offset = global_domain.GetDimensionStartIndex(Dimension::Lon);

    int block_id{0};

    /* Create all blocks */
    for (DomainIndex time_iter = 0; time_iter < num_blocks_per_dim[Dimension::Time]; ++time_iter)
    {
        for (DomainIndex lev_iter = 0; lev_iter < num_blocks_per_dim[Dimension::Lev]; ++lev_iter)
        {
            for (DomainIndex lat_iter = 0; lat_iter < num_blocks_per_dim[Dimension::Lat]; ++lat_iter)
            {
                for (DomainIndex lon_iter = 0; lon_iter < num_blocks_per_dim[Dimension::Lon]; ++lon_iter, ++block_id)
                {
                    Hyperslab block(
                        (is_lon_considered ? lon_offset + lon_iter * block_size_per_dim : 0), (is_lon_considered ? block_size_per_dim : 1),
                        (is_lat_considered ? lat_offset + lat_iter * block_size_per_dim : 0), (is_lat_considered ? block_size_per_dim : 1),
                        (is_lev_considered ? lev_offset + lev_iter * block_size_per_dim : 0), (is_lev_considered ? block_size_per_dim : 1),
                        (is_time_considered ? time_offset + time_iter * block_size_per_dim : 0), (is_time_considered ? block_size_per_dim : 1)
                    );

                    /* Create a new compression block */
                    blocks.emplace_back(block_id, std::move(block), layout);
                }
            }
        }
    }

    /* After all blocks have been defined, we need to get the blocks' corresponding data */
    for (auto block_iter = blocks.begin(); block_iter != blocks.end(); ++block_iter)
    {
        block_iter->SetData(var.GetDataFromHyperslab(block_iter->GetHyperslab()));
    }

    return blocks;
}


template <typename T>
void
XorValues(std::vector<T>& values)
{
    // Idea should be: Build XOR values from the coarsened data try to erase as many bits
    // from the end position with the remaining (is this really possible;
    // only instantaneously when xoring otherwise it would change all succesive xor values?)

    //DO we try to alter first unequal bit in the suffixes with the remaining error (maybe not)

    //For data that is negative and positive we should remove the sign bit and store is seperately (rle encoded)

    int32_t val{0};
    int32_t xorval{0};

    for (size_t index = 1; index < values.size(); ++index)
    {
        std::memcpy(&val, static_cast<void*>(&values[index]), 4);
        cmc_debug_msg(std::bitset<32>(val), " fuer Elem: ", index, ", original value: ", values[index]);

        uint32_t v1{0}, v2{0};
        std::memcpy(&v1, static_cast<void*>(&values[index]), 4);
        std::memcpy(&v2, static_cast<void*>(&values[index - 1]), 4);

        uint32_t xval = (v1 ^ v2);
        cmc_debug_msg(std::bitset<32>(xval), " fuer Elem: ", index);
    }
}

template<typename T>
std::vector<CompressionValue<sizeof(T)>>
TransformCompressionToByteValues(std::vector<T>&& values)
{
    const int N = sizeof(T);

    std::vector<CompressionValue<N>> byte_values;
    byte_values.reserve(values.size());

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        byte_values.emplace_back(*iter);
    }
    
    return byte_values;
}

template<typename T>
struct PrefixData
{
    constexpr static int N = sizeof(T);

    PrefixData() = delete;
    PrefixData(const size_t size_hint_num_elems)
    {
        levels.reserve(size_hint_num_elems);
        byte_values.reserve(size_hint_num_elems);
        bit_indications.Reserve((2 * size_hint_num_elems) / bit_vector::kCharBit + 1);
        coarsening_indication.Reserve(size_hint_num_elems / bit_vector::kCharBit + 1);

        //cmc_debug_msg("In prefix data: N = ", N, " und sizof T: ", sizeof(T));
    };

    PrefixData(std::vector<LevelInteger>&& value_levels, std::vector<CompressionValue<N>>&& values)
    : levels{std::move(value_levels)}, byte_values(std::move(values)) {
        bit_indications.Reserve(2 * levels.size());
        coarsening_indication.Reserve(levels.size());
    };

    std::vector<LevelInteger> levels;
    std::vector<CompressionValue<N>> byte_values;
    bit_map::BitMap bit_indications;
    bit_map::BitMap coarsening_indication;
};

constexpr bool kIndicatePrefixHasBeenExtracted = true;
constexpr bool kIndicateNoPrefix = false;
constexpr bool kIndicateCoarsening = true;
constexpr bool kIndicateElementStaysUnchanged = false;


template<typename T>
class PrefixPyramid
{
public:
    constexpr static int N = sizeof(T);

    PrefixPyramid() = delete;
    PrefixPyramid(const int block_id, const Hyperslab& domain, std::vector<LevelInteger>&& levels, std::vector<T>&& values)
    : block_id_{block_id}, domain_{domain} {
        const auto max_element_iter = std::max_element(levels.begin(), levels.end());
        #ifdef CMC_ENABLE_DEBUG
        cmc_assert(max_element_iter != levels.end());
        //cmc_debug_msg("PyramidPrefix constructor with max level: ", static_cast<int>(*max_element_iter));
        //cmc_debug_msg("Size of t is: ", sizeof(T));
        #endif
        //for (auto iter = levels.begin(); iter != levels.end(); ++iter)
        //{
        //    cmc_debug_msg("levels: ", static_cast<int>(*iter));
        //}
        initial_max_level_ = (max_element_iter != levels.end() ? *max_element_iter : kInvalidLevel);
        //cmc_debug_msg("Initial max level is: ", static_cast<int>(initial_max_level_ + 1));
        levelwise_prefixes_.reserve(initial_max_level_ + 1);
        
        levelwise_prefixes_.emplace_back(std::move(levels), TransformCompressionToByteValues<T>(std::move(values)));
    };

    PrefixPyramid(const int block_id, const Hyperslab& domain, std::vector<LevelInteger>&& levels, std::vector<CompressionValue<sizeof(T)>>&& byte_values)
    : block_id_{block_id}, domain_{domain} {
        const auto max_element_iter = std::max_element(levels.begin(), levels.end());
        cmc_assert(max_element_iter != levels.end());
        initial_max_level_ = (max_element_iter != levels.end() ? *max_element_iter : kInvalidLevel);
        levelwise_prefixes_.reserve(initial_max_level_ + 1);
        levelwise_prefixes_.emplace_back(std::move(levels), std::move(byte_values));
    };

    void ExtractPrefixes();
    const std::vector<CompressionValue<N>>& GetInitialByteValues() const { return levelwise_prefixes_.front().byte_values;}
    const std::vector<LevelInteger>& GetInitialLevels() const {return levelwise_prefixes_.front().levels;}
    
    //std::vector<uint8_t> EncodeAndSerializeData() const;
    bit_vector::BitVector EncodeAndSerializeData() const;

private:

    LevelInteger GetInitialMaxReferenceLevel() const {return initial_max_level_;};
    int GetDimensionality() const {return domain_.GetDimensionality();}

    const std::vector<CompressionValue<N>>& GetPreviousByteValues() const {cmc_assert(levelwise_prefixes_.size() > 1); return (*std::prev(levelwise_prefixes_.end(), 2)).byte_values;}
    const std::vector<LevelInteger>& GetPreviousLevels() const {cmc_assert(levelwise_prefixes_.size() > 1); return (*std::prev(levelwise_prefixes_.end(), 2)).levels;}

    void EmplaceBackPrefixData(const size_t size_hint_num_elems) {levelwise_prefixes_.emplace_back(size_hint_num_elems);}

    void PushBackPrefix(CompressionValue<sizeof(T)>&& prefix) {levelwise_prefixes_.back().byte_values.push_back(std::move(prefix));}
    void PushBackPrefix(const CompressionValue<sizeof(T)>& prefix) {levelwise_prefixes_.back().byte_values.push_back(prefix);}
    void PushBackLevel(const LevelInteger level) {levelwise_prefixes_.back().levels.push_back(level);}
    void PushBackIndicationBit(const bool IsBitSet) {
        levelwise_prefixes_.back().bit_indications.AppendBit(IsBitSet);
    }

    void PushBackCoarseningBit(const bool IsBitSet) {
        PushBackIndicationBit(IsBitSet);
        levelwise_prefixes_.back().coarsening_indication.AppendBit(IsBitSet);
    }

    void PushBackPrefixIndicationBit(const bool IsBitSet)
    {
        PushBackIndicationBit(IsBitSet);
    }

    void TrimPreviousByteValues(const size_t index, const size_t num_siblings, const std::vector<size_t>& coarse_value_indices)
    {
        /* Get the previous byte values/prefixes (second to last) */
        std::vector<CompressionValue<sizeof(T)>>& prev_prefixes = (*std::prev(levelwise_prefixes_.end(), 2)).byte_values;

        const CompressionValue<sizeof(T)>& prefix = levelwise_prefixes_.back().byte_values.back();

        for (size_t elem_index = 0; elem_index < num_siblings; ++elem_index)
        {
            /* If the values have been coarsened from the beginning onwards, we just trim the previous prefix.
             * However, if the it is the first time a coarser value runs into the extarction, we will store the suffix on
             * suffix byte vector in the pyramid */
            bool is_coarsened_from_the_beginning_onwards{true};

            if (levelwise_prefixes_.size() > 2)
            {
                if (not (*std::prev(levelwise_prefixes_.end(), 2)).bit_indications.IsBitSet(2 * (index + elem_index) + 1))
                {
                    is_coarsened_from_the_beginning_onwards = false;
                }
            }

            if (is_coarsened_from_the_beginning_onwards)
            {
                prev_prefixes[index + elem_index].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());
                //offset_coarse_value_index += num_siblings;
                LevelInteger init_level = GetInitialCoarseValueLevel(coarse_value_indices[index + elem_index]);
                //cmc_assert(init_level == initial_max_level_);
                const LevelInteger current_ref_level = initial_max_level_ - (static_cast<LevelInteger>((levelwise_prefixes_.size())) - 2);
                //coarse_value_index_ += GetInitialValueCoverage(num_siblings, init_level, current_ref_level) + (current_ref_level == init_level ? 0 : 0);
                //cmc_debug_msg("Coarsened from beginning onwards: index: ", index, " und elem_index: ", elem_index, ", init_level: ", static_cast<int>(init_level), ", current_ref_level: ", static_cast<int>(current_ref_level));
            } else
            {
                //cmc_debug_msg("Reset in suffix array;", " coarse value index: ", coarse_value_indices[index + elem_index]);
                /* If it is the first time a coarser value is extracted, we assign an empty compression value to the previous prefixes and store
                 *  the suffix in the initial byte vector at the bottom of the prefix pxramid */
                levelwise_prefixes_.front().byte_values[coarse_value_indices[index + elem_index]] = prev_prefixes[index + elem_index];
                levelwise_prefixes_.front().byte_values[coarse_value_indices[index + elem_index]].SetFrontBit(sizeof(T) * CHAR_BIT - prefix.GetTrailBit());
                prev_prefixes[index + elem_index] = CompressionValue<sizeof(T)>();
                //++coarse_value_index_;
            }

            /* Check if the prefix is empty, if so we will erase the indicator bit */
            if (levelwise_prefixes_.size() > 2)
            {
            if (prev_prefixes[index + elem_index].IsEmpty())
            {
                bit_map::BitMap& indications = (*std::prev(levelwise_prefixes_.end(), 2)).bit_indications;
                //cmc_debug_msg("Clear Bit an pos: ", 2 * (index + elem_index) + 1, "indications size: ", indications.size());
                indications.ClearBit(2 * (index + elem_index) + 1);
                
                bit_map::BitMap& coarsening_indications = (*std::prev(levelwise_prefixes_.end(), 2)).coarsening_indication;
                coarsening_indications.ClearBit(index + elem_index);
            }
            }
        }
    }

    void ResetPreviousByteValue(const size_t index)
    {
        /* Get the previous byte values/prefixes (second to last) */
        std::vector<CompressionValue<sizeof(T)>>& prev_prefixes = (*std::prev(levelwise_prefixes_.end(), 2)).byte_values;

        /* Exchange the prefix with an empty one */
        prev_prefixes[index] = CompressionValue<sizeof(T)>();
    }

    LevelInteger GetInitialCoarseValueLevel(const size_t pos){cmc_assert(pos < levelwise_prefixes_.front().levels.size()); return levelwise_prefixes_.front().levels[pos];}

    size_t GetSizeOfCurrentByteValues() const  {return levelwise_prefixes_.back().byte_values.size();}

    size_t GetCurrentCoarseValueIndex() const {return coarse_value_index_;}
    void AddToInitialValuePos(const size_t summand) {coarse_value_index_ += summand;} 
    void ResetInitialValuePos() {coarse_value_index_ = 0;}

    std::vector<uint16_t> GetPrefixLengthFrequency() const;

    /* Member variables */
    int block_id_;
    Hyperslab domain_;
    LevelInteger initial_max_level_;
    std::vector<PrefixData<T>> levelwise_prefixes_;
    size_t coarse_value_index_{0};
};

inline bool
IsPrefixExtractionProgressing(const LevelInteger current_ref_level)
{
    return (current_ref_level > 0);
}

template<typename T>
std::pair<bool, CompressionValue<sizeof(T)>>
EvaluateCommonPrefix(const std::vector<CompressionValue<sizeof(T)>>& prefixes)
{
    /* Initialize the first prefix with the first value */
    CompressionValue<sizeof(T)> prefix = prefixes.front();

    if (prefix.IsEmpty()) {return std::make_pair(false, CompressionValue<sizeof(T)>());}

    /* Iterate over all leftover values */
    for (auto pref_iter = ++(prefixes.begin()); pref_iter != prefixes.end(); ++pref_iter)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix, *pref_iter);

        /* Check if there is a prefix */
        if (prefix.IsEmpty()) {return  std::make_pair(false, CompressionValue<sizeof(T)>());}
    }

    /* If the function arrives here, we do have found a common prefix of all values */
    return std::make_pair(true, prefix);
}

template<typename T>
void
PrefixPyramid<T>::ExtractPrefixes()
{
    const LevelInteger initial_max_refinement_level = GetInitialMaxReferenceLevel();

    const int dimensionality = GetDimensionality();
    const size_t num_siblings = GetNumSiblings(dimensionality);

    std::vector<LevelInteger> levels(num_siblings);
    std::vector<CompressionValue<sizeof(T)>> values(num_siblings);

    LevelInteger current_ref_level = initial_max_refinement_level;

    /* At the beginning the coarse value indices directly correspond to the vector indices of the coarsened values */
    std::vector<size_t> coarse_value_indices(GetSizeOfCurrentByteValues() + 1);
    std::iota(coarse_value_indices.begin(), coarse_value_indices.end(), 0);
    std::vector<size_t> coarse_value_indices_new;
    coarse_value_indices_new.reserve(coarse_value_indices.size());

    while (IsPrefixExtractionProgressing(current_ref_level))
    {
        int cv_idx_accessor = 0;
        /* We emplace a new PrefixData struct for this iteration of prefix extraction */
        /* We allocate a loose upper bound of memory */
        const size_t size_hint = GetSizeOfCurrentByteValues();
        EmplaceBackPrefixData(size_hint);

        /* Get the previous level and value data */
        const std::vector<LevelInteger>& previous_levels = GetPreviousLevels();
        const std::vector<CompressionValue<sizeof(T)>>& previous_byte_values = GetPreviousByteValues();

    #if 0
        int iiii = 0;
        size_t val_count = 0;
        for (auto liter = previous_levels.begin(); liter != previous_levels.end(); ++liter, ++iiii)
        {
            cmc_debug_msg("Index: ", iiii, ", hat level: ", static_cast<int>(*liter));
            val_count += GetInitialValueCoverage(num_siblings, 4, *liter);
        }
        cmc_debug_msg("Damit kommt man auf Num elemente: ", val_count);
    #endif

    
        coarse_value_indices_new.push_back(0);

        /* Iterate over all values and try to coarsen all element families */
        for (size_t index = 0; index < previous_byte_values.size();)
        {
            //cmc_debug_msg("Before Xetraction, index: ", index, " Coarse Value Pos: ", GetCurrentCoarseValueIndex());
            const bool is_family = GetNextValuesToEvaluate(current_ref_level, values, levels, index, previous_byte_values, previous_levels);
            
            if (is_family)
            {
                /* If it is a family, we are going to extract a prefix */
                auto [is_common_prefix_present, prefix] = EvaluateCommonPrefix<T>(values);

                if (is_common_prefix_present)
                {
                    //cmc_debug_msg("Prefix extraction takes place");
                    /* Store the prefix and set the indication bits */
                    PushBackLevel(levels.front() - 1);
                    PushBackPrefix(std::move(prefix));
                    PushBackCoarseningBit(kIndicateCoarsening);
                    PushBackPrefixIndicationBit(kIndicatePrefixHasBeenExtracted);

                    /* When we extract a prefix, we need to update the last byte values accordingly */
                    /* We need to trim the previous prefixes by the extracted common prefix */
                    TrimPreviousByteValues(index, num_siblings, coarse_value_indices);

                    /* Update the index by the process family of elements */
                    index += num_siblings;

                    coarse_value_indices_new.push_back(coarse_value_indices[cv_idx_accessor + num_siblings]);
                    cv_idx_accessor += num_siblings;

                    //for (size_t elem_id = 0; elem_id < num_siblings; ++elem_id)
                    //{
                    //    const LevelInteger initial_level = GetInitialCoarseValueLevel(GetCurrentCoarseValueIndex());
                    //    AddToInitialValuePos(GetInitialValueCoverage(num_siblings, initial_level, current_ref_level) + (initial_level == current_ref_level ? 0 : 1));
                    //}
                } else
                {
                    //cmc_debug_msg("NOOOOO Prefix extraction takes place");
                    /* In case we do not have a common prefix to extract */
                    PushBackLevel(levels.front() - 1);
                    PushBackPrefix(CompressionValue<sizeof(T)>());
                    PushBackCoarseningBit(kIndicateCoarsening);
                    PushBackPrefixIndicationBit(kIndicateNoPrefix);

                    /* Update the index by the processed family of elements */
                    index += num_siblings;

                    for (size_t elem_id = 0; elem_id < num_siblings; ++elem_id)
                    {
                        const LevelInteger initial_level = GetInitialCoarseValueLevel(GetCurrentCoarseValueIndex());
                        AddToInitialValuePos(GetInitialValueCoverage(num_siblings, initial_level, current_ref_level) + (initial_level == current_ref_level ? 0 : 0));
                    }

                    //coarse_value_indices[cv_idx_accessor + 1] = coarse_value_indices[cv_idx_accessor] + num_siblings;
                    //++cv_idx_accessor;

                    coarse_value_indices_new.push_back(coarse_value_indices[cv_idx_accessor + num_siblings]);
                    cv_idx_accessor += num_siblings;
                }
            } else
            {
                /* In case we do not have a family, the element stays unchanged and we drag the value along
                 * until it can be coarsened until the siblings are on the same refinement level */
                PushBackLevel(levels.front());
                PushBackPrefix(values.front());
                PushBackCoarseningBit(kIndicateElementStaysUnchanged);
                PushBackPrefixIndicationBit(kIndicateNoPrefix);
                
                /* Since we copied the whole value as a prefix, we need to erase it from the previous prefixes */
                if (current_ref_level < initial_max_refinement_level)
                {
                    ResetPreviousByteValue(index);
                }
                /* Update the index by a single element */
                ++index;

                AddToInitialValuePos(1);

                //coarse_value_indices[cv_idx_accessor + 1] = coarse_value_indices[cv_idx_accessor] + 1;
                //++cv_idx_accessor;

                coarse_value_indices_new.push_back(coarse_value_indices[cv_idx_accessor + 1]);
                ++cv_idx_accessor;
            }
            //cmc_debug_msg("After Xetraction, index: ", index, " Coarse Value Pos: ", GetCurrentCoarseValueIndex());
        }

        std::swap(coarse_value_indices, coarse_value_indices_new);
        coarse_value_indices_new.clear();

        //if (current_ref_level == initial_max_refinement_level)
        //{
        //    const std::vector<CompressionValue<sizeof(T)>>& extra_init_byte_vals = GetPreviousByteValues();
        //
        //    int ii = 0;
        //    for (auto viter = extra_init_byte_vals.begin(); viter != extra_init_byte_vals.end(); ++viter, ++ii)
        //    {
        //        cmc_debug_msg("Index: ", ii, ", Num Significant Bits: ", viter->GetCountOfSignificantBits());
        //    }
        //}
        //cmc_debug_msg("After level: ", current_ref_level, " Num Elements: ", GetCurrentByteValues().size());

        /* Update the reference level such that the extraction continues on the next coarser level */
        --current_ref_level;
        ResetInitialValuePos();

        //cmc_debug_msg("After extraction num data vals: ", GetSizeOfCurrentByteValues());
    }

    //const std::vector<CompressionValue<sizeof(T)>>& suffixes = GetInitialByteValues();
    //const std::vector<LevelInteger>& suff_level = GetInitialLevels();
    //int ii = 0;
    //for (auto viter = suffixes.begin(); viter != suffixes.end(); ++viter, ++ii)
    //{
    //    cmc_debug_msg("Index: ", ii, ", level: ", static_cast<int>(suff_level[ii]), " Num Significant Bits: ", viter->GetCountOfSignificantBits());
    //}
}

template<typename T>
std::vector<uint16_t>
PrefixPyramid<T>::GetPrefixLengthFrequency() const
{
    std::vector<uint16_t> frequency(sizeof(T) * bit_map::kCharBit, uint16_t{0});

    /* Iterate over all levels of extraction */
    for (auto lvl_iter = levelwise_prefixes_.begin(); lvl_iter != levelwise_prefixes_.end(); ++lvl_iter)
    {
        /* Iterate over all prefixes on this level */
        for (auto val_iter = lvl_iter->byte_values.begin(); val_iter != lvl_iter->byte_values.end(); ++val_iter)
        {
            /* Get the count of significant bits, i.e. the length of the prefix */
            const uint16_t pref_length = static_cast<uint16_t>(val_iter->GetCountOfSignificantBits());

            /* Increment the counter for the given prefix length */
            ++frequency[pref_length];
        }
    }

    return frequency;
}

#if 0
class ArithmeticCoder
{
public:
    ArithmeticCoder() = default;

private:

};

template<typename T>
ArithmeticCoder
PrefixPyramid<T>::CreateArithmeticCoder()
{
    std::vector<uint16_t> frequency = GetPrefixLengthFrequency();

    
}
#endif

#if 0
template<typename T>
std::vector<uint8_t>
PrefixPyramid<T>::EncodeAndSerializeData() const
{
    std::vector<uint16_t> frequency = GetPrefixLengthFrequency();

    huffman::HuffmanTree huffman_encoder(frequency);

    huffman_encoder.EncodeSymbols(std::vector<uint8_t>());


    return std::vector<uint8_t>();   
}
#endif

static size_t byte_counter = 0;

/* We encode the block data as follows: [huffman_tree], [refinement_bit, Encoded Prefix Length, Prefix], ... */
template<typename T>
bit_vector::BitVector
PrefixPyramid<T>::EncodeAndSerializeData() const
{
    std::vector<uint16_t> frequency = GetPrefixLengthFrequency();

    huffman::HuffmanTree<uint16_t> huffman_encoder(frequency);

    //huffman_encoder.EncodeSymbols(std::vector<uint8_t>());

    bit_vector::BitVector encoded_byte_stream;
    encoded_byte_stream.Reserve(levelwise_prefixes_.front().byte_values.size());

    /* Serialize the Huffman tree and append it in the encoded stream */
    const size_t symbol_size = 5;
    bit_vector::BitVector serialized_huffman_tree = huffman_encoder.SerializeHuffmanTree(symbol_size);
    byte_counter += serialized_huffman_tree.size_bits();

    //cmc_debug_msg("Bit Vector will be appended to serialized stream");
    encoded_byte_stream.AppendBits(serialized_huffman_tree);

    for (auto lvl_iter = levelwise_prefixes_.rbegin(); lvl_iter != levelwise_prefixes_.rend(); ++lvl_iter)
    {
        size_t index{0};

        for (auto val_iter = lvl_iter->byte_values.begin(); val_iter != lvl_iter->byte_values.end(); ++val_iter, ++index)
        {
            //cmc_debug_msg("We are on level: ", std::distance(levelwise_prefixes_.rbegin(), lvl_iter), " with index: ", index);

            /* For each prefix level we set the indication bit, on the suffix level, we do not have any indication bits */
            if (lvl_iter != std::prev(levelwise_prefixes_.rend()))
            {
                /* Check if a coarsening bit is set in the bitmap */
                const bool is_coarsening_bit_set = lvl_iter->coarsening_indication.IsBitSet(index);
                
                /* Set the coarsening/refinement indication in the stream */
                encoded_byte_stream.AppendBit(is_coarsening_bit_set);
            }


            /* Get the count of significant bits of the current prefix in order to evaluate whether there are prefix bits */
            const uint8_t num_bits_prefix = static_cast<uint8_t>(val_iter->GetCountOfSignificantBits());

            /* Encode the prefix length with the huffman tree */
            const auto [code, code_length] = huffman_encoder.EncodeSymbol(num_bits_prefix);
            if (code_length > 0)
            {
                byte_counter += code_length;
                /* The code length is only encoded when there is an actual code. A code of length zero occurs, when the whole block is 
                 * is coarsened to only a single value (the optimal code in that case is of length zero) */
                encoded_byte_stream.AppendBits(code, static_cast<int>(code_length));
            }

            if (num_bits_prefix != 0)
            {
                /* If there is a prefix to encode, we get the prefix and set it in the encoded stream */
                const std::vector<uint8_t> prefix = val_iter->GetSignificantBitsInBigEndianOrdering();
                encoded_byte_stream.AppendBits(prefix, num_bits_prefix);
            }
        }
    }

    /* We need to store the number of prefix levels as well in order to perform the correct iterations during decompression */
    const uint8_t num_prefix_levels = levelwise_prefixes_.size() - 1;
    //return std::make_pair(encoded_byte_stream, num_prefix_levels);
    //cmc_debug_msg("Accumulated Huffman tree bits: ", byte_counter);

    return encoded_byte_stream;
}


template<typename T>
void
ExtractPrefixes(PrefixPyramid<T>& prefix_pyramid)
{
    const LevelInteger initial_max_refinement_level = prefix_pyramid.GetInitialMaxReferenceLevel();

    const int dimensionality = prefix_pyramid.GetDimensionality();
    const size_t num_siblings = GetNumSiblings(dimensionality);

    std::vector<LevelInteger> levels(num_siblings);
    std::vector<CompressionValue<sizeof(T)>> values(num_siblings);

    LevelInteger current_ref_level = initial_max_refinement_level;

    /* At the beginning the coarse value indices directly correspond to the vector indices of the coarsened values */
    std::vector<size_t> coarse_value_indices(prefix_pyramid.GetSizeOfCurrentByteValues() + 1);
    std::iota(coarse_value_indices.begin(), coarse_value_indices.end(), 0);
    std::vector<size_t> coarse_value_indices_new;
    coarse_value_indices_new.reserve(coarse_value_indices.size());

    while (IsPrefixExtractionProgressing(current_ref_level))
    {
        int cv_idx_accessor = 0;
        /* We emplace a new PrefixData struct for this iteration of prefix extraction */
        /* We allocate a loose upper bound of memory */
        const size_t size_hint = prefix_pyramid.GetSizeOfCurrentByteValues();
        prefix_pyramid.EmplaceBackPrefixData(size_hint);

        /* Get the previous level and value data */
        const std::vector<LevelInteger>& previous_levels = prefix_pyramid.GetPreviousLevels();
        const std::vector<CompressionValue<sizeof(T)>>& previous_byte_values = prefix_pyramid.GetPreviousByteValues();

        int iiii = 0;
        size_t val_count = 0;
        for (auto liter = previous_levels.begin(); liter != previous_levels.end(); ++liter, ++iiii)
        {
            //cmc_debug_msg("Index: ", iiii, ", hat level: ", static_cast<int>(*liter));
            val_count += GetInitialValueCoverage(num_siblings, 4, *liter);
        }
        //cmc_debug_msg("Damit kommt man auf Num elemente: ", val_count);

        coarse_value_indices_new.push_back(0);

        /* Iterate over all values and try to coarsen all element families */
        for (size_t index = 0; index < previous_byte_values.size();)
        {
            //cmc_debug_msg("Before Xetraction, index: ", index, " Coarse Value Pos: ", prefix_pyramid.GetCurrentCoarseValueIndex());
            const bool is_family = GetNextValuesToEvaluate(current_ref_level, values, levels, index, previous_byte_values, previous_levels);
            
            if (is_family)
            {
                /* If it is a family, we are going to extract a prefix */
                auto [is_common_prefix_present, prefix] = EvaluateCommonPrefix<T>(values);

                if (is_common_prefix_present)
                {
                    //cmc_debug_msg("Prefix extraction takes place");
                    /* Store the prefix and set the indication bits */
                    prefix_pyramid.PushBackLevel(levels.front() - 1);
                    prefix_pyramid.PushBackPrefix(std::move(prefix));
                    prefix_pyramid.PushBackCoarseningBit(kIndicateCoarsening);
                    prefix_pyramid.PushBackPrefixIndicationBit(kIndicatePrefixHasBeenExtracted);

                    /* When we extract a prefix, we need to update the last byte values accordingly */
                    /* We need to trim the previous prefixes by the extracted common prefix */
                    prefix_pyramid.TrimPreviousByteValues(index, num_siblings, coarse_value_indices);

                    /* Update the index by the process family of elements */
                    index += num_siblings;

                    coarse_value_indices_new.push_back(coarse_value_indices[cv_idx_accessor + num_siblings]);
                    cv_idx_accessor += num_siblings;

                    //for (size_t elem_id = 0; elem_id < num_siblings; ++elem_id)
                    //{
                    //    const LevelInteger initial_level = prefix_pyramid.GetInitialCoarseValueLevel(prefix_pyramid.GetCurrentCoarseValueIndex());
                    //    prefix_pyramid.AddToInitialValuePos(GetInitialValueCoverage(num_siblings, initial_level, current_ref_level) + (initial_level == current_ref_level ? 0 : 1));
                    //}
                } else
                {
                    //cmc_debug_msg("NOOOOO Prefix extraction takes place");
                    /* In case we do not have a common prefix to extract */
                    prefix_pyramid.PushBackLevel(levels.front() - 1);
                    prefix_pyramid.PushBackPrefix(CompressionValue<sizeof(T)>());
                    prefix_pyramid.PushBackCoarseningBit(kIndicateCoarsening);
                    prefix_pyramid.PushBackPrefixIndicationBit(kIndicateNoPrefix);

                    /* Update the index by the process family of elements */
                    index += num_siblings;

                    for (size_t elem_id = 0; elem_id < num_siblings; ++elem_id)
                    {
                        const LevelInteger initial_level = prefix_pyramid.GetInitialCoarseValueLevel(prefix_pyramid.GetCurrentCoarseValueIndex());
                        prefix_pyramid.AddToInitialValuePos(GetInitialValueCoverage(num_siblings, initial_level, current_ref_level) + (initial_level == current_ref_level ? 0 : 0));
                    }

                    //coarse_value_indices[cv_idx_accessor + 1] = coarse_value_indices[cv_idx_accessor] + num_siblings;
                    //++cv_idx_accessor;

                    coarse_value_indices_new.push_back(coarse_value_indices[cv_idx_accessor + num_siblings]);
                    cv_idx_accessor += num_siblings;
                }
            } else
            {
                /* In case we do not have a family, the element stays unchanged and we drag the value along
                 * until it can be coarsened until the siblings are on the same refinement level */
                prefix_pyramid.PushBackLevel(levels.front());
                prefix_pyramid.PushBackPrefix(values.front());
                prefix_pyramid.PushBackCoarseningBit(kIndicateElementStaysUnchanged);
                prefix_pyramid.PushBackPrefixIndicationBit(kIndicateNoPrefix);
                
                /* Since we copied the whole value as a prefix, we need to erase it from the previous prefixes */
                if (current_ref_level < initial_max_refinement_level)
                {
                    prefix_pyramid.ResetPreviousByteValue(index);
                }
                /* Update the index by a single element */
                ++index;

                prefix_pyramid.AddToInitialValuePos(1);

                //coarse_value_indices[cv_idx_accessor + 1] = coarse_value_indices[cv_idx_accessor] + 1;
                //++cv_idx_accessor;

                coarse_value_indices_new.push_back(coarse_value_indices[cv_idx_accessor + 1]);
                ++cv_idx_accessor;
            }
            //cmc_debug_msg("After Xetraction, index: ", index, " Coarse Value Pos: ", prefix_pyramid.GetCurrentCoarseValueIndex());
        }

        std::swap(coarse_value_indices, coarse_value_indices_new);
        coarse_value_indices_new.clear();

        //if (current_ref_level == initial_max_refinement_level)
        //{
        //    const std::vector<CompressionValue<sizeof(T)>>& extra_init_byte_vals = prefix_pyramid.GetPreviousByteValues();
        //
        //    int ii = 0;
        //    for (auto viter = extra_init_byte_vals.begin(); viter != extra_init_byte_vals.end(); ++viter, ++ii)
        //    {
        //        cmc_debug_msg("Index: ", ii, ", Num Significant Bits: ", viter->GetCountOfSignificantBits());
        //    }
        //}
        //cmc_debug_msg("After level: ", current_ref_level, " Num Elements: ", prefix_pyramid.GetCurrentByteValues().size());

        /* Update the reference level such that the extraction continues on the next coarser level */
        --current_ref_level;
        prefix_pyramid.ResetInitialValuePos();

        //cmc_debug_msg("After extraction num data vals: ", prefix_pyramid.GetSizeOfCurrentByteValues());
    }

    #if 0
    const std::vector<CompressionValue<sizeof(T)>>& suffixes = prefix_pyramid.GetInitialByteValues();
    const std::vector<LevelInteger>& suff_level = prefix_pyramid.GetInitialLevels();
    int ii = 0;
    for (auto viter = suffixes.begin(); viter != suffixes.end(); ++viter, ++ii)
    {
        cmc_debug_msg("Index: ", ii, ", level: ", static_cast<int>(suff_level[ii]), " Num Significant Bits: ", viter->GetCountOfSignificantBits());
    }
    #endif
}

template <typename T>
auto TryPrediction(const std::vector<double>& errors, const std::vector<T>& coarsened_values)
-> std::enable_if_t<std::is_unsigned_v<T>, void>
{
    cmc_debug_msg("UNSIGNED passiert nichts");
}
template <typename T>
auto TryPrediction(const std::vector<double>& errors, const std::vector<T>& coarsened_values)
-> std::enable_if_t<std::is_signed_v<T>, void>
{
    cmc_assert(coarsened_values.size() > 4);
    cmc_assert(errors.size() == coarsened_values.size());

    std::vector<T> last_vals{coarsened_values[0], coarsened_values[1], coarsened_values[2], coarsened_values[3]};

    size_t index = 4;

    int num_fits = 0;

    for (auto val_iter = std::next(coarsened_values.begin(), 4); val_iter != coarsened_values.end(); ++val_iter, ++index)
    {
        const T const_fit = last_vals.back();

        const T lin_fit = 2 * last_vals[3] - last_vals[2];

        const T quad_fit = 3 * last_vals[3] - 3 * last_vals[2] + last_vals[1];

        const double const_fit_err = static_cast<double>(std::abs(coarsened_values[index] - const_fit));
        const double lin_fit_err = static_cast<double>(std::abs(coarsened_values[index] - lin_fit));
        const double quad_fit_err = static_cast<double>(std::abs(coarsened_values[index] - quad_fit));

        double err;
        double val;

        if (const_fit_err <= lin_fit_err && const_fit_err <= quad_fit_err)
        {
            err = const_fit_err;
            val = const_fit;
            //cmc_debug_msg("Constant fit with error: ", err, " and permitted error: ", errors[index]);
        } else if (lin_fit_err <= const_fit_err && lin_fit <= quad_fit_err)
        {
            err = lin_fit_err;
            val = lin_fit;
            //cmc_debug_msg("Linear fit with error: ", err, " and permitted error: ", errors[index]);
        } else if (quad_fit_err <= const_fit_err && quad_fit_err <= lin_fit_err)
        {
            err = lin_fit_err;
            val = lin_fit;
            //cmc_debug_msg("Quad fit with error: ", err, " and permitted error: ", errors[index]);
        } else
        {
            //cmc_debug_msg("No error case has been ecvaluated!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!");
        }

        if (err <= errors[index])
        {
            //cmc_debug_msg("Fit CAN be applied\n");
            last_vals[0] = last_vals[1];
            last_vals[1] = last_vals[2];
            last_vals[2] = last_vals[3];
            last_vals[3] = val;

            ++num_fits;
        } else
        {
            //cmc_debug_msg("Fit CANNOT be applied\n");
            last_vals[0] = last_vals[1];
            last_vals[1] = last_vals[2];
            last_vals[2] = last_vals[3];
            last_vals[3] = coarsened_values[index];
        }
    }

    cmc_debug_msg("Num fits with remaining error: ", num_fits);
}


template<typename T>
CompressionValue<sizeof(T)>
GetMaximumTailToggledValue(const VectorView<T>& initial_values, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value, const PermittedError& error_threshold) 
{
    bool is_toogling_progressing = true;
    CompressionValue<sizeof(T)> toggled_value = initial_serialized_value;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_toogling_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<sizeof(T)> save_previous_value = toggled_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        toggled_value.ToggleTailUntilNextUnsetBit();
        const T reinterpreted_value = toggled_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        if (not IsErrorCompliant(initial_values, reinterpreted_value, error_threshold, missing_value))
        {
            /* Revert the changes to the value */
            toggled_value = save_previous_value;
            is_toogling_progressing = false;
        }

        ++iteration_count;
    }

    return toggled_value;
}

template<typename T>
CompressionValue<sizeof(T)>
GetMaximumTailClearedValue(const VectorView<T>& initial_values, const CompressionValue<sizeof(T)>& initial_serialized_value, const T& missing_value, const PermittedError& error_threshold) 
{
    bool is_clearing_progressing = true;
    CompressionValue<sizeof(T)> cleared_value = initial_serialized_value;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<sizeof(T)> save_previous_value = cleared_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        cleared_value.ClearNextSetBitFromTail();
        const T reinterpreted_value = cleared_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        if (not IsErrorCompliant(initial_values, reinterpreted_value, error_threshold, missing_value))
        {
            /* Revert the changes to the value */
            cleared_value = save_previous_value;
            is_clearing_progressing = false;
        }

        ++iteration_count;
    }

    return cleared_value;
}


template<typename T>
void
PerformTailTruncation(const CompressionBlock<T>& compression_block, const std::vector<LevelInteger>& levels, std::vector<CompressionValue<sizeof(T)>>& coarsened_values, const T missing_value, const PermittedError& error_threshold)
{
    const int dimensionality = compression_block.GetDimensionality();
    const size_t num_siblings = GetNumSiblings(dimensionality);
    const size_t initial_reference_level = static_cast<size_t>(GetInitialRefinementLevel(dimensionality));

    const std::vector<T>& initial_values = compression_block.GetInitialData();
    
    size_t initial_val_index{0};
    size_t data_accessor{0};

    //cmc_debug_msg("\n\n Size of levels: ", levels.size(), "\n\n");
    for (auto level_iter = levels.begin(); level_iter != levels.end(); ++level_iter, ++data_accessor)
    {
        /* Get the covered domain of the coarse value */
        const size_t size_intial_value_coverage = GetInitialValueCoverage(num_siblings, initial_reference_level, static_cast<size_t>(*level_iter));
        
        /* Create a view over the domain */
        const VectorView<T> initial_values_view(initial_values.data() + initial_val_index, size_intial_value_coverage);

        /* Create truncated values */
        const CompressionValue<sizeof(T)> toggled_value = GetMaximumTailToggledValue(initial_values_view, coarsened_values[data_accessor], missing_value, error_threshold);
        const CompressionValue<sizeof(T)> cleared_value = GetMaximumTailClearedValue(initial_values_view, coarsened_values[data_accessor], missing_value, error_threshold);
        
        /* Update the index */
        initial_val_index += size_intial_value_coverage;

        /* Check which approach leads to more zero bits at the end */
        const int num_toogled_trailing_zeros = toggled_value.GetNumberTrailingZeros();
        const int num_cleared_trailing_zeros = cleared_value.GetNumberTrailingZeros();

        /* Replace the initial value with the transformed one */
        if (num_toogled_trailing_zeros >= num_cleared_trailing_zeros)
        {
            /* If the toggling approach has been more successfull */
            coarsened_values[data_accessor] = toggled_value;
        } else
        {
            /* If the clearing approach has been more successfull */
            coarsened_values[data_accessor] = cleared_value;
        }

        /* Update the trail bit count for the new value */
        coarsened_values[data_accessor].UpdateTrailBitCount();

    }
}


template<typename T>
void CompressLosslessly(const InputVariable<T>& var)
{
    std::vector<CompressionBlock<T>> compression_blocks = GetCompressionBlocks(var);
    cmc_debug_msg("Number of compression blocks: ", compression_blocks.size());

    std::vector<uint8_t> serialized_variable;
    serialized_variable.reserve(compression_blocks.size() * sizeof(T) * 200);

    /* Iterate over all gathered compression blocks */
    for (auto block_iter = compression_blocks.begin(); block_iter != compression_blocks.end(); ++block_iter)
    {
        //TryPrediction(GetErrorsForCoarsenedValues(*block_iter, levels, coarsened_values, var.GetMissingValue(), error_thresholds), coarsened_values);
        if (not block_iter->IsMortonOrdered())
        {
            block_iter->OrderMortonCurveCompliant();
        }

        const int dimensionality = block_iter->GetDimensionality();
        const size_t num_siblings = GetNumSiblings(dimensionality);
        const size_t initial_refinement_level = static_cast<size_t>(GetInitialRefinementLevel(dimensionality));
        
        const size_t num_elems = block_iter->size();

        std::vector<LevelInteger> data_levels(num_elems, initial_refinement_level);
        std::vector<T> initial_values = block_iter->GetInitialData();

        std::vector<CompressionValue<sizeof(T)>> byte_vals = TransformCompressionToByteValues(std::move(initial_values));
        //cmc_debug_msg("\n\n VOr tail truncation\n");

        PrefixPyramid<T> prefix_pyramid(block_iter->GetBlockID(), block_iter->GetHyperslab(), std::move(data_levels), std::move(byte_vals));

        prefix_pyramid.ExtractPrefixes();

        //cmc_debug_msg("Encode and serialize data is called");
        bit_vector::BitVector serialized_data = prefix_pyramid.EncodeAndSerializeData();

        size_t encoded_block_size = serialized_data.size();
        size_t encoded_block_size_bits = serialized_data.size_bits();

        //cmc_debug_msg("Size of serialized data: ", encoded_block_size, " und in bits: ", encoded_block_size_bits);

        /* Currently, just add two placeholders */
        serialized_variable.emplace_back(0);
        serialized_variable.emplace_back(0);
        
        /* Insert the serialized block */
        std::copy_n(serialized_data.data(), encoded_block_size, std::back_inserter(serialized_variable)); 

        //std::vector<CompressionValue<sizeof(T)>> byte_vals2 = prefix_pyramid.GetInitialByteValues();
        //PerformTailTruncation(*block_iter, levels_backup, byte_vals2, var.GetMissingValue(), error_thresholds);

        if constexpr (std::is_same_v<T, float>)
        {
            #if 0
            //const std::vector<T>& cvals = cvals_backup;
            const std::vector<CompressionValue<sizeof(T)>>& compr_vals = prefix_pyramid.GetInitialByteValues();

            for (size_t jjjj = 0; jjjj < cvals_backup.size(); ++jjjj)
            {
                uint32_t vval;
                std::memcpy(&vval, static_cast<void*>(&cvals_backup[jjjj]), 4);
                cmc_debug_msg(std::bitset<32>(vval), " fuer Elem: ", jjjj, ", original value: ", cvals_backup[jjjj], " und num signioficat bits: ", compr_vals[jjjj].GetCountOfSignificantBits());

            }
            #endif

            #if 0
            std::vector<CompressionValue<sizeof(T)>> cvvvals = TransformCompressionToByteValues(std::move(cvals_backup2));

            cmc_debug_msg("\n\n VOr tail truncation\n");

            PerformTailTruncation(*block_iter, levels_backup, cvvvals, var.GetMissingValue(), error_thresholds);

            cmc_debug_msg("\n\n Nach tail truncation\n\n");
            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                T reinterpreted_value = cvvvals[jjjj].template ReinterpretDataAs<T>();
                uint32_t vval;
                std::memcpy(&vval, static_cast<void*>(&reinterpreted_value), 4);
                cmc_debug_msg(std::bitset<32>(vval), " fuer Elem: ", jjjj, ", real value: ", reinterpreted_value, " und num signioficat bits: ", cvvvals[jjjj].GetCountOfSignificantBits(), " und num trailing bits: ", cvvvals[jjjj].GetNumberTrailingZeros());

            }

            #endif

            //const std::vector<CompressionValue<sizeof(T)>>& cvvvals = byte_vals2; //prefix_pyramid.GetInitialByteValues();

            #if 0
            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                T reinterpreted_value = cvvvals[jjjj].template ReinterpretDataAs<T>();
                uint32_t vval;
                std::memcpy(&vval, static_cast<void*>(&reinterpreted_value), 4);
                cmc_debug_msg(std::bitset<32>(vval), " fuer Elem: ", jjjj, ", real value: ", reinterpreted_value, " und num significant bits: ", cvvvals[jjjj].GetCountOfSignificantBits(), " und num trailing bits: ", cvvvals[jjjj].GetNumberTrailingZeros());

            }

            #endif


            #if 0
            cmc_debug_msg("\n\nPrefLengths:\n");

            std::vector<int> freq(32, 0);

            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                int sigbits = cvvvals[jjjj].GetCountOfSignificantBits();
                std::cout << sigbits << ", ";
                freq[sigbits] += 1;
            }
            cmc_debug_msg("end");

            cmc_debug_msg("Frequency Prefix Lengths: ");
            int index = 0;
            for (auto fiter = freq.begin(); fiter != freq.end(); ++fiter, ++index)
            {
                if (*fiter > 0)
                {
                    cmc_debug_msg("Length: ", index, ", Freq: ", *fiter);
                }
            }
            #endif

            #if 0
            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                int sigbits = cvvvals[jjjj].GetCountOfSignificantBits();
                cmc_debug_msg(std::bitset<8>(sigbits));
            }
            #endif

        }


    }

    cmc_debug_msg("The whole variable has been serialized into: ", serialized_variable.size(), " bytes.");
}


template<typename T>
void Compress(const InputVariable<T>& var)
{
    //cmc_assert(var.IsValid());

    std::vector<CompressionBlock<T>> compression_blocks = GetCompressionBlocks(var);
    cmc_debug_msg("Number of compression blocks: ", compression_blocks.size());

    std::vector<uint8_t> serialized_variable;
    serialized_variable.reserve(compression_blocks.size() * sizeof(T) * 200);

    /* Iterate over all gathered compression blocks */
    for (auto block_iter = compression_blocks.begin(); block_iter != compression_blocks.end(); ++block_iter)
    {
        PermittedError error_thresholds = GetPermittedError();

        /* Coarsen the values as much as possible */
        auto [levels, coarsened_values] = CoarsenValuesAdaptively<T>(var, *block_iter, error_thresholds);
        
        //std::vector<T> cvals_backup = coarsened_values;
        //std::vector<T> cvals_backup2 = coarsened_values;
        //std::vector<LevelInteger> levels_backup = levels;

        //cmc_debug_msg("\n\n\nAfter coarsen values:");
        //for (auto iter = levels.begin(); iter != levels.end(); ++iter)
        //{
        //    cmc_debug_msg("levels: ", static_cast<int>(*iter));
        //}
        //cmc_debug_msg("\n\n\nBefore get error");
        /* Obtain the errors that hhave been introduced by the adaptive coarsening process */
        std::vector<double> errors = GetErrorsAfterCoarsening(*block_iter, levels, coarsened_values, var.GetMissingValue(), error_thresholds);

        //TryPrediction(GetErrorsForCoarsenedValues(*block_iter, levels, coarsened_values, var.GetMissingValue(), error_thresholds), coarsened_values);

        //std::vector<CompressionValue<sizeof(T)>> cvvvals = TransformCompressionToByteValues(std::move(cvals_backup2));

        std::vector<CompressionValue<sizeof(T)>> byte_vals = TransformCompressionToByteValues(std::move(coarsened_values));
        //cmc_debug_msg("\n\n VOr tail truncation\n");

        PerformTailTruncation(*block_iter, levels, byte_vals, var.GetMissingValue(), error_thresholds);
        //cmc_debug_msg("Perform Tail Truncation abgeschlossen");

        
        //for (size_t jjjj = 0; jjjj < byte_vals.size(); ++jjjj)
        //{
        //    T reinterpreted_value = byte_vals[jjjj].template ReinterpretDataAs<T>();
        //    uint32_t vval;
        //    std::memcpy(&vval, static_cast<void*>(&reinterpreted_value), 4);
        //    cmc_debug_msg(std::bitset<32>(vval), " fuer Elem: ", jjjj, ", real value: ", reinterpreted_value, " und num significant bits: ", byte_vals[jjjj].GetCountOfSignificantBits(), " und num trailing bits: ", byte_vals[jjjj].GetNumberTrailingZeros());
        //}

        //std::exit(1);
        PrefixPyramid<T> prefix_pyramid(block_iter->GetBlockID(), block_iter->GetHyperslab(), std::move(levels), std::move(byte_vals));

        //cmc_debug_msg("Before prefix pyramid is built");
        //cmc_debug_msg("Before extract prefixes");
        //XorValues(coarsened_values);
       
       //PrefixPyramid<T> prefix_pyramid(block_iter->GetBlockID(), block_iter->GetHyperslab(), std::move(levels), std::move(coarsened_values));

       
        prefix_pyramid.ExtractPrefixes();

        //cmc_debug_msg("Encode and serialize data is called");
        bit_vector::BitVector serialized_data = prefix_pyramid.EncodeAndSerializeData();

        size_t encoded_block_size = serialized_data.size();
        size_t encoded_block_size_bits = serialized_data.size_bits();

        //cmc_debug_msg("Size of serialized data: ", encoded_block_size, " und in bits: ", encoded_block_size_bits);

        /* Currently, just add two placeholders */
        serialized_variable.emplace_back(0);
        serialized_variable.emplace_back(0);
        
        /* Insert the serialized block */
        std::copy_n(serialized_data.data(), encoded_block_size, std::back_inserter(serialized_variable)); 

        //std::vector<CompressionValue<sizeof(T)>> byte_vals2 = prefix_pyramid.GetInitialByteValues();
        //PerformTailTruncation(*block_iter, levels_backup, byte_vals2, var.GetMissingValue(), error_thresholds);

        if constexpr (std::is_same_v<T, float>)
        {
            #if 0
            //const std::vector<T>& cvals = cvals_backup;
            const std::vector<CompressionValue<sizeof(T)>>& compr_vals = prefix_pyramid.GetInitialByteValues();

            for (size_t jjjj = 0; jjjj < cvals_backup.size(); ++jjjj)
            {
                uint32_t vval;
                std::memcpy(&vval, static_cast<void*>(&cvals_backup[jjjj]), 4);
                cmc_debug_msg(std::bitset<32>(vval), " fuer Elem: ", jjjj, ", original value: ", cvals_backup[jjjj], " und num signioficat bits: ", compr_vals[jjjj].GetCountOfSignificantBits());

            }
            #endif

            #if 0
            std::vector<CompressionValue<sizeof(T)>> cvvvals = TransformCompressionToByteValues(std::move(cvals_backup2));

            cmc_debug_msg("\n\n VOr tail truncation\n");

            PerformTailTruncation(*block_iter, levels_backup, cvvvals, var.GetMissingValue(), error_thresholds);

            cmc_debug_msg("\n\n Nach tail truncation\n\n");
            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                T reinterpreted_value = cvvvals[jjjj].template ReinterpretDataAs<T>();
                uint32_t vval;
                std::memcpy(&vval, static_cast<void*>(&reinterpreted_value), 4);
                cmc_debug_msg(std::bitset<32>(vval), " fuer Elem: ", jjjj, ", real value: ", reinterpreted_value, " und num signioficat bits: ", cvvvals[jjjj].GetCountOfSignificantBits(), " und num trailing bits: ", cvvvals[jjjj].GetNumberTrailingZeros());

            }

            #endif

            //const std::vector<CompressionValue<sizeof(T)>>& cvvvals = byte_vals2; //prefix_pyramid.GetInitialByteValues();

            #if 0
            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                T reinterpreted_value = cvvvals[jjjj].template ReinterpretDataAs<T>();
                uint32_t vval;
                std::memcpy(&vval, static_cast<void*>(&reinterpreted_value), 4);
                cmc_debug_msg(std::bitset<32>(vval), " fuer Elem: ", jjjj, ", real value: ", reinterpreted_value, " und num significant bits: ", cvvvals[jjjj].GetCountOfSignificantBits(), " und num trailing bits: ", cvvvals[jjjj].GetNumberTrailingZeros());

            }

            #endif


            #if 0
            cmc_debug_msg("\n\nPrefLengths:\n");

            std::vector<int> freq(32, 0);

            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                int sigbits = cvvvals[jjjj].GetCountOfSignificantBits();
                std::cout << sigbits << ", ";
                freq[sigbits] += 1;
            }
            cmc_debug_msg("end");

            cmc_debug_msg("Frequency Prefix Lengths: ");
            int index = 0;
            for (auto fiter = freq.begin(); fiter != freq.end(); ++fiter, ++index)
            {
                if (*fiter > 0)
                {
                    cmc_debug_msg("Length: ", index, ", Freq: ", *fiter);
                }
            }
            #endif

            #if 0
            for (size_t jjjj = 0; jjjj < cvvvals.size(); ++jjjj)
            {
                int sigbits = cvvvals[jjjj].GetCountOfSignificantBits();
                cmc_debug_msg(std::bitset<8>(sigbits));
            }
            #endif

        }


    }

    cmc_debug_msg("The whole variable has been serialized into: ", serialized_variable.size(), " bytes.");
}

#if 0

01000011011110011 011110000111111 fuer Elem: 42, original value: 249.735 und num signioficat bits: 17
01000011011110011 010111100100010 fuer Elem: 43, original value: 249.684 und num signioficat bits: 17
01000011011110011 111010100111010 fuer Elem: 44, original value: 249.958 und num signioficat bits: 17
01000011011110011 011101010101111 fuer Elem: 45, original value: 249.729 und num signioficat bits: 17

als bit plane 0010 1111 1011 1101 1110 0101 0110 0001 0000 1111 1010 1011 1001 1111 1001
normal memory 0111 1000 0111 1110 1011 1100 1000 1011 1010 1001 1101 0011 1010 1010 1111

              {0,1} x 3Bit für Laenge 5,6,7,8,9,10,11,12

Überlegung: Enkodiert man die Länge des Suffix oder die Anzahl Bits die gekürzt werden? Für kleinere erlaubte Fehler bietet sich letzteres vielleicht eher an
Längster Datentyp der betrachet wird: 64 bit -> diese Länge könnte man in 6 Bits kodieren, aber vermutlich reichen meistens 32 Bit aus, wenn nicht sogar 16

-> Vielleicht wäre eine Kürzung um maximal 16 Bit realistisch -> Kodierung in 4 Bit pro Wert, oder sogar nur 3 Bit Kodierung -> 0 bis 6 Bits Kürzung und 7 als flag für längere Kürzung
#endif


struct PerformPatchCompression
{
public:
    void operator()(const InputVariable<int8_t>& var) {
        Compress<int8_t>(var);
    }
    void operator()(const InputVariable<char>& var) {
        Compress<char>(var);
    }
    void operator()(const InputVariable<int16_t>& var) {
        Compress<int16_t>(var);
    }
    void operator()(const InputVariable<int32_t>& var) {
        Compress<int32_t>(var);
    }
    void operator()(const InputVariable<float>& var) {
        Compress<float>(var);
    }
    void operator()(const InputVariable<double>& var) {
        Compress<double>(var);
    }
    void operator()(const InputVariable<uint8_t>& var) {
        Compress<uint8_t>(var);
    }
    void operator()(const InputVariable<uint16_t>& var) {
        Compress<uint16_t>(var);
    }
    void operator()(const InputVariable<uint32_t>& var) {
        Compress<uint32_t>(var);
    }
    void operator()(const InputVariable<int64_t>& var) {
        Compress<int64_t>(var);
    }
    void operator()(const InputVariable<uint64_t>& var) {
        Compress<uint64_t>(var);
    }
};


void
Compressor::Compress()
{
    /* Iterate over all variables */
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        /* Apply the scaling and offset of the data */
        var_iter->ApplyScalingAndOffset();

        /* Peform the blocked compression for this variable */
        std::visit(PerformPatchCompression(), var_iter->GetInternalVariant());
    }
}

}

}
