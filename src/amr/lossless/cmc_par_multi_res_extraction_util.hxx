#ifndef CMC_PAR_MULTI_RES_EXTRACTION_UTIL_HXX
#define CMC_PAR_MULTI_RES_EXTRACTION_UTIL_HXX

#include "utilities/cmc_bits.hxx"
#include "utilities/cmc_bits_vector.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"

#include <array>
#include <limits>
#include <execution>

#include <bitset>

namespace cmc::par::lossless::multi_res
{

template<int32_t DIM>
concept Dimension = (DIM >= 1 && DIM <= 4);

/* Set the number of adjacent data points that will be coarsened given a certain dimensionality */
template<int32_t DIM>
constexpr int kPackSize;
template<>
constexpr int kPackSize<1> = 2;
template<>
constexpr int kPackSize<2> = 4;
template<>
constexpr int kPackSize<3> = 8;
template<>
constexpr int kPackSize<4> = 16;

//TODO: Make t8code dependent
/* Set the maximum number of children elements per dimension */
template<int32_t DIM>
constexpr int kNumMaxChildrenElements;
template<>
constexpr int kNumMaxChildrenElements<1> = 2;
template<>
constexpr int kNumMaxChildrenElements<2> = 4;
template<>
constexpr int kNumMaxChildrenElements<3> = 8;

template<typename T>
concept ArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T>);

template<typename T>
concept OneByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 1));

template<typename T>
concept TwoByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 2));

template<typename T>
concept FourByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 4));

template<typename T>
concept EightByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 8));

template<typename T>
concept UnsignedIntegerType = (std::is_unsigned_v<T> && std::is_integral_v<T>);

using OneByteResidualType = uint8_t;

using TwoByteResidualType = uint16_t;

using FourByteResidualType = uint32_t;

using EightByteResidualType = uint64_t;

template<OneByteArithmeticType T>
constexpr inline
uint8_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint8_t>(value);
}

template<TwoByteArithmeticType T>
constexpr inline
uint16_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint16_t>(value);
}

template<FourByteArithmeticType T>
constexpr inline
uint32_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint32_t>(value);
}

template<EightByteArithmeticType T>
constexpr inline
uint64_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint64_t>(value);
}

using SymbolType = uint8_t;

inline constexpr int32_t kMaxPresentElementLevelUnknown = -1;

constexpr inline SymbolType kResidualSignumIndication = 0x80;

template<ArithmeticType T>
constexpr inline SymbolType kProcessEndSymbol = kResidualSignumIndication + sizeof(T) * cmc::bits::kCharBit;

template<UnsignedIntegerType T>
constexpr inline SymbolType
CreateEntropySymbol(const bool is_approx_greater, const T residual)
{
    /* Check if the residual is zero, in this case we only store a single entropy code for +/- 0 */
    if (cmc::bits::GetLZC(residual) == sizeof(T) * cmc::bits::kCharBit) [[unlikely]]
    {
        return static_cast<SymbolType>(cmc::bits::GetLZC(residual));
    }

    return ((SymbolType{is_approx_greater} << 7) | static_cast<SymbolType>(cmc::bits::GetLZC(residual)));
}

constexpr inline bool
IsApproximationGreater(const SymbolType entropy_symbol)
{
    return (entropy_symbol >> 7);
}

template<ArithmeticType T>
constexpr inline int
MapEntropySymbolToArrayIndex(const SymbolType entropy_symbol)
{
    return (entropy_symbol & SymbolType{0x7F}) + (entropy_symbol >> 7) * (sizeof(T) * cmc::bits::kCharBit + 1);
}

template<ArithmeticType T>
constexpr inline SymbolType
MapArrayIndexToEntropySymbol(int array_idx)
{
    return static_cast<SymbolType>(array_idx + (array_idx > sizeof(T) * cmc::bits::kCharBit ? kResidualSignumIndication - sizeof(T) * cmc::bits::kCharBit - 1 : 0));
}

template<ArithmeticType T>
constexpr inline int
GetNumEntropySymbols()
{
    return 2 * (sizeof(T) * cmc::bits::kCharBit + 1);
}

constexpr inline SymbolType
GetLZCFromEntropySymbol(const SymbolType entropy_symbol)
{
    return (entropy_symbol & SymbolType{0x7F});
}

template<ArithmeticType T>
constexpr inline void
AddProcessEndSymbol(std::array<uint64_t, GetNumEntropySymbols<T>()>& entropy_symbols_frequency, const uint64_t num_local_proc_end_symbols)
{
    /* We store the process-end-symbol in the last array entry */
    entropy_symbols_frequency[sizeof(T) * cmc::bits::kCharBit + kResidualSignumIndication] = num_local_proc_end_symbols;
}

template<typename T, int32_t DIM>
struct ElemEncodingData;

template<typename T, int32_t DIM>
requires Dimension<DIM> && OneByteArithmeticType<T>
struct ElemEncodingData<T, DIM>
{
    uint8_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
    std::array<OneByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && TwoByteArithmeticType<T>
struct ElemEncodingData<T, DIM>
{
    uint16_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
    std::array<TwoByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && FourByteArithmeticType<T>
struct ElemEncodingData<T, DIM>
{
    uint32_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
    std::array<FourByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && EightByteArithmeticType<T>
struct ElemEncodingData<T, DIM>
{
    uint64_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
    std::array<EightByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<ArithmeticType T, int32_t N>
constexpr inline T
ComputeArithmeticMean(const std::array<T, N>& values)
{
    T sum = static_cast<T>(0);
    for (int idx{0}; idx < N; ++idx)
    {
        sum += values[idx];
    }
    return sum / static_cast<T>(N);
}

template<ArithmeticType T, int32_t N>
constexpr inline T
ComputeMidRange(const std::array<T, N>& values)
{
    static_assert(N >= 2);

    T min = std::numeric_limits<T>::max();
    /* Find the minimum */
    for (int idx{0}; idx < N; ++idx)
    {
        if (min > values[idx])
        {
            min = values[idx];
        }
    }

    T max = std::numeric_limits<T>::lowest();
    /* Find the maximum */
    for (int idx{0}; idx < N; ++idx)
    {
        if (max < values[idx])
        {
            max = values[idx];
        }
    }

    /* Compute the mid-range */
    return ((max / 2) + (min / 2));
}

template<ArithmeticType T>
inline T
ComputeArithmeticMean(const std::span<T> values)
{
    cmc_assert(values.size() >= 2);

    T sum = static_cast<T>(0);
    for (int idx{0}; idx < values.size(); ++idx)
    {
        sum += values[idx];
    }
    return sum / static_cast<T>(values.size());
}

template<ArithmeticType T>
inline T
ComputeMidRange(const std::span<T> values)
{
    cmc_assert(values.size() >= 2);

    T min{std::numeric_limits<T>::max()};
    /* Find the minimum */
    for (int idx{0}; idx < values.size(); ++idx)
    {
        if (min > values[idx])
        {
            min = values[idx];
        }
    }

    T max{std::numeric_limits<T>::lowest()};
    /* Find the maximum */
    for (int idx{0}; idx < values.size(); ++idx)
    {
        if (max < values[idx])
        {
            max = values[idx];
        }
    }

    /* Compute the mid-range */
    return ((max / 2) + (min / 2));
}

/**
 * Compute the DIM-th root of N which resembles the amount of data points per single direction.
 * In case, there is no integral DIM-th root of N, the function returns zero, since we cannot work 
 * with such structures currently.
 */
template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int
ComputeNumDataPerDim()
{
    static_assert(N >= 1);

    for (uint64_t i{1}; i <= N; ++i)
    {
        uint64_t exp{1};
        for (int dim{0}; dim < DIM; ++dim)
        {
            exp *= i;
        }

        if (exp == N)
        {
            return i;
        }
    }

    return 0;
}

#if 0

/** 
 * Compute the number of compression levels for the amount of data points N given in a single direction.
 * This gives the amount of levels on which residuals and entropy codes exist, e.g. 3x3x3 data points in 3D
 * gives two levels; 3x3x3 -> 2x2x2 (-> afterwards the coarse value has been reached 1x1x1) */
template<int32_t N>
constexpr int
ComputeNumIntraCompressionLevels()
{
    if constexpr (N == 1)
    {
        return 0;
    }

    int exp = 1;
    for (int i{1}; i < N; ++i)
    {
        exp *= 2;
        if (exp >= N)
        {
            return i;
        }
    }

    return 0;
}

/** Compute the number of compression levels for the amount of overall data points N given on the element with the corersponding diemnsionality.
 * This gives the amount of levels on which residuals and entropy codes exist, e.g. 8x8x8 data points in 3D
 * gives three levels; 8x8x8 -> 4x4x4 -> 2x2x2 (-> afterwards the coarse value has been reached 1x1x1) */
template<int32_t DIM, int32_t N>
constexpr int
ComputeNumIntraCompressionLevels()
{
    /* Compute the data per dimension */
    constexpr int32_t num_data_per_dim = ComputeNumDataPerDim<DIM, N>();

    /* Compute the number of compression levels */
    return ComputeNumIntraCompressionLevels<num_data_per_dim>();
}

template<int32_t DIM, int32_t N>
constexpr int
ComputeNumIntraPredictors()
{
    if constexpr (N == 1)
    {
        return 0;
    }

    /* Compute the data per dimension */
    constexpr int32_t num_data_per_dim = ComputeNumDataPerDim<DIM, N>();
    
    int num_intra_predictors{0};
    int current_data_per_dim = num_data_per_dim;

    if constexpr (num_data_per_dim > 0)
    {
        /* In case we are able to work with the certain amount of data */
        constexpr int num_compression_lvls_ = ComputeNumIntraCompressionLevels<DIM, N>();

        for (int clvl{0}; clvl < num_compression_lvls_; ++clvl)
        {
            /* Compute the coarse number of elements per dimension on this level */
            const int num_lvl_data_per_dim = current_data_per_dim / 2 + (current_data_per_dim % 2 != 0 ? 1 : 0);

            int num_level_entropy_codes{1};
            /* Determine the number of entropy codes on this level per element */
            for (int i{0}; i < DIM; ++i)
            {
                num_level_entropy_codes *= num_lvl_data_per_dim;
            }

            /* Update the overall count by this level */
            num_intra_predictors += num_level_entropy_codes;

            /* Update the number of data per dimension */
            current_data_per_dim = num_lvl_data_per_dim;
        }

        /* And add the last coarse elem value as predictor as well */
        ++num_intra_predictors;
    } else
    {
        static_assert(false, "Currently, the amount of data points N needs to be an integral DIM-th root; e.g. in 3D: N=1,8,27,64,...");
    }

    return num_intra_predictors;
}

/**
 * Compute the overall number of entropy codes/residuals per elemnent based on the dimensionality DIM
 * and the number of data points N per element.
 */
//TODO: Incorrect computation since we only have one incomplete pack per level since we iterate linearily through the data
//and make use of an SFC ordering at this position
template<int32_t DIM, int32_t N>
constexpr inline int32_t
ComputeNumPyramidalCodes()
{
    /* Compute the data per dimension */
    constexpr int32_t num_data_per_dim = ComputeDataPerDim<DIM, N>();
    
    int num_entropy_codes_per_elem{N};
    int current_data_per_dim = num_data_per_dim;

    if constexpr (num_data_per_dim > 0)
    {
        /* In case of we are able to work with the certain amount of data */
        constexpr int num_compression_lvls_ = ComputeNumIntraCompressionLevels<DIM, N>();

        for (int clvl{0}; clvl < num_compression_lvls_; ++clvl)
        {
            /* Compute the coarse number of elements per dimension on this level */
            const int num_lvl_data_per_dim = current_data_per_dim / 2 + (current_data_per_dim % 2 != 0 ? 1 : 0);

            int num_level_entropy_codes{1};
            /* Determine the number of entropy codes on this level per element */
            for (int i{0}; i < DIM; ++i)
            {
                num_level_entropy_codes *= num_lvl_data_per_dim;
            }

            /* Update the overall count by this level */
            num_entropy_codes_per_elem += num_level_entropy_codes;

            /* Update the number of data per dimension */
            current_data_per_dim = num_lvl_data_per_dim;
        }
    } else
    {
        static_assert(false, "Currently, the amount of data points N needs to be an integral DIM-th root; e.g. in 3D: N=1,8,27,64,...");
    }

    return num_entropy_codes_per_elem;
}

#else

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int
ComputeNumIntraCompressionLevels()
{
    static_assert(N > 0, "There needs to be data given (N > 0) on the element.");

    if constexpr (N == 1)
    {
        return 0;
    }

    int32_t exp = 1;
    for (int i{1}; i <= N; ++i)
    {
        exp *= kPackSize<DIM>;
        if (exp >= N)
        {
            return i;
        }
    }

    return 0;
}

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int
ComputeNumIntraPredictors()
{
    static_assert(N > 0, "There needs to be data given (N > 0) on the element.");

    if constexpr (N == 1)
    {
        return 0;
    }

    int num_intra_predictors{0};
    int32_t current_lvl_num_data{N};

    /* In case we are able to work with the certain amount of data */
    constexpr int num_compression_lvls_ = ComputeNumIntraCompressionLevels<DIM, N>();

    for (int clvl{0}; clvl < num_compression_lvls_; ++clvl)
    {
        /* Compute the coarse number of elements on this level */
        const int num_lvl_data = current_lvl_num_data / kPackSize<DIM> + (current_lvl_num_data % kPackSize<DIM> != 0 ? 1 : 0);

        /* Update the overall count by this level */
        num_intra_predictors += num_lvl_data;

        /* Update the number of data */
        current_lvl_num_data = num_lvl_data;
    }

    return num_intra_predictors;
}

/**
 * Compute the overall number of entropy codes/residuals per element based on the dimensionality DIM
 * and the number of data points N per element.
 */
template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr inline int32_t
ComputeNumPyramidalCodes()
{
    static_assert(N > 0, "There needs to be data given (N > 0) on the element.");
    
    if constexpr (N == 1)
    {
        return 0;
    }

    int num_entropy_codes_per_elem{0};
    int32_t current_lvl_num_data{N};

    /* In case of we are able to work with the certain amount of data */
    constexpr int num_compression_lvls_ = ComputeNumIntraCompressionLevels<DIM, N>();

    for (int clvl{0}; clvl < num_compression_lvls_; ++clvl)
    {
        /* Update the number of entropy codes on this level */
        num_entropy_codes_per_elem += current_lvl_num_data;

        /* Update the number of data for the next level */
        current_lvl_num_data = current_lvl_num_data / kPackSize<DIM> + (current_lvl_num_data % kPackSize<DIM> != 0 ? 1 : 0);
    }

    return num_entropy_codes_per_elem;
}

#endif

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int kNumIntraCompressionLevels = ComputeNumIntraCompressionLevels<DIM, N>();

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int kNumPyramidalCodes = ComputeNumPyramidalCodes<DIM, N>();

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int kNumIntraPredictors = ComputeNumIntraPredictors<DIM, N>();

template<typename T, int32_t DIM, int32_t N>
struct IntraElementCoding;

template<typename T, int32_t DIM, int32_t N>
requires Dimension<DIM> && OneByteArithmeticType<T>
struct IntraElementCoding<T, DIM, N>
{
    T GetCoarseValuePredictor() const {return coarse_value_predictor;}
    
    T coarse_value_predictor{};
    std::array<OneByteResidualType, kNumPyramidalCodes<DIM, N>> residuals{};
    std::array<SymbolType, kNumPyramidalCodes<DIM, N>> entropy_codes{};
};

template<typename T, int32_t DIM, int32_t N>
requires Dimension<DIM> && TwoByteArithmeticType<T>
struct IntraElementCoding<T, DIM, N>
{
    T GetCoarseValuePredictor() const {return coarse_value_predictor;}
    
    T coarse_value_predictor{};
    std::array<TwoByteResidualType, kNumPyramidalCodes<DIM, N>> residuals{};
    std::array<SymbolType, kNumPyramidalCodes<DIM, N>> entropy_codes{};
};

template<typename T, int32_t DIM, int32_t N>
requires Dimension<DIM> && FourByteArithmeticType<T>
struct IntraElementCoding<T, DIM, N>
{
    T GetCoarseValuePredictor() const {return coarse_value_predictor;}
    
    T coarse_value_predictor{};
    std::array<FourByteResidualType, kNumPyramidalCodes<DIM, N>> residuals{};
    std::array<SymbolType, kNumPyramidalCodes<DIM, N>> entropy_codes{};
};

template<typename T, int32_t DIM, int32_t N>
requires Dimension<DIM> && EightByteArithmeticType<T>
struct IntraElementCoding<T, DIM, N>
{
    T GetCoarseValuePredictor() const {return coarse_value_predictor;}
    
    T coarse_value_predictor{};
    std::array<EightByteResidualType, kNumPyramidalCodes<DIM, N>> residuals{};
    std::array<SymbolType, kNumPyramidalCodes<DIM, N>> entropy_codes{};
};


template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr inline std::array<int32_t, kNumIntraCompressionLevels<DIM, N>>
ComputeNumCodesPerLevel()
{
    static_assert(N > 0, "There needs to be data given (N > 0) on the element.");

    /* Allocate the output array */
    std::array<int32_t, kNumIntraCompressionLevels<DIM, N>> elems_per_level{};

    int32_t current_lvl_num_data{N};

    /* Iterate over all intra element compression levels */
    for (int clvl{0}; clvl < kNumIntraCompressionLevels<DIM, N>; ++clvl)
    {
        /* Update the number of entropy codes on this level */
        elems_per_level[clvl] = current_lvl_num_data;

        /* Update the number of data for the next level */
        current_lvl_num_data = current_lvl_num_data / kPackSize<DIM> + (current_lvl_num_data % kPackSize<DIM> != 0 ? 1 : 0);
    }

    return elems_per_level;
}

/**
 *  A constexpr for-loop construct
 */
template <auto START, auto END, auto INCREMENT, class FUNC>
constexpr void
for_constexpr(FUNC&& f)
{
    if constexpr (START < END)
    {
        f(std::integral_constant<decltype(START), START>());
        for_constexpr<START + INCREMENT, END, INCREMENT>(f);
    }
}

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr inline int
ComputeNumFullPacks()
{
	return N / kPackSize<DIM>;
}

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr inline int
ComputeNumFullPacks(const int lvl_iteration)
{
    int32_t num_current_data{N};
    for (int i{0}; i < lvl_iteration; ++i)
    {
        num_current_data = num_current_data / kPackSize<DIM> + (num_current_data % kPackSize<DIM> != 0 ? 1 : 0);
    }
	return num_current_data / kPackSize<DIM>;
}

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr inline int
ComputeIntraCompressionLevelIncompletePackSize(const int lvl_iteration)
{
    int num_current_data = N;
    for (int i{0}; i < lvl_iteration; ++i)
    {
        num_current_data = num_current_data / kPackSize<DIM> + (num_current_data % kPackSize<DIM> != 0 ? 1 : 0);
    }
	return num_current_data % kPackSize<DIM>;
}

template<ArithmeticType T, int32_t N>
inline std::array<T, N+2>
CreatePredictors(const std::array<T, N>& init_data)
{
    /* Utilize each individual value as a predictor */
    std::array<T, N + 2> predictors;
    std::copy_n(init_data.begin(), N, predictors.begin());

    /* Append the arithmetic mean as a predictor */
    predictors[N] = ComputeArithmeticMean<T, N>(init_data);

    /* Append the mid-range as a predictor */
    predictors[N + 1] = ComputeMidRange<T, N>(init_data);

    return predictors;
}

template<OneByteArithmeticType T, int32_t N>
inline std::tuple<T, std::array<OneByteResidualType, N>, std::array<SymbolType, N>>
PerformDefaultIntraElemLevelCompression(const std::array<T, N> init_data)
{
    static_assert(N >= 1);

    if constexpr (N == 1)
    {
        /* In case a 2:1 coarsening is used and the overall number is uneven, we just return the same value as extracted */
        return std::make_tuple(init_data[0], std::array<OneByteResidualType, 1>{0}, std::array<SymbolType, 1>{CreateEntropySymbol(false, OneByteResidualType{0})});
    } else
    {
        /* Create all predictors */
        const std::array<T, N+2> predictors = CreatePredictors<T, N>(init_data);

        int current_lzc{-1};
        T lzc_maximizing_predictor{};
        std::array<OneByteResidualType, N> residuals{};
        std::array<SymbolType, N> entropy_codes{};

        for (int pred_idx{0}; pred_idx < N + 2; ++pred_idx)
        {
            std::array<OneByteResidualType, N> current_residuals{};
            std::array<SymbolType, N> current_entropy_codes{};

            /* This computation has a symmetrical part (in those cases an init value is used as a predictor), but we are computing it fully currently */
            /* Check the predictor for all element values */
            for (int elem_idx{0}; elem_idx < N; ++elem_idx)
            {
                /* Compute the residual */
                const auto [is_approx_greater, residual] = cmc::bits::ComputeIntegerResidual<T>(predictors[pred_idx], predictors[elem_idx]);
                
                /* Store the residual */
                current_residuals[elem_idx] = residual;

                /* Store the entropy code */
                current_entropy_codes[elem_idx] = CreateEntropySymbol(is_approx_greater, residual);
            }

            /* Compute the cumulative leading zero count */
            const int cumulative_lzc = std::transform_reduce(std::execution::par_unseq, current_residuals.cbegin(), current_residuals.cend(), static_cast<int>(0),
                                                            std::plus<>{}, [](auto res){return cmc::bits::GetLZC(res);});

            /* If the predictor maximizes the LZC, we store it */
            if (current_lzc < cumulative_lzc)
            {
                current_lzc = cumulative_lzc;
                lzc_maximizing_predictor = predictors[pred_idx];
                std::copy_n(current_residuals.cbegin(), N, residuals.begin());
                std::copy_n(current_entropy_codes.cbegin(), N, entropy_codes.begin());
            }
        }

        return std::make_tuple(lzc_maximizing_predictor, residuals, entropy_codes);
    }
}

template<TwoByteArithmeticType T, int32_t N>
inline std::tuple<T, std::array<TwoByteResidualType, N>, std::array<SymbolType, N>>
PerformDefaultIntraElemLevelCompression(const std::array<T, N> init_data)
{
    static_assert(N >= 1);

    if constexpr (N == 1)
    {
        /* In case a 2:1 coarsening is used and the overall number is uneven, we just return the same value as extracted */
        return std::make_tuple(init_data[0], std::array<TwoByteResidualType, 1>{0}, std::array<SymbolType, 1>{CreateEntropySymbol(false, TwoByteResidualType{0})});
    } else
    {
        /* Create all predictors */
        const std::array<T, N+2> predictors = CreatePredictors<T, N>(init_data);

        int current_lzc{-1};
        T lzc_maximizing_predictor{};
        std::array<TwoByteResidualType, N> residuals{};
        std::array<SymbolType, N> entropy_codes{};

        for (int pred_idx{0}; pred_idx < N + 2; ++pred_idx)
        {
            std::array<TwoByteResidualType, N> current_residuals{};
            std::array<SymbolType, N> current_entropy_codes{};

            /* This computation has a symmetrical part (in those cases an init value is used as a predictor), but we are computing it fully currently */
            /* Check the predictor for all element values */
            for (int elem_idx{0}; elem_idx < N; ++elem_idx)
            {
                /* Compute the residual */
                const auto [is_approx_greater, residual] = cmc::bits::ComputeIntegerResidual<T>(predictors[pred_idx], predictors[elem_idx]);
                
                /* Store the residual */
                current_residuals[elem_idx] = residual;

                /* Store the entropy code */
                current_entropy_codes[elem_idx] = CreateEntropySymbol(is_approx_greater, residual);
            }

            /* Compute the cumulative leading zero count */
            const int cumulative_lzc = std::transform_reduce(std::execution::par_unseq, current_residuals.cbegin(), current_residuals.cend(), static_cast<int>(0),
                                                            std::plus<>{}, [](auto res){return cmc::bits::GetLZC(res);});

            /* If the predictor maximizes the LZC, we store it */
            if (current_lzc < cumulative_lzc)
            {
                current_lzc = cumulative_lzc;
                lzc_maximizing_predictor = predictors[pred_idx];
                std::copy_n(current_residuals.cbegin(), N, residuals.begin());
                std::copy_n(current_entropy_codes.cbegin(), N, entropy_codes.begin());
            }
        }

        return std::make_tuple(lzc_maximizing_predictor, residuals, entropy_codes);
    }
}

template<FourByteArithmeticType T, int32_t N>
inline std::tuple<T, std::array<FourByteResidualType, N>, std::array<SymbolType, N>>
PerformDefaultIntraElemLevelCompression(const std::array<T, N> init_data)
{
    static_assert(N >= 1);

    if constexpr (N == 1)
    {
        /* In case a 2:1 coarsening is used and the overall number is uneven, we just return the same value as extracted */
        return std::make_tuple(init_data[0], std::array<FourByteResidualType, 1>{0}, std::array<SymbolType, 1>{CreateEntropySymbol(false, FourByteResidualType{0})});
    } else
    {

        /* Create all predictors */
        const std::array<T, N+2> predictors = CreatePredictors<T, N>(init_data);

        int current_lzc{-1};
        T lzc_maximizing_predictor{};
        std::array<FourByteResidualType, N> residuals{};
        std::array<SymbolType, N> entropy_codes{};

        for (int pred_idx{0}; pred_idx < N + 2; ++pred_idx)
        {
            std::array<FourByteResidualType, N> current_residuals{};
            std::array<SymbolType, N> current_entropy_codes{};

            /* This computation has a symmetrical part (in those cases an init value is used as a predictor), but we are computing it fully currently */
            /* Check the predictor for all element values */
            for (int elem_idx{0}; elem_idx < N; ++elem_idx)
            {
                /* Compute the residual */
                const auto [is_approx_greater, residual] = cmc::bits::ComputeIntegerResidual<T>(predictors[pred_idx], predictors[elem_idx]);
                
                /* Store the residual */
                current_residuals[elem_idx] = residual;

                /* Store the entropy code */
                current_entropy_codes[elem_idx] = CreateEntropySymbol(is_approx_greater, residual);
            }

            /* Compute the cumulative leading zero count */
            const int cumulative_lzc = std::transform_reduce(std::execution::par_unseq, current_residuals.cbegin(), current_residuals.cend(), static_cast<int>(0),
                                                            std::plus<>{}, [](auto res){return cmc::bits::GetLZC(res);});

            /* If the predictor maximizes the LZC, we store it */
            if (current_lzc < cumulative_lzc)
            {
                current_lzc = cumulative_lzc;
                lzc_maximizing_predictor = predictors[pred_idx];
                std::copy_n(current_residuals.cbegin(), N, residuals.begin());
                std::copy_n(current_entropy_codes.cbegin(), N, entropy_codes.begin());
            }
        }

        return std::make_tuple(lzc_maximizing_predictor, residuals, entropy_codes);
    }
}

template<EightByteArithmeticType T, int32_t N>
inline std::tuple<T, std::array<EightByteResidualType, N>, std::array<SymbolType, N>>
PerformDefaultIntraElemLevelCompression(const std::array<T, N> init_data)
{
    static_assert(N >= 1);

    if constexpr (N == 1)
    {
        /* In case a 2:1 coarsening is used and the overall number is uneven, we just return the same value as extracted */
        return std::make_tuple(init_data[0], std::array<EightByteResidualType, 1>{0}, std::array<SymbolType, 1>{CreateEntropySymbol(false, EightByteResidualType{0})});
    } else
    {
        /* Create all predictors */
        const std::array<T, N+2> predictors = CreatePredictors<T, N>(init_data);

        int current_lzc{-1};
        T lzc_maximizing_predictor{};
        std::array<EightByteResidualType, N> residuals{};
        std::array<SymbolType, N> entropy_codes{};

        for (int pred_idx{0}; pred_idx < N + 2; ++pred_idx)
        {
            std::array<EightByteResidualType, N> current_residuals{};
            std::array<SymbolType, N> current_entropy_codes{};

            /* This computation has a symmetrical part (in those cases an init value is used as a predictor), but we are computing it fully currently */
            /* Check the predictor for all element values */
            for (int elem_idx{0}; elem_idx < N; ++elem_idx)
            {
                /* Compute the residual */
                const auto [is_approx_greater, residual] = cmc::bits::ComputeIntegerResidual<T>(predictors[pred_idx], predictors[elem_idx]);
                
                /* Store the residual */
                current_residuals[elem_idx] = residual;

                /* Store the entropy code */
                current_entropy_codes[elem_idx] = CreateEntropySymbol(is_approx_greater, residual);
            }

            /* Compute the cumulative leading zero count */
            const int cumulative_lzc = std::transform_reduce(std::execution::par_unseq, current_residuals.cbegin(), current_residuals.cend(), static_cast<int>(0),
                                                            std::plus<>{}, [](auto res){return cmc::bits::GetLZC(res);});

            /* If the predictor maximizes the LZC, we store it */
            if (current_lzc < cumulative_lzc)
            {
                current_lzc = cumulative_lzc;
                lzc_maximizing_predictor = predictors[pred_idx];
                std::copy_n(current_residuals.cbegin(), N, residuals.begin());
                std::copy_n(current_entropy_codes.cbegin(), N, entropy_codes.begin());
            }
        }

        return std::make_tuple(lzc_maximizing_predictor, residuals, entropy_codes);
    }
}

/** Gather the values in a linear fashion from the data based on the extraction key defined by kPackSize<DIM> */
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline std::array<T, kPackSize<DIM>>
GetValues(const std::span<T> data, const int offset)
{
    static_assert(kPackSize<DIM> > 1, "The kPackSize<DIM> for this dimension needs to be larger than one in order to gather the values for an extraction.");
    static_assert(kPackSize<DIM> <= N, "The kPackSize<DIM> cannot be larger than the amount of data N per element.");

    if constexpr (kPackSize<DIM> == 2)
    {
        return std::array<T, kPackSize<DIM>>{data[offset], data[offset + 1]};
    } else if constexpr (kPackSize<DIM> == 4)
    {
        return std::array<T, kPackSize<DIM>>{data[offset], data[offset + 1], data[offset + 2], data[offset + 3]};
    } else if constexpr (kPackSize<DIM> == 8)
    {
        return std::array<T, kPackSize<DIM>>{data[offset], data[offset + 1], data[offset + 2], data[offset + 3],
                                             data[offset + 4], data[offset + 5], data[offset + 6], data[offset + 7]};
    } else if constexpr (kPackSize<DIM> == 16)
    {
        return std::array<T, kPackSize<DIM>>{data[offset], data[offset + 1], data[offset + 2], data[offset + 3],
                                             data[offset + 4], data[offset + 5], data[offset + 6], data[offset + 7],
                                             data[offset + 8], data[offset + 9], data[offset + 10], data[offset + 11],
                                             data[offset + 12], data[offset + 13], data[offset + 14], data[offset + 15]};
    } else
    {
        std::array<T, kPackSize<DIM>> values{};
        std::copy_n(&data[offset], kPackSize<DIM>, values.data());
        return values;
    }
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
IntraElementCoding<T, DIM, N>
ComputeElementEncoding(const std::span<T> init_data)
{
    /* We always coarsen kPackSize values at a time */
	constexpr int num_compression_lvls = ComputeNumIntraCompressionLevels<DIM, N>();

    int current_start_lvl_idx{0};
    int current_end_lvl_idx{0};
    int access_idx{0};

    /* Allocate intra-element predictors */
    std::array<T, kNumIntraPredictors<DIM, N>> predictors{};
    int pred_idx{0};

    /* Allocate the intra element coding struct */
	IntraElementCoding<T, DIM, N> elem_coding;

    /* Iterate over all compression levels */
    for_constexpr<0, num_compression_lvls, 1>([&](auto LVL_IDX)
    {
        /* During the first iteration, we need to get the data from the span */
        if constexpr (LVL_IDX == 0)
        {
            /* Determine the number of full packs */
            constexpr int num_full_packs = ComputeNumFullPacks<DIM, N>(LVL_IDX);

            /* Iterate over all full packs */
            for (int pack_idx{0}; pack_idx < num_full_packs; ++pack_idx)
            {
                /* Determine current start index */
                const int elem_idx = kPackSize<DIM> * pack_idx;

                /* Encode this sub-element */
                const auto [predictor, residuals, entropy_codes] = PerformDefaultIntraElemLevelCompression<T, kPackSize<DIM>>(GetValues<T, DIM, N>(init_data, elem_idx));

                /* Store the predictor */
                predictors[pack_idx] = predictor;

                /* Store the residuals */
                std::copy_n(residuals.data(), kPackSize<DIM>, &(elem_coding.residuals[elem_idx]));                                

                /* Store the entropy codes */
                std::copy_n(entropy_codes.data(), kPackSize<DIM>, &(elem_coding.entropy_codes[elem_idx]));
            }

            /* Potentially, compute the encoding of an incomplete pack */
            constexpr int num_elems_incomplete_pack = ComputeIntraCompressionLevelIncompletePackSize<DIM, N>(LVL_IDX);
            if constexpr (num_elems_incomplete_pack > 0)
            {
                /* Compute coarsening for the incomplete pack */
                std::array<T, num_elems_incomplete_pack> incomplete_pack_data;
                std::copy_n(init_data.data() + num_full_packs * kPackSize<DIM>, num_elems_incomplete_pack, incomplete_pack_data.data());
                
                /* Encode this incomplete sub-element */
                const auto [predictor, residuals, entropy_codes] = PerformDefaultIntraElemLevelCompression<T, num_elems_incomplete_pack>(incomplete_pack_data);
                
                /* Store the predictor */
                predictors[num_full_packs] = predictor;
                
                /* Store the residuals */
                std::copy_n(residuals.data(), num_elems_incomplete_pack, &(elem_coding.residuals[num_full_packs * kPackSize<DIM>]));                                

                /* Store the entropy codes from the incomplete pack */
                std::copy_n(entropy_codes.data(), num_elems_incomplete_pack, &(elem_coding.entropy_codes[num_full_packs * kPackSize<DIM>]));
            }

            /* Compute the offsets */
            pred_idx = num_full_packs + (num_elems_incomplete_pack > 0 ? 1 : 0);
            access_idx = num_full_packs * kPackSize<DIM> + num_elems_incomplete_pack;

            /* Update the end counter for the next iteration level */
            current_end_lvl_idx = pred_idx;
        }
        /* In all succeeding iterations, we need the recently computed data */
        else
        {
            /* Compute the number of packs that are full concerning a coarsening */
            constexpr int num_full_packs = ComputeNumFullPacks<DIM, N>(LVL_IDX);

            /* Coarsen the full packs of data */
            for (int pack_idx{0}; pack_idx < num_full_packs; ++pack_idx, ++pred_idx, access_idx += kPackSize<DIM>)
            {
                /* Compute the current start position */
                const int elem_idx = current_start_lvl_idx + kPackSize<DIM> * pack_idx;

                /* Create a view on the values to be coarsened*/
                const std::span<T> pred_values(predictors.data(), pred_idx);

                /* Encode this sub-element */
                const auto [predictor, residuals, entropy_codes] = PerformDefaultIntraElemLevelCompression<T, kPackSize<DIM>>(GetValues<T, DIM, N>(pred_values, elem_idx));

                /* Store the predictor */
                predictors[pred_idx] = predictor;

                /* Store the residuals */
                std::copy_n(residuals.data(), kPackSize<DIM>, &(elem_coding.residuals[access_idx]));                                

                /* Store the entropy codes*/
                std::copy_n(entropy_codes.data(), kPackSize<DIM>, &(elem_coding.entropy_codes[access_idx]));
            }

            /* Handle potential incomplete pack */
            constexpr int num_elems_incomplete_pack = ComputeIntraCompressionLevelIncompletePackSize<DIM, N>(LVL_IDX);
            
            if constexpr (num_elems_incomplete_pack > 0)
            {
                /* Get the data of the incomplete pack */
                std::array<T, num_elems_incomplete_pack> incomplete_pack_data;
                std::copy_n(predictors.data() + current_start_lvl_idx + num_full_packs * kPackSize<DIM>, num_elems_incomplete_pack, incomplete_pack_data.data());
                
                /* Encode this incomplete sub-element */
                const auto [predictor, residuals, entropy_codes] = PerformDefaultIntraElemLevelCompression<T, num_elems_incomplete_pack>(incomplete_pack_data);
                
                /* Store the predictor */
                predictors[pred_idx] = predictor;
                ++pred_idx;
                
                /* Store the residuals */
                std::copy_n(residuals.data(), num_elems_incomplete_pack, &(elem_coding.residuals[access_idx]));                                

                /* Store the entropy codes */
                std::copy_n(entropy_codes.data(), num_elems_incomplete_pack, &(elem_coding.entropy_codes[access_idx]));
                
                /* Update the accessor */
                access_idx += num_elems_incomplete_pack;
            }
            
            /* Update the start and end for the next compression level */
            current_start_lvl_idx = current_end_lvl_idx;
            current_end_lvl_idx = pred_idx;
        }
    });

    /* Set the coarse level predictor */
    elem_coding.coarse_value_predictor = predictors.back();

    /* Return the element encoding specifications */
    return elem_coding;
}

/**
 * We store the computed encoding level-wise in reverse (from caorse to fine) and encode all entropy codes and afterwards all residuals
 */
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline void
PerformElementEncoding(cmc::bits::vector& encoding, const cmc::entropy_coding::huffman::HuffmanCoder<SymbolType>& huffman_coder,
                       const IntraElementCoding<T, DIM, N>& elem_coding)
{
    /* We always coarsen kPackSize values at a time */
	constexpr int num_compression_lvls = ComputeNumIntraCompressionLevels<DIM, N>();

    /* Get the number of values per intra compression level */
    constexpr std::array<int32_t, num_compression_lvls> num_elems_per_level = ComputeNumCodesPerLevel<DIM, N>();

    /* Get the overall number of entropy codes/residuals */
    constexpr int32_t num_all_codes = ComputeNumPyramidalCodes<DIM, N>();
    static_assert(num_all_codes >= 1);

    int32_t offset{num_all_codes};

    /* Iterate over all compression levels */
    for_constexpr<0, num_compression_lvls, 1>([&](auto LVL_IDX)
    {
        /* Store the data level-wise in reverse */
        constexpr int32_t lvl_idx = num_compression_lvls - 1 - LVL_IDX;

        /* Get the number of elements on this level */
        constexpr int32_t num_elems_on_level = num_elems_per_level[lvl_idx];

        /* Compute the encoding offset */
        offset -= num_elems_on_level;

        /* Get a view on the entropy codes */
        const std::span<const SymbolType> entropy_codes(elem_coding.entropy_codes.begin() + offset, num_elems_on_level);

        /* Iterate over the symbols and residuals on this level and encode and append them to the encoded bits::vector */
        for (int32_t idx{0}; idx < num_elems_on_level; ++idx)
        {
            /* Encode all entropy codes */
            const cmc::entropy_coding::huffman::HuffmanCode code = huffman_coder.EncodeSymbol(entropy_codes[idx]);

            /* Store the bit sequence of the code */
            encoding.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
        }

        for (int32_t idx{0}; idx < num_elems_on_level; ++idx)
        {
            /* Store all significant bits of the residuals, we do not need to store the implicit one bit that succeeds the leading zeros */
            const int lzc = GetLZCFromEntropySymbol(entropy_codes[idx]);

            if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
            {
                encoding.AppendBits(elem_coding.residuals[offset + idx], lzc + 1, 0);
            }
        }
    });
}

template<FourByteArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::array<T, N>
PerformElementDecoding(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const T coarse_predictor)
{
    /* Allocate memory for the predictors */
    std::array<FourByteResidualType, kNumIntraPredictors<DIM, N> + N> predictors{};

    /* Store the coarse level predictor */
    predictors[0] = std::bit_cast<FourByteResidualType>(coarse_predictor);
    int pred_idx{1};

    /* We always coarsen kPackSize values at a time */
	constexpr int num_compression_lvls = ComputeNumIntraCompressionLevels<DIM, N>();

    /* Get the number of values per intra compression level */
    constexpr std::array<int32_t, num_compression_lvls> num_elems_per_level = ComputeNumCodesPerLevel<DIM, N>();

    /* Get the overall number of entropy codes/residuals */
    constexpr int32_t num_all_codes = ComputeNumPyramidalCodes<DIM, N>();
    static_assert(num_all_codes >= 1);

    int32_t offset{num_all_codes};
    int pred_lvl_offset{0};

    /* Iterate over all compression levels */
    for_constexpr<0, num_compression_lvls, 1>([&](auto LVL_IDX)
    {
        /* Store the data level-wise in reverse */
        constexpr int32_t lvl_idx = num_compression_lvls - 1 - LVL_IDX;

        /* Get the number of elements on this level */
        constexpr int32_t num_elems_on_level = num_elems_per_level[lvl_idx];

        /* Compute the encoding offset */
        offset -= num_elems_on_level;

        /* Allocate all entropy codes on this level */
        std::array<SymbolType, num_elems_on_level> lvl_entropy_codes;
        
        /* Gather all entropy codes for this level */
        for (int32_t idx{0}; idx < num_elems_on_level; ++idx)
        {
            /* Decode the next entropy symbol */
            lvl_entropy_codes[idx] = stream_decoder.DecodeNextEntropySymbol();
        }

        /* Determine the number of full packs */
        constexpr int num_full_packs = ComputeNumFullPacks<DIM, N>(lvl_idx);

        /* Iterate over all full packs */
        for (int pack_idx{0}; pack_idx < num_full_packs; ++pack_idx)
        {
            /* Get the predictor for this pack of data on this level */
            const FourByteResidualType predictor = predictors[pred_lvl_offset + pack_idx];

            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = pack_idx * kPackSize<DIM>;

            /* Get all residuals on this level and apply them */
            for (int idx{0}; idx < kPackSize<DIM>; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const FourByteResidualType residual = stream_decoder.GetNextBitSequence<FourByteResidualType>(residual_length) | (FourByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const FourByteResidualType residual = (FourByteResidualType{0x70000000} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Potentially, handle incomplete pack on this level */
        constexpr int num_elems_incomplete_pack = ComputeIntraCompressionLevelIncompletePackSize<DIM, N>(lvl_idx);

        if constexpr (num_elems_incomplete_pack > 0)
        {
            /* Get the predictor for this pack of data on this level */
            const FourByteResidualType predictor = predictors[pred_lvl_offset + num_full_packs];
            
            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = num_full_packs * kPackSize<DIM>;
            
            for (int idx{0}; idx < num_elems_incomplete_pack; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const FourByteResidualType residual = stream_decoder.GetNextBitSequence<FourByteResidualType>(residual_length) | (FourByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const FourByteResidualType residual = (FourByteResidualType{0x70000000} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const FourByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Update the offset for the next level predictors */
        pred_lvl_offset += num_full_packs + (num_elems_incomplete_pack > 0 ? 1 : 0);
    });

    /* Extract the the leaf data */
    std::array<T, N> decompressed_elem_data;
    for (int idx{0}; idx < N; ++idx)
    {
        decompressed_elem_data[idx] = std::bit_cast<T>(predictors[kNumIntraPredictors<DIM, N> + idx]);
    }

    return decompressed_elem_data;
}

template<EightByteArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::array<T, N>
PerformElementDecoding(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const T coarse_predictor)
{
    /* Allocate memory for the predictors */
    std::array<EightByteResidualType, kNumIntraPredictors<DIM, N> + N> predictors{};

    /* Store the coarse level predictor */
    predictors[0] = std::bit_cast<EightByteResidualType>(coarse_predictor);
    int pred_idx{1};

    /* We always coarsen kPackSize values at a time */
	constexpr int num_compression_lvls = ComputeNumIntraCompressionLevels<DIM, N>();

    /* Get the number of values per intra compression level */
    constexpr std::array<int32_t, num_compression_lvls> num_elems_per_level = ComputeNumCodesPerLevel<DIM, N>();

    /* Get the overall number of entropy codes/residuals */
    constexpr int32_t num_all_codes = ComputeNumPyramidalCodes<DIM, N>();
    static_assert(num_all_codes >= 1);

    int32_t offset{num_all_codes};
    int pred_lvl_offset{0};

    /* Iterate over all compression levels */
    for_constexpr<0, num_compression_lvls, 1>([&](auto LVL_IDX)
    {
        /* Store the data level-wise in reverse */
        constexpr int32_t lvl_idx = num_compression_lvls - 1 - LVL_IDX;

        /* Get the number of elements on this level */
        constexpr int32_t num_elems_on_level = num_elems_per_level[lvl_idx];

        /* Compute the encoding offset */
        offset -= num_elems_on_level;

        /* Allocate all entropy codes on this level */
        std::array<SymbolType, num_elems_on_level> lvl_entropy_codes;
        
        /* Gather all entropy codes for this level */
        for (int32_t idx{0}; idx < num_elems_on_level; ++idx)
        {
            /* Decode the next entropy symbol */
            lvl_entropy_codes[idx] = stream_decoder.DecodeNextEntropySymbol();
        }

        /* Determine the number of full packs */
        constexpr int num_full_packs = ComputeNumFullPacks<DIM, N>(lvl_idx);

        /* Iterate over all full packs */
        for (int pack_idx{0}; pack_idx < num_full_packs; ++pack_idx)
        {
            /* Get the predictor for this pack of data on this level */
            const EightByteResidualType predictor = predictors[pred_lvl_offset + pack_idx];

            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = pack_idx * kPackSize<DIM>;

            /* Get all residuals on this level and apply them */
            for (int idx{0}; idx < kPackSize<DIM>; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const EightByteResidualType residual = stream_decoder.GetNextBitSequence<EightByteResidualType>(residual_length) | (EightByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const EightByteResidualType residual = (EightByteResidualType{0x7000000000000000} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Potentially, handle incomplete pack on this level */
        constexpr int num_elems_incomplete_pack = ComputeIntraCompressionLevelIncompletePackSize<DIM, N>(lvl_idx);
        if constexpr (num_elems_incomplete_pack > 0)
        {
            /* Get the predictor for this pack of data on this level */
            const EightByteResidualType predictor = predictors[pred_lvl_offset + num_full_packs];
            
            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = num_full_packs * kPackSize<DIM>;

            for (int idx{0}; idx < num_elems_incomplete_pack; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const EightByteResidualType residual = stream_decoder.GetNextBitSequence<EightByteResidualType>(residual_length) | (EightByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const EightByteResidualType residual = (EightByteResidualType{0x7000000000000000} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const EightByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Update the offset for the next level predictors */
        pred_lvl_offset += num_full_packs + (num_elems_incomplete_pack > 0 ? 1 : 0);
    });

    /* Extract the the leaf data */
    std::array<T, N> decompressed_elem_data;
    for (int idx{0}; idx < N; ++idx)
    {
        decompressed_elem_data[idx] = std::bit_cast<T>(predictors[kNumIntraPredictors<DIM, N> + idx]);
    }

    return decompressed_elem_data;
}

template<TwoByteArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::array<T, N>
PerformElementDecoding(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const T coarse_predictor)
{
    /* Allocate memory for the predictors */
    std::array<TwoByteResidualType, kNumIntraPredictors<DIM, N> + N> predictors{};

    /* Store the coarse level predictor */
    predictors[0] = std::bit_cast<TwoByteResidualType>(coarse_predictor);
    int pred_idx{1};

    /* We always coarsen kPackSize values at a time */
	constexpr int num_compression_lvls = ComputeNumIntraCompressionLevels<DIM, N>();

    /* Get the number of values per intra compression level */
    constexpr std::array<int32_t, num_compression_lvls> num_elems_per_level = ComputeNumCodesPerLevel<DIM, N>();

    /* Get the overall number of entropy codes/residuals */
    constexpr int32_t num_all_codes = ComputeNumPyramidalCodes<DIM, N>();
    static_assert(num_all_codes >= 1);

    int32_t offset{num_all_codes};
    int pred_lvl_offset{0};

    /* Iterate over all compression levels */
    for_constexpr<0, num_compression_lvls, 1>([&](auto LVL_IDX)
    {
        /* Store the data level-wise in reverse */
        constexpr int32_t lvl_idx = num_compression_lvls - 1 - LVL_IDX;

        /* Get the number of elements on this level */
        constexpr int32_t num_elems_on_level = num_elems_per_level[lvl_idx];

        /* Compute the encoding offset */
        offset -= num_elems_on_level;

        /* Allocate all entropy codes on this level */
        std::array<SymbolType, num_elems_on_level> lvl_entropy_codes;
        
        /* Gather all entropy codes for this level */
        for (int32_t idx{0}; idx < num_elems_on_level; ++idx)
        {
            /* Decode the next entropy symbol */
            lvl_entropy_codes[idx] = stream_decoder.DecodeNextEntropySymbol();
        }

        /* Determine the number of full packs */
        constexpr int num_full_packs = ComputeNumFullPacks<DIM, N>(lvl_idx);

        /* Iterate over all full packs */
        for (int pack_idx{0}; pack_idx < num_full_packs; ++pack_idx)
        {
            /* Get the predictor for this pack of data on this level */
            const TwoByteResidualType predictor = predictors[pred_lvl_offset + pack_idx];

            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = pack_idx * kPackSize<DIM>;

            /* Get all residuals on this level and apply them */
            for (int idx{0}; idx < kPackSize<DIM>; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const TwoByteResidualType residual = stream_decoder.GetNextBitSequence<TwoByteResidualType>(residual_length) | (TwoByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const TwoByteResidualType residual = (TwoByteResidualType{0x7000} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Potentially, handle incomplete pack on this level */
        constexpr int num_elems_incomplete_pack = ComputeIntraCompressionLevelIncompletePackSize<DIM, N>(lvl_idx);
        if constexpr (num_elems_incomplete_pack > 0)
        {
            /* Get the predictor for this pack of data on this level */
            const TwoByteResidualType predictor = predictors[pred_lvl_offset + num_full_packs];
            
            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = num_full_packs * kPackSize<DIM>;

            for (int idx{0}; idx < num_elems_incomplete_pack; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const TwoByteResidualType residual = stream_decoder.GetNextBitSequence<TwoByteResidualType>(residual_length) | (TwoByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const TwoByteResidualType residual = (TwoByteResidualType{0x7000} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const TwoByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Update the offset for the next level predictors */
        pred_lvl_offset += num_full_packs + (num_elems_incomplete_pack > 0 ? 1 : 0);
    });

    /* Extract the the leaf data */
    std::array<T, N> decompressed_elem_data;
    for (int idx{0}; idx < N; ++idx)
    {
        decompressed_elem_data[idx] = std::bit_cast<T>(predictors[kNumIntraPredictors<DIM, N> + idx]);
    }

    return decompressed_elem_data;
}

template<OneByteArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::array<T, N>
PerformElementDecoding(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const T coarse_predictor)
{
    /* Allocate memory for the predictors */
    std::array<OneByteResidualType, kNumIntraPredictors<DIM, N> + N> predictors{};

    /* Store the coarse level predictor */
    predictors[0] = std::bit_cast<OneByteResidualType>(coarse_predictor);
    int pred_idx{1};

    /* We always coarsen kPackSize values at a time */
	constexpr int num_compression_lvls = ComputeNumIntraCompressionLevels<DIM, N>();

    /* Get the number of values per intra compression level */
    constexpr std::array<int32_t, num_compression_lvls> num_elems_per_level = ComputeNumCodesPerLevel<DIM, N>();

    /* Get the overall number of entropy codes/residuals */
    constexpr int32_t num_all_codes = ComputeNumPyramidalCodes<DIM, N>();
    static_assert(num_all_codes >= 1);

    int32_t offset{num_all_codes};
    int pred_lvl_offset{0};

    /* Iterate over all compression levels */
    for_constexpr<0, num_compression_lvls, 1>([&](auto LVL_IDX)
    {
        /* Store the data level-wise in reverse */
        constexpr int32_t lvl_idx = num_compression_lvls - 1 - LVL_IDX;

        /* Get the number of elements on this level */
        constexpr int32_t num_elems_on_level = num_elems_per_level[lvl_idx];

        /* Compute the encoding offset */
        offset -= num_elems_on_level;

        /* Allocate all entropy codes on this level */
        std::array<SymbolType, num_elems_on_level> lvl_entropy_codes;
        
        /* Gather all entropy codes for this level */
        for (int32_t idx{0}; idx < num_elems_on_level; ++idx)
        {
            /* Decode the next entropy symbol */
            lvl_entropy_codes[idx] = stream_decoder.DecodeNextEntropySymbol();
        }

        /* Determine the number of full packs */
        constexpr int num_full_packs = ComputeNumFullPacks<DIM, N>(lvl_idx);

        /* Iterate over all full packs */
        for (int pack_idx{0}; pack_idx < num_full_packs; ++pack_idx)
        {
            /* Get the predictor for this pack of data on this level */
            const OneByteResidualType predictor = predictors[pred_lvl_offset + pack_idx];

            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = pack_idx * kPackSize<DIM>;

            /* Get all residuals on this level and apply them */
            for (int idx{0}; idx < kPackSize<DIM>; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const OneByteResidualType residual = stream_decoder.GetNextBitSequence<OneByteResidualType>(residual_length) | (OneByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const OneByteResidualType residual = (OneByteResidualType{0x70} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Potentially, handle incomplete pack on this level */
        constexpr int num_elems_incomplete_pack = ComputeIntraCompressionLevelIncompletePackSize<DIM, N>(lvl_idx);
        if constexpr (num_elems_incomplete_pack > 0)
        {
            /* Get the predictor for this pack of data on this level */
            const OneByteResidualType predictor = predictors[pred_lvl_offset + num_full_packs];

            /* Compute the offset in this level's entropy codes */
            const int entropy_idx_offset = num_full_packs * kPackSize<DIM>;

            for (int idx{0}; idx < num_elems_incomplete_pack; ++idx)
            {
                /* Determine the leading zero count */
                const int lzc = GetLZCFromEntropySymbol(lvl_entropy_codes[entropy_idx_offset + idx]);

                if (lzc < sizeof(T) * cmc::bits::kCharBit - 1) [[likely]]
                {
                    /* Compute the length of the significant residual bits */
                    const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

                    /* We obtain the residual and add the implicit one bit  */
                    const OneByteResidualType residual = stream_decoder.GetNextBitSequence<OneByteResidualType>(residual_length) | (OneByteResidualType{1} << residual_length);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                } else
                {
                    /* Compute the residual in case we do not need to extract a bit-sequence */
                    const OneByteResidualType residual = (OneByteResidualType{0x70} >> lzc);

                    /* We create the residual applied value */
                    if (IsApproximationGreater(lvl_entropy_codes[entropy_idx_offset + idx]))
                    {
                        /* Subtract the residual from the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerSubtraction(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    } else
                    {
                        /* Add the residual to the prediction */
                        const OneByteResidualType value = cmc::bits::IntegerAddition(predictor, residual);

                        /* Store the residual applied value */
                        predictors[pred_idx] = value;
                        ++pred_idx;
                    }
                }
            }
        }

        /* Update the offset for the next level predictors */
        pred_lvl_offset += num_full_packs + (num_elems_incomplete_pack > 0 ? 1 : 0);
    });

    /* Extract the the leaf data */
    std::array<T, N> decompressed_elem_data;
    for (int idx{0}; idx < N; ++idx)
    {
        decompressed_elem_data[idx] = std::bit_cast<T>(predictors[kNumIntraPredictors<DIM, N> + idx]);
    }

    return decompressed_elem_data;
}

/***** Specialized implementations for certain setups *****/
#if 0
/*** Start: Dimension: 2; Num Points per Element: 4 ***/
template<ArithmeticType T, int32_t DIM, int32_t N>
requires (DIM == 2 && N == 4)
IntraElementCoding<T, 2, 4>
ComputeElementEncoding(const std::span<T> init_data)
{
	/* Implementation of 2D elements with four data points per element */
	static_assert(false, "This implementation is not yet available");
    //...
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires (DIM == 2 && N == 4)
inline void
PerformElementEncoding(cmc::bits::vector& encoding, const cmc::entropy_coding::huffman::HuffmanCoder<SymbolType>& huffman_coder,
                       const IntraElementCoding<T, D, N>& elem_coding)
{
    /* Implementation of 2D elements with four data points per element */
	static_assert(false, "This implementation is not yet available");
    //...
}
/*** End: Dimension: 2; Num Points per Element: 4 ***/
#endif
/**********************************************************/

}

#endif /* !CMC_PAR_MULTI_RES_EXTRACTION_UTIL_HXX */
