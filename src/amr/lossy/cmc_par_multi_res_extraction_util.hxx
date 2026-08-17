#ifndef CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_UTIL_HXX
#define CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_UTIL_HXX

#include "utilities/cmc_bits.hxx"
#include "utilities/cmc_bits_vector.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"
#include "mpi/cmc_mpi.hxx"
#include "utilities/cmc_error_domain.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <array>
#include <limits>
#include <execution>
#include <cfloat>
#include <string>

namespace cmc::par::lossy::multi_res
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

using VarInfoType = uint32_t;
using SizeType = uint64_t;
constexpr SizeType kNumCharsVariableName = 256;
inline const MPI_Datatype MPI_SIZE_TYPE = MPI_UINT64_T;

constexpr int kMaxPossibleInitialRefinementLevel = 64;

constexpr int kRootRank = 0;
constexpr int kTagMeshEncoding = 1000;

struct PartitionInfo
{
    PartitionInfo() = default;
    PartitionInfo(const uint64_t num_elems_, const uint64_t num_bytes_encoding_)
    : num_elems{num_elems_}, num_bytes_encoding{num_bytes_encoding_} {}
    PartitionInfo(const uint64_t num_elems_)
    : num_elems{num_elems_} {}

    uint64_t num_elems{0};
    uint64_t num_bytes_encoding{0};
};

inline void
CreatePartitionInfoMPIType(MPI_Datatype* partition_info_mpi_type)
{
    /* Define the properties of the custom 'LevelOffset' data type */
    constexpr int num_fields = 2;
    int array_of_blocklengths[] = {1,1};
    MPI_Aint array_of_displacements[num_fields];
    array_of_displacements[0] = offsetof(PartitionInfo, num_elems);
    array_of_displacements[1] = offsetof(PartitionInfo, num_bytes_encoding);
    MPI_Datatype array_of_types[] = {MPI_UINT64_T, MPI_UINT64_T};

    const int ret_val_struct = MPI_Type_create_struct(num_fields, array_of_blocklengths, array_of_displacements,
                                                      array_of_types, partition_info_mpi_type); 
    MPICheckError(ret_val_struct);
    const int ret_val_type_commit = MPI_Type_commit(partition_info_mpi_type);
    MPICheckError(ret_val_type_commit); 
}

struct LevelPartition
{
    LevelPartition() = default;
    LevelPartition(const SizeType elem_offset_, const SizeType coding_byte_offset)
    : elem_offset{elem_offset_}, coding_byte_offset{coding_byte_offset} {}
    
    SizeType elem_offset{0};
    SizeType coding_byte_offset{0};
};

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

/** Entropy Symbol Definitions **/
using SymbolType = uint16_t;

constexpr inline SymbolType kPredictionWithinBound = 0;

constexpr inline int kResidualMaxDeviationFactor = 255;
static_assert(kResidualMaxDeviationFactor % 2 != 0 && kResidualMaxDeviationFactor > 0 && kResidualMaxDeviationFactor < std::numeric_limits<SymbolType>::max() - 2);

constexpr inline int kEntropySymbolBinSignSwitch = (kResidualMaxDeviationFactor - 1) / 2;

constexpr inline SymbolType kFlagUnpredictable = kResidualMaxDeviationFactor;

constexpr inline SymbolType kProcessEndSymbol = kResidualMaxDeviationFactor + 1;

inline SymbolType
CreateEntropySymbolFromQuantizationBin(const SymbolType bin, const bool is_prediction_greater)
{
    cmc_assert(bin > static_cast<SymbolType>(0) && bin <= static_cast<SymbolType>(kEntropySymbolBinSignSwitch));
    return (is_prediction_greater ? kEntropySymbolBinSignSwitch + bin : bin);
}

inline std::pair<bool, int>
GetQuantizationBinFromEntropySymbol(const SymbolType symbol)
{
    cmc_assert(symbol != kProcessEndSymbol && symbol != kFlagUnpredictable);
    const bool is_prediction_greater = (symbol > kEntropySymbolBinSignSwitch);
    const int bin = (is_prediction_greater ? symbol - kEntropySymbolBinSignSwitch : symbol);
    return std::make_pair(is_prediction_greater, bin);
}

template<ArithmeticType T>
constexpr inline int
MapEntropySymbolToArrayIndex(const SymbolType entropy_symbol)
{
    return static_cast<int>(entropy_symbol);
}

template<ArithmeticType T>
constexpr inline SymbolType
MapArrayIndexToEntropySymbol(int array_idx)
{
    return static_cast<SymbolType>(array_idx);
}

constexpr inline int
GetNumEntropySymbols()
{
    return kResidualMaxDeviationFactor + 2;
}

constexpr inline void
AddProcessEndSymbol(std::array<uint64_t, GetNumEntropySymbols()>& entropy_symbols_frequency, const uint64_t num_local_proc_end_symbols)
{
    /* We store the process-end-symbol in the last array entry */
    entropy_symbols_frequency[GetNumEntropySymbols() - 1] = num_local_proc_end_symbols;
}
/** END OF Entropy Symbol Definitions **/


template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
struct ElemEncodingData
{
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> quantization_bins{};
    std::array<T, kNumMaxChildrenElements<DIM>> unpredictable_values{};
    uint32_t num_elements;
};

inline constexpr int32_t kMaxPresentElementLevelUnknown = -1;


template<ArithmeticType T>
inline T
GetAbsResidual(const T& value1, const T& value2)
{
    if constexpr (std::is_signed_v<T>)
    {
        return std::fabs(value2 - value1);
    } else
    {
        return (value1 > value2 ? value1 - value2 : value2 - value1);
    }
}

template <ArithmeticType T>
inline void
WriteDataToVTK(t8_forest_t mesh, const std::vector<T>& data, const std::string variable_name, const std::string file_prefix)
{
    cmc_assert(data.size() >= static_cast<size_t>(t8_forest_get_local_num_leaf_elements(mesh)));
    std::vector<double> double_data;
    double_data.reserve(data.size());

    for (int idx{0}; idx < t8_forest_get_local_num_leaf_elements(mesh); ++idx)
    {
        double_data.push_back(static_cast<double>(data[idx]));
    }

    t8_vtk_data_field_t vtk_data[1];
    snprintf (vtk_data[0].description, BUFSIZ, "%s", variable_name.c_str());
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = double_data.data();

    t8_forest_write_vtk_ext (mesh, file_prefix.c_str(), 1, 1, 1, 1, 0, 0, 0, 1, vtk_data);
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

#endif /* !CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_UTIL_HXX */
