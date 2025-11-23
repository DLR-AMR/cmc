#ifndef CMC_LOSSLESS_PATCH_BYTE_COMPRESSION_VARIABLE_HXX
#define CMC_LOSSLESS_PATCH_BYTE_COMPRESSION_VARIABLE_HXX

#include "cmc.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_iface_patch_compression_variable.hxx"

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <type_traits>

namespace cmc::patch::lossless
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

constexpr size_t kDimReductionFactor = 2;

/**
 * @brief A struct holding the data for an extraction process.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct ExtractionData
{
    ExtractionData(const CompressionValue<T>& coarse_val, std::vector<CompressionValue<T>>&& fine_vals)
    : coarse_value(coarse_val), fine_values(std::move(fine_vals)) {};
    ExtractionData(CompressionValue<T>&& coarse_val, std::vector<CompressionValue<T>>&& fine_vals)
    : coarse_value(std::move(coarse_val)), fine_values(std::move(fine_vals)) {};
    
    CompressionValue<T> coarse_value;
    std::vector<CompressionValue<T>> fine_values;
};

struct kTag3D {};

/**
 * @brief The Interface/Template for a variable that performs lossless compression on the serialized data
 * in a byte-/bit-wise fashion. The compression algorithm is fixed and may be specialized in a derived class
 * in order to fit the compression for the given needs.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T, size_t Dim>
class AbstractPatchByteCompressionVariable : public IPatchCompressionVariable<T>
{
public:
    void Compress() override;

    virtual ~AbstractPatchByteCompressionVariable(){}

    void MoveEncodedEntropyCodesInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_entropy_codes) override;
    void MoveEncodedDataInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_data) override;

    const std::vector<std::vector<uint8_t>>& GetEncodedEntropyCodes() const override {return buffered_entropy_codes_;}
    const std::vector<std::vector<uint8_t>>& GetEncodedData() const override {return buffered_encoded_data_;}

    const std::string& GetName() const override {return name_;}
    size_t Size() const override {return data_.size();}
    const std::vector<CompressionValue<T>>& GetData() const {return data_;}
    size_t GetDimensionality() const override {return Dim;}
    const std::vector<size_t>& GetInitialDimensionLengths() const override {return dim_lengths_pyramid_.front();}
    const std::vector<std::vector<size_t>>& GetDimensionLengthPyramid() const override {return dim_lengths_pyramid_;}
    DataLayout GetInitialDataLayout() const override {return init_data_layout_;}
    GeoDomain GetInitialDomain() const override {return init_domain_;}
    int GetNumCompressionIterations() const {return num_compression_lvls_;}

    virtual CompressionSchema GetCompressionSchema() const = 0;

    static struct kTag3D{} tag3D;
    static struct kTag2D{} tag2D;
    static struct kTag1D{} tag1D;
protected:
    AbstractPatchByteCompressionVariable() = delete;
    AbstractPatchByteCompressionVariable(input::Var& input_variable)
    {
        this->SetupInputVariableForCompression(input_variable);
    };

    virtual void PreCompressionProcessing([[maybe_unused]] std::vector<CompressionValue<T>>& initial_data){}
    virtual void InitializeExtractionIteration() {}
    virtual void FinalizeExtractionIteration() {}
    virtual void CompleteExtraction([[maybe_unused]] const std::vector<size_t>& fine_dim_lengths, [[maybe_unused]] const size_t kDimReductionFactor, [[maybe_unused]] const std::vector<CompressionValue<T>>& fine_vals, [[maybe_unused]] const std::vector<CompressionValue<T>>& coarse_vals) {}
    virtual ExtractionData<T> PerformExtraction(const std::vector<CompressionValue<T>>& patch) = 0;
    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const = 0;
    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const = 0;

    void SetName(const std::string& name) {name_ = name;}
    void SetData(const std::vector<T>& initial_data);
    void SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data);
    void SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data);
    
private:
    void Compress(kTag3D);
    void SetupInputVariableForCompression(input::Var& input_variable);

    std::string name_; //!< The name of the variable
    std::vector<CompressionValue<T>> data_; //!< The current data of the variable 
    std::vector<CompressionValue<T>> data_new_; //!< A helper variable for the adaptation

    std::vector<std::vector<uint8_t>> buffered_entropy_codes_; //!< Level-wise storage of the entropy codes, e.g. LZC or prefix lengths
    std::vector<std::vector<uint8_t>> buffered_encoded_data_; //!< Level-wise storage for the encoded data
    std::vector<std::vector<size_t>> dim_lengths_pyramid_;

    GeoDomain init_domain_;
    DataLayout init_data_layout_;
    int num_compression_lvls_{0};
};


template <typename T, size_t Dim>
void
AbstractPatchByteCompressionVariable<T, Dim>::Compress()
{
    if constexpr (Dim == 3)
    {
        cmc_debug_msg("Lossless 3D Compression.");
        this->Compress(AbstractPatchByteCompressionVariable<T, Dim>::tag3D);
    } else
    {
        cmc_err_msg("Unsupported variable's dimensionality (Dim = ", Dim, ").");
    }
}


template <typename T, size_t Dim>
void
AbstractPatchByteCompressionVariable<T, Dim>::SetupInputVariableForCompression(input::Var& input_var)
{
    cmc_debug_msg("The InputVariable will be set up for compression.");

    //const std::vector<T>& values = std::get<input::Variable<T>>(input_var.GetInternalVariant()).GetDataForReading();
    //FILE* file_out = fopen("direct_output_data_bin_reader.cmc", "wb");
    //fwrite(values.data(), sizeof(T), values.size(), file_out);
    //fclose(file_out);
    //cmc_err_msg("End here");

    /* Potentially, apply scaling values and offsets if defined (potentially, the datatype changes by this call) */
    input_var.ApplyScalingAndOffset();
    
    if (not std::holds_alternative<input::Variable<T>>(input_var.GetInternalVariant()))
    {
        cmc_err_msg("The data type of the input variable differs from the patch-based byte compression variable.");
    }

    /* Get the name of the variable */
    name_ = input_var.GetName();

    /* Get the actual variable from the input var wrapper */
    input::Variable<T> input_variable = std::get<input::Variable<T>>(input_var.GetInternalVariant());

    /* Get the global domain of the variable */
    init_domain_ = input_variable.GetGlobalDomain();

    /* Get the initial data layout */
    init_data_layout_ = input_variable.GetInitialDataLayout();

    cmc_assert(GetDimensionalityOfDataLayout(init_data_layout_) == init_domain_.GetDimensionality());

    /* Check whether the dimensionality of the patch variable is consisten with the init global domain */
    if (static_cast<size_t>(init_domain_.GetDimensionality()) != Dim)
    {
        cmc_err_msg("The dimensionality of the domain does not coincide with the dimensionality of the patch variable.");
    }

    /* Get the dimension id vector of the layout */
    std::vector<Dimension> dims = GetDimensionVectorFromLayout(init_data_layout_);

    /* Store the initial dim_lengths in the pyramid */
    dim_lengths_pyramid_.emplace_back();

    /* Get the length of the dimensions and store them in the pyramid */
    for (auto dim_iter = dims.begin(); dim_iter != dims.end(); ++dim_iter)
    {
        dim_lengths_pyramid_.back().push_back(init_domain_.GetDimensionLength(*dim_iter));
    }

    /* Get the largest dimension length */
    const size_t max_dim_length = static_cast<size_t>(init_domain_.GetLargestDimensionLength());

    /* Determine the number of compression iterations */
    num_compression_lvls_ = static_cast<int>(std::ceil(std::log2(max_dim_length) + std::numeric_limits<double>::epsilon()));

    //FILE* file_out = fopen("initial_input_prefix_amr_data.cmc", "wb");
    //fwrite(input_variable.GetDataForReading().data(), sizeof(T), input_variable.GetDataForReading().size(), file_out);
    //fclose(file_out);

    /* Get the data from the variable and store it as CompressionValues */
    this->SetData(input_variable.GetDataForReading());
    
    cmc_debug_msg("The setup of the InputVariable for compression has been successfull.");
}

template <typename T>
inline CompressionValue<T>
GetValue(const std::vector<CompressionValue<T>>& data, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength,  const int kLevLength)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    return data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon];
}

template <typename T>
inline void
SetValue(std::vector<CompressionValue<T>>& data, const CompressionValue<T>& value, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength,  const int kLevLength)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon] = value;
}

template <typename T, size_t Dim>
void
AbstractPatchByteCompressionVariable<T, Dim>::Compress(AbstractPatchByteCompressionVariable<T, Dim>::kTag3D)
{
    /* Potentially, create a pre-compression processing step */
    this->PreCompressionProcessing(data_);

    cmc_debug_msg("Lossless patch-based compression of variable ", this->name_, " starts...");

    cmc_assert(dim_lengths_pyramid_.size() == 1);
    cmc_assert(dim_lengths_pyramid_.front().size() == 3);

    /* Store the intial dimension sizes */
    size_t lev_length = dim_lengths_pyramid_.front()[0];
    size_t lat_length = dim_lengths_pyramid_.front()[1];
    size_t lon_length = dim_lengths_pyramid_.front()[2];

    cmc_debug_msg("Number of compression iterations to perform: ", num_compression_lvls_);
    
    /* Perform the iterative compression steps up until the root level */
    for (int lvl_idx{0}; lvl_idx < num_compression_lvls_; ++lvl_idx)
    {
        /* Define a reference for the ease of notation */
        const std::vector<size_t>& dim_lengths = dim_lengths_pyramid_.back();

        cmc_assert(dim_lengths_pyramid_.back().size() == static_cast<size_t>(3));
        cmc_debug_msg("A coarsening iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions: ", dim_lengths[0], ", ", dim_lengths[1], ", ", dim_lengths[2]);
        
        /* Intitialize the extraction iteration */
        this->InitializeExtractionIteration();

        /* Allocate the coarse level */
        data_new_.reserve(((dim_lengths[0] / kDimReductionFactor) + 1) * ((dim_lengths[1] / kDimReductionFactor) + 1) * ((dim_lengths[2] / kDimReductionFactor) + 1));

        /* Iterate over patches */
        for (size_t lev = 0; lev < dim_lengths[0]; lev += kDimReductionFactor)
        {
            for (size_t lat = 0; lat < dim_lengths[1]; lat += kDimReductionFactor)
            {
                for (size_t lon = 0; lon < dim_lengths[2]; lon += kDimReductionFactor)
                {
                    /* Allocate avector for the values of the current patch */
                    std::vector<CompressionValue<T>> patch;
                    patch.reserve(3 * kDimReductionFactor);

                    /* Gather the values for this patch */
                    for (size_t lev_idx = 0; lev_idx < kDimReductionFactor; ++lev_idx)
                    {
                        for (size_t lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                        {
                            for (size_t lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                            {
                                if (lev + lev_idx >= dim_lengths[0] || lat + lat_idx >= dim_lengths[1] || lon + lon_idx >= dim_lengths[2])
                                {
                                    continue;
                                } else
                                {
                                    patch.push_back(GetValue<T>(data_, lev + lev_idx, lat + lat_idx, lon + lon_idx, dim_lengths[2], dim_lengths[1], dim_lengths[0]));
                                }
                            }
                        }
                    }

                    /* Perform an extraction operation on this patch of values */
                    const ExtractionData<T> extracted_values = this->PerformExtraction(patch);

                    /* Re-Assign the fine values */
                    int extracted_val_idx{0};
                    for (size_t lev_idx = 0; lev_idx < kDimReductionFactor; ++lev_idx)
                    {
                        for (size_t lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                        {
                            for (size_t lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                            {
                                if (lev + lev_idx >= dim_lengths[0] || lat + lat_idx >= dim_lengths[1] || lon + lon_idx >= dim_lengths[2])
                                {
                                    continue;
                                } else
                                {
                                    /* Update the fine value accordingly */
                                    SetValue<T>(data_, extracted_values.fine_values[extracted_val_idx], lev + lev_idx, lat + lat_idx, lon + lon_idx, dim_lengths[2], dim_lengths[1], dim_lengths[0]);
                                    ++extracted_val_idx;
                                }
                            }
                        }
                    }

                    /* Store the extracted information */
                    data_new_.push_back(extracted_values.coarse_value);
                }
            }
        }

        /* Potentially, perform a step after the coarse extraction */
        this->CompleteExtraction(dim_lengths, kDimReductionFactor, data_, data_new_);

        /* Encode the data from this compression step */
        auto [encoded_entropy_codes, encoded_data] = this->EncodeLevelData(data_);
        buffered_entropy_codes_.push_back(std::move(encoded_entropy_codes));
        buffered_encoded_data_.push_back(std::move(encoded_data));

        /* Switch to the coarser data for the next extraction iteration */
        data_ = std::move(data_new_);
        data_new_.clear();

        /* Finalize the extraction iteration */
        this->FinalizeExtractionIteration();

        /* Update the dimension sizes */
        lev_length = (lev_length % kDimReductionFactor == 0 ? (lev_length / kDimReductionFactor) : (lev_length / kDimReductionFactor) + 1);
        lat_length = (lat_length % kDimReductionFactor == 0 ? (lat_length / kDimReductionFactor) : (lat_length / kDimReductionFactor) + 1);
        lon_length = (lon_length % kDimReductionFactor == 0 ? (lon_length / kDimReductionFactor) : (lon_length / kDimReductionFactor) + 1);
        
        /* Store the dimensions for the next iteration */
        dim_lengths_pyramid_.emplace_back();
        dim_lengths_pyramid_.back().push_back(lev_length);
        dim_lengths_pyramid_.back().push_back(lat_length);
        dim_lengths_pyramid_.back().push_back(lon_length);

        cmc_debug_msg("The coarsening iteration is finished.");
    }

    /* At last, we need to encode the root level data */
    auto [encoded_root_entropy_codes, encoded_root_data] = this->EncodeRootLevelData(data_);
    buffered_entropy_codes_.push_back(std::move(encoded_root_entropy_codes));
    buffered_encoded_data_.push_back(std::move(encoded_root_data));

    cmc_debug_msg("Compression of variable ", this->name_, " is finished.");
}

template <typename T, size_t Dim>
inline void
AbstractPatchByteCompressionVariable<T, Dim>::SetData(const std::vector<T>& initial_data)
{
    data_.reserve(initial_data.size());

    for (auto val_iter = initial_data.begin(); val_iter != initial_data.end(); ++val_iter)
    {
        data_.emplace_back(CompressionValue<T>(*val_iter));
    }
}

template <typename T, size_t Dim>
inline void
AbstractPatchByteCompressionVariable<T, Dim>::SetData(const std::vector<CompressionValue<T>>& initial_data)
{
    data_ = initial_data;
}

template <typename T, size_t Dim>
inline void
AbstractPatchByteCompressionVariable<T, Dim>::SetData(std::vector<CompressionValue<T>>&& initial_data)
{
    data_ = std::move(initial_data);
}


template <typename T, size_t Dim>
inline void
AbstractPatchByteCompressionVariable<T, Dim>::MoveEncodedEntropyCodesInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_entropy_codes)
{
    vec_to_hold_encoded_levelwise_entropy_codes = std::move(buffered_entropy_codes_);
    buffered_entropy_codes_ = std::vector<std::vector<uint8_t>>();
}

template <typename T, size_t Dim>
inline void
AbstractPatchByteCompressionVariable<T, Dim>::MoveEncodedDataInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_data) 
{
    vec_to_hold_encoded_levelwise_data = std::move(buffered_encoded_data_);
    buffered_encoded_data_ = std::vector<std::vector<uint8_t>>();
}

}


#endif /* !CMC_LOSSLESS_PATCH_BYTE_COMPRESSION_VARIABLE_HXX */
