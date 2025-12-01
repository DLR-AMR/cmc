#ifndef CMC_PATCH_BYTE_DECOMPRESSION_VARIABLE_HXX
#define CMC_PATCH_BYTE_DECOMPRESSION_VARIABLE_HXX

#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_entropy_coder.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_iface_patch_decompression_variable.hxx"

namespace cmc::patch::decompression
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

constexpr size_t kDimReductionFactor = 2;

/**
 * @brief The Interface/Template for a variable that performs lossless compression on the serialized data
 * in a byte-/bit-wise fashion. The compression algorithm is fixed and may be specialized in a derived class
 * in order to fit the compression for the given needs.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T, size_t Dim>
class AbstractPatchByteDecompressionVariable : public IPatchDecompressionVariable<T>
{
public:
    void Decompress() override;
    
    virtual ~AbstractPatchByteDecompressionVariable(){}

    std::vector<T> GetDecompressedData() const override;
    const std::string& GetName() const override {return name_;}
    size_t Size() const override {return data_.size();}
    size_t GetDimensionality() const {return Dim;}
    const std::vector<size_t>& GetInitialDimensionLengths() const override {return dim_lengths_pyramid_.back();}
    const std::vector<std::vector<size_t>>& GetDimensionLengthPyramid() const override {return dim_lengths_pyramid_;}
    DataLayout GetInitialDataLayout() const {return init_data_layout_;}
    GeoDomain GetInitialDomain() const override {return init_domain_;}
    int GetNumDecompressionIterations() const {return num_decompression_lvls_;}

    static struct kTag4D{} tag4D;
    static struct kTag3D{} tag3D;
    static struct kTag2D{} tag2D;
    static struct kTag1D{} tag1D;

protected:
    AbstractPatchByteDecompressionVariable() = delete;
    AbstractPatchByteDecompressionVariable(std::vector<uint8_t>&& encoded_data_byte_stream, const std::string& name, std::vector<std::vector<size_t>>&& dim_lengths_pyramid,
                                              const GeoDomain& init_domain, const DataLayout init_data_layout, const int num_decompression_lvls)
    : encoded_data_byte_stream_(std::move(encoded_data_byte_stream)), name_(name), dim_lengths_pyramid_(std::move(dim_lengths_pyramid)),
      init_domain_(init_domain), init_data_layout_{init_data_layout}, num_decompression_lvls_{num_decompression_lvls}
    {
        this->SetupDecompression();
    }

    void SetName(const std::string& name) {name_ = name;};
    void SetInitDataValue(const CompressionValue<T>& init_value);
    const std::vector<uint8_t>& GetEncodedDataStream() const {return encoded_data_byte_stream_;}

    virtual CompressionValue<T> GetInitDecompressionValue() = 0;
    virtual CompressionValue<T> RefineNextValue(const CompressionValue<T>& coarse_value) = 0;
    virtual void InitializeDecompressionIteration() {}
    virtual void FinalizeDecompressionIteration() {}
    
private:
    void SetupDecompression();
    void Decompress(kTag1D);
    void Decompress(kTag2D);
    void Decompress(kTag3D);
    void Decompress(kTag4D);

    std::vector<CompressionValue<T>> data_; //!< The current data of the variable 
    std::vector<CompressionValue<T>> data_new_; //!< A helper variable for the adaptation

    const std::vector<uint8_t> encoded_data_byte_stream_; //!< The encoded byte stream of the variable
    
    std::string name_; //!< The name of the variable
    std::vector<std::vector<size_t>> dim_lengths_pyramid_;
    GeoDomain init_domain_;
    DataLayout init_data_layout_;
    int num_decompression_lvls_{0};

    bool has_been_decompressed_{false};
};

template <typename T, size_t Dim>
void
AbstractPatchByteDecompressionVariable<T, Dim>::SetupDecompression()
{
    if (dim_lengths_pyramid_.back().size() != Dim)
    {
        cmc_err_msg("The dimensionality of the patch decompression variable does not coincide with the compressed data's dimensionality.");
    }
    /* The decompression always starts with a single element */
    data_.emplace_back();
}

template <typename T, size_t Dim>
void
AbstractPatchByteDecompressionVariable<T, Dim>::SetInitDataValue(const CompressionValue<T>& init_value)
{
    data_ = std::vector<CompressionValue<T>>{init_value};
}

template <typename T, size_t Dim>
std::vector<T>
AbstractPatchByteDecompressionVariable<T, Dim>::GetDecompressedData() const 
{
    if (not has_been_decompressed_)
    {
        cmc_warn_msg("The data has not yet been decompressed. It needs to be decompressed before it is actual usable (call Decompress()).");
    }

    /* Revert to the actual data type */
    std::vector<T> data;
    data.reserve(data_.size());
    
    /* Reinterpret the values */
    for (auto citer = data_.begin(); citer != data_.end(); ++citer)
    {
        data.push_back(citer->template ReinterpretDataAs<T>());
    }

    return data;
}

template <typename T, size_t Dim>
void
AbstractPatchByteDecompressionVariable<T, Dim>::Decompress()
{
    if constexpr (Dim == 1)
    {
        cmc_debug_msg("Lossless 1D Decompression");
        this->Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::tag1D);
    } else if constexpr (Dim == 2)
    {
        cmc_debug_msg("Lossless 2D Decompression");
        this->Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::tag2D);
    } else if constexpr (Dim == 3)
    {
        cmc_debug_msg("Lossless 3D Decompression");
        this->Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::tag3D);
    }
    else if constexpr (Dim == 4)
    {
        cmc_debug_msg("Lossless 4D Decompression");
        this->Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::tag4D);
    } else
    {
        cmc_err_msg("Unsupported variable's dimensionality (Dim = ", Dim, ").");
    }
}


/****** 4D Decompression ******/

template <typename T>
inline CompressionValue<T>
GetCoarseValue(const std::vector<CompressionValue<T>>& coarse_data, const std::vector<size_t>& coarse_lvl_dims, const size_t next_lvl_time, const size_t next_lvl_lev, const size_t next_lvl_lat, const size_t next_lvl_lon)
{
    cmc_assert(coarse_lvl_dims.size() == static_cast<size_t>(4));
    constexpr int kTimeID = 0;
    constexpr int kLevID = 1;
    constexpr int kLatID = 2;
    constexpr int kLonID = 3;

    const size_t coarse_lvl_time = next_lvl_time / kDimReductionFactor;
    const size_t coarse_lvl_lev = next_lvl_lev / kDimReductionFactor;
    const size_t coarse_lvl_lat = next_lvl_lat / kDimReductionFactor;
    const size_t coarse_lvl_lon = next_lvl_lon / kDimReductionFactor;

    cmc_assert(coarse_lvl_time < coarse_lvl_dims[kTimeID]);
    cmc_assert(coarse_lvl_lev < coarse_lvl_dims[kLevID]);
    cmc_assert(coarse_lvl_lat < coarse_lvl_dims[kLatID]);
    cmc_assert(coarse_lvl_lon < coarse_lvl_dims[kLonID]);
    cmc_assert(coarse_lvl_time * coarse_lvl_dims[kLevID] * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lev * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon < coarse_data.size());

    return coarse_data[coarse_lvl_time * coarse_lvl_dims[kLevID] * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lev * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon];
}

template <typename T, size_t Dim>
void
AbstractPatchByteDecompressionVariable<T, Dim>::Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::kTag4D)
{
    cmc_debug_msg("Decompression of variable ", this->name_, " starts...");

    constexpr int kDim = 4;
    constexpr int kTimeID = 0;
    constexpr int kLevID = 1;
    constexpr int kLatID = 2;
    constexpr int kLonID = 3;

    /* Decode the root level value */
    data_ = std::vector<CompressionValue<T>>{this->GetInitDecompressionValue()};

    cmc_assert(data_.size() == static_cast<size_t>(1));

    /* Perform the iterative decompression */
    for (int decompression_iter{0}; decompression_iter < num_decompression_lvls_ - 1; ++decompression_iter)
    {
        cmc_debug_msg("A decompression iteration is initialized (step: ", decompression_iter,").");
        this->InitializeDecompressionIteration();

        /* Get the dimensions of this and the next level */
        const std::vector<size_t>& curr_lvl_dims = dim_lengths_pyramid_[decompression_iter];
        const std::vector<size_t>& next_lvl_dims = dim_lengths_pyramid_[decompression_iter + 1];

        /* Allocate the data_new_*/
        data_new_.reserve(next_lvl_dims[kTimeID] * next_lvl_dims[kLevID] * next_lvl_dims[kLatID] * next_lvl_dims[kLonID]);

        /* Iterate over the finer level */
        for (size_t time = 0; time < next_lvl_dims[kTimeID]; ++time)
        {
            for (size_t lev = 0; lev < next_lvl_dims[kLevID]; ++lev)
            {
                for (size_t lat = 0; lat < next_lvl_dims[kLatID]; ++lat)
                {
                    for (size_t lon = 0; lon < next_lvl_dims[kLonID]; ++lon)
                    {
                        /* Get the corresponding coarse level value */
                        const CompressionValue<T> coarse_val = GetCoarseValue<T>(data_, curr_lvl_dims, time, lev, lat, lon);

                        /* Refine the value */
                        const CompressionValue<T> fine_val = this->RefineNextValue(coarse_val);

                        /* Store the refined value */
                        data_new_.push_back(fine_val);
                    }
                }
            }
        }

        /* Switch to the new refined data for the next iteration */
        data_ = std::move(data_new_);
        data_new_.clear();

        this->FinalizeDecompressionIteration();
        cmc_debug_msg("The decompression iteration is finished.");
    }

    has_been_decompressed_ = true;
    cmc_debug_msg("Decompression of variable ", this->name_, " is finished.");
}

/******************************/


/****** 3D Decompression ******/

template <typename T>
inline CompressionValue<T>
GetCoarseValue(const std::vector<CompressionValue<T>>& coarse_data, const std::vector<size_t>& coarse_lvl_dims, const size_t next_lvl_lev, const size_t next_lvl_lat, const size_t next_lvl_lon)
{
    cmc_assert(coarse_lvl_dims.size() == static_cast<size_t>(3));
    const size_t coarse_lvl_lev = next_lvl_lev / kDimReductionFactor;
    const size_t coarse_lvl_lat = next_lvl_lat / kDimReductionFactor;
    const size_t coarse_lvl_lon = next_lvl_lon / kDimReductionFactor;

    cmc_assert(coarse_lvl_lev < coarse_lvl_dims[0]);
    cmc_assert(coarse_lvl_lev * coarse_lvl_dims[1] * coarse_lvl_dims[2] + coarse_lvl_lat * coarse_lvl_dims[2] + coarse_lvl_lon < coarse_data.size());

    return coarse_data[coarse_lvl_lev * coarse_lvl_dims[1] * coarse_lvl_dims[2] + coarse_lvl_lat * coarse_lvl_dims[2] + coarse_lvl_lon];
}

template <typename T, size_t Dim>
void
AbstractPatchByteDecompressionVariable<T, Dim>::Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::kTag3D)
{
    cmc_debug_msg("Decompression of variable ", this->name_, " starts...");

    /* Decode the root level value */
    data_ = std::vector<CompressionValue<T>>{this->GetInitDecompressionValue()};

    cmc_assert(data_.size() == static_cast<size_t>(1));

    /* Perform the iterative decompression */
    for (int decompression_iter{0}; decompression_iter < num_decompression_lvls_ - 1; ++decompression_iter)
    {
        cmc_debug_msg("A decompression iteration is initialized (step: ", decompression_iter,").");
        this->InitializeDecompressionIteration();

        /* Get the dimensions of this and the next level */
        const std::vector<size_t>& curr_lvl_dims = dim_lengths_pyramid_[decompression_iter];
        const std::vector<size_t>& next_lvl_dims = dim_lengths_pyramid_[decompression_iter + 1];

        /* Allocate the data_new_*/
        data_new_.reserve(next_lvl_dims[0] * next_lvl_dims[1] * next_lvl_dims[2]);

        /* Iterate over the finer level */
        for (size_t lev = 0; lev < next_lvl_dims[0]; ++lev)
        {
            for (size_t lat = 0; lat < next_lvl_dims[1]; ++lat)
            {
                for (size_t lon = 0; lon < next_lvl_dims[2]; ++lon)
                {
                    /* Get the corresponding coarse level value */
                    const CompressionValue<T> coarse_val = GetCoarseValue<T>(data_, curr_lvl_dims, lev, lat, lon);

                    /* Refine the value */
                    const CompressionValue<T> fine_val = this->RefineNextValue(coarse_val);

                    /* Store the refined value */
                    data_new_.push_back(fine_val);
                }
            }
        }

        /* Switch to the new refined data for the next iteration */
        data_ = std::move(data_new_);
        data_new_.clear();

        this->FinalizeDecompressionIteration();
        cmc_debug_msg("The decompression iteration is finished.");
    }

    has_been_decompressed_ = true;
    cmc_debug_msg("Decompression of variable ", this->name_, " is finished.");
}

/******************************/


/****** 2D Decompression ******/

template <typename T>
inline CompressionValue<T>
GetCoarseValue(const std::vector<CompressionValue<T>>& coarse_data, const std::vector<size_t>& coarse_lvl_dims, const size_t next_lvl_lat, const size_t next_lvl_lon)
{
    cmc_assert(coarse_lvl_dims.size() == static_cast<size_t>(2));
    constexpr int kLatID = 0;
    constexpr int kLonID = 1;

    const size_t coarse_lvl_lat = next_lvl_lat / kDimReductionFactor;
    const size_t coarse_lvl_lon = next_lvl_lon / kDimReductionFactor;

    cmc_assert(coarse_lvl_lat < coarse_lvl_dims[kLatID]);
    cmc_assert(coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon < coarse_data.size());

    return coarse_data[coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon];
}

template <typename T, size_t Dim>
void
AbstractPatchByteDecompressionVariable<T, Dim>::Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::kTag2D)
{
    constexpr int kDim = 2;
    constexpr int kLatID = 0;
    constexpr int kLonID = 1;

    cmc_debug_msg("Decompression of variable ", this->name_, " starts...");

    /* Decode the root level value */
    data_ = std::vector<CompressionValue<T>>{this->GetInitDecompressionValue()};

    cmc_assert(data_.size() == static_cast<size_t>(1));

    /* Perform the iterative decompression */
    for (int decompression_iter{0}; decompression_iter < num_decompression_lvls_ - 1; ++decompression_iter)
    {
        cmc_debug_msg("A decompression iteration is initialized (step: ", decompression_iter,").");
        this->InitializeDecompressionIteration();

        /* Get the dimensions of this and the next level */
        const std::vector<size_t>& curr_lvl_dims = dim_lengths_pyramid_[decompression_iter];
        const std::vector<size_t>& next_lvl_dims = dim_lengths_pyramid_[decompression_iter + 1];

        /* Allocate the data_new_*/
        data_new_.reserve(next_lvl_dims[kLatID] * next_lvl_dims[kLonID]);

        /* Iterate over the finer level */
        for (size_t lat = 0; lat < next_lvl_dims[kLatID]; ++lat)
        {
            for (size_t lon = 0; lon < next_lvl_dims[kLonID]; ++lon)
            {
                /* Get the corresponding coarse level value */
                const CompressionValue<T> coarse_val = GetCoarseValue<T>(data_, curr_lvl_dims, lat, lon);

                /* Refine the value */
                const CompressionValue<T> fine_val = this->RefineNextValue(coarse_val);

                /* Store the refined value */
                data_new_.push_back(fine_val);
            }
        }

        /* Switch to the new refined data for the next iteration */
        data_ = std::move(data_new_);
        data_new_.clear();

        this->FinalizeDecompressionIteration();
        cmc_debug_msg("The decompression iteration is finished.");
    }

    has_been_decompressed_ = true;
    cmc_debug_msg("Decompression of variable ", this->name_, " is finished.");
}

/******************************/


/****** 1D Decompression ******/

template <typename T>
inline CompressionValue<T>
GetCoarseValue(const std::vector<CompressionValue<T>>& coarse_data, const std::vector<size_t>& coarse_lvl_dims, const size_t next_lvl_lon)
{
    cmc_assert(coarse_lvl_dims.size() == static_cast<size_t>(2));
    constexpr int kLonID = 0;

    const size_t coarse_lvl_lon = next_lvl_lon / kDimReductionFactor;

    cmc_assert(next_lvl_lon < coarse_lvl_dims[kLonID]);
    cmc_assert(coarse_lvl_lon < coarse_data.size());

    return coarse_data[coarse_lvl_lon];
}

template <typename T, size_t Dim>
void
AbstractPatchByteDecompressionVariable<T, Dim>::Decompress(AbstractPatchByteDecompressionVariable<T, Dim>::kTag1D)
{
    constexpr int kDim = 1;
    constexpr int kLonID = 0;

    cmc_debug_msg("Decompression of variable ", this->name_, " starts...");

    /* Decode the root level value */
    data_ = std::vector<CompressionValue<T>>{this->GetInitDecompressionValue()};

    cmc_assert(data_.size() == static_cast<size_t>(1));

    /* Perform the iterative decompression */
    for (int decompression_iter{0}; decompression_iter < num_decompression_lvls_ - 1; ++decompression_iter)
    {
        cmc_debug_msg("A decompression iteration is initialized (step: ", decompression_iter,").");
        this->InitializeDecompressionIteration();

        /* Get the dimensions of this and the next level */
        const std::vector<size_t>& curr_lvl_dims = dim_lengths_pyramid_[decompression_iter];
        const std::vector<size_t>& next_lvl_dims = dim_lengths_pyramid_[decompression_iter + 1];

        /* Allocate the data_new_*/
        data_new_.reserve(next_lvl_dims[kLonID]);

        /* Iterate over the finer level */
        for (size_t lon = 0; lon < next_lvl_dims[kLonID]; ++lon)
        {
            /* Get the corresponding coarse level value */
            const CompressionValue<T> coarse_val = GetCoarseValue<T>(data_, curr_lvl_dims, lon);

            /* Refine the value */
            const CompressionValue<T> fine_val = this->RefineNextValue(coarse_val);

            /* Store the refined value */
            data_new_.push_back(fine_val);
        }

        /* Switch to the new refined data for the next iteration */
        data_ = std::move(data_new_);
        data_new_.clear();

        this->FinalizeDecompressionIteration();
        cmc_debug_msg("The decompression iteration is finished.");
    }

    has_been_decompressed_ = true;
    cmc_debug_msg("Decompression of variable ", this->name_, " is finished.");
}

/******************************/

}


#endif /* !CMC_PATCH_BYTE_DECOMPRESSION_VARIABLE_HXX */
