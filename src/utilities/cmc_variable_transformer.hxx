#ifndef CMC_VARIABLE_TRANSFORMER_HXX
#define CMC_VARIABLE_TRANSFORMER_HXX

#include "utilities/cmc_variable_transformer_forward.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "t8code/cmc_t8_data_variables.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "t8code/cmc_t8_byte_variable.hxx"

namespace cmc
{
    
class TransformerInputToCompressionVariable
{
public:
    Var operator()(InputVar& transform_input_variable);

private:
    template<typename T> void MoveData(Variable<T>& destination, InputVariable<T>& source);
};


template<typename T>
void
TransformerInputToCompressionVariable::MoveData(Variable<T>& destination, InputVariable<T>& source)
{
    destination.name_ = std::move(source.name_);
    destination.global_domain_ = std::move(source.global_domain_);
    destination.attributes_.missing_value_ = source.missing_value_;
    destination.attributes_.add_offset_ = source.add_offset_;
    destination.attributes_.scale_factor_ = source.scale_factor_;
    destination.attributes_.global_context_information_ = source.global_context_information_;
    destination.attributes_.has_split_dimension_ = source.has_split_dimension_;
    destination.attributes_.initial_data_layout_ = source.initial_layout_;
    destination.attributes_.pre_compression_layout_ = source.pre_compression_layout_;
    destination.utilities_.interpolate_ = source.interpolation_;
    destination.utilities_.tracking_option_ = source.inaccuracy_tracking_;
    destination.data_ = std::move(source.data_);

    //Test:save the data
    //if (std::is_same_v<T, char>)
    #if 0
    const int num_elems = destination.data_.size();
    std::vector<short> converted_data;
    converted_data.reserve(destination.data_.size());
    const short missin_val = -32767;
    for (auto data_iter = destination.data_.begin(); data_iter != destination.data_.end(); ++data_iter)
    {
        if (!ApproxCompare(static_cast<short>(*data_iter), missin_val))
        {
            converted_data.push_back(static_cast<short>(*data_iter));
        }
    }
    std::string fname = "era5_initial_t_sfc_ordered_" + std::to_string(destination.attributes_.global_context_information_);
    FILE* file = fopen(fname.c_str(), "wb");
    fwrite(converted_data.data(), sizeof(T), converted_data.size(), file);
    fclose(file);
    #if 0
    const int num_elems = destination.data_.size();
    std::vector<float> converted_data;
    converted_data.reserve(destination.data_.size());
    const float missin_val = -32767.0;
    for (auto data_iter = destination.data_.begin(); data_iter != destination.data_.end(); ++data_iter)
    {
        if (!ApproxCompare(static_cast<float>(*data_iter), missin_val))
        {
            converted_data.push_back(static_cast<float>(*data_iter));
        }
    }
    std::string fname = "era5_initial_t_sfc_ordered_" + std::to_string(destination.attributes_.global_context_information_);
    FILE* file = fopen(fname.c_str(), "wb");
    fwrite(converted_data.data(), sizeof(T), converted_data.size(), file);
    fclose(file);
    #endif
    //converted_data.reserve(num_elems);
    //const float add_offset = 0.0970683052537568;
    //const float scale_factor = 2.96247040388686e-06;
    //for (int j = 0; j < num_elems; ++j)
    //{
    //    converted_data.push_back(scale_factor * static_cast<float>(local_data[j]) + add_offset);
    //}
    //FILE* file = fopen("era5_land_tp_lon_lat.bin", "wb");
    //fwrite(converted_data.data(), sizeof(float), num_elems, file);
    //fclose(file);
    //std::exit(1);
    #endif

    #if 0
    const float prevmval = -9.0E+33;
    const float mval = 65536.0;
    destination.attributes_.missing_value_ = mval;
    for (auto data_iter = destination.data_.begin(); data_iter != destination.data_.end(); ++data_iter)
    {
        if (static_cast<float>(*data_iter) <= -66000)
        {
            *data_iter = mval;
        }
    }

    #endif
    #if 0
    //Test for now, move all values into negative domain
    const float missin_val = -32767.0;
    //for (auto data_iter = destination.data_.begin(); data_iter != destination.data_.end(); ++data_iter)
    //{
    //    if (!ApproxCompare(static_cast<float>(*data_iter), missin_val, static_cast<float>(1.0)))
    //    {
    //        *data_iter *= -1.0;
    //    }
    //}
    const int num_elems = destination.data_.size();
    std::vector<float> converted_data;
    converted_data.reserve(destination.data_.size());
    //const float missin_val = -32767.0;
    for (auto data_iter = destination.data_.begin(); data_iter != destination.data_.end(); ++data_iter)
    {
        if (!ApproxCompare(static_cast<float>(*data_iter), missin_val, static_cast<float>(1.0)))
        {
            converted_data.push_back(static_cast<float>(*data_iter));
        }
    }
    std::string fname = "era5_initial_t_sfc_ordered_" + std::to_string(destination.attributes_.global_context_information_);
    FILE* file = fopen(fname.c_str(), "wb");
    fwrite(converted_data.data(), sizeof(T), converted_data.size(), file);
    fclose(file);
    #endif
}


class TransformerCompressionToOutputVariable
{
public:
    OutputVar operator()(Var& transform_compression_variable);

private:
    template<typename T> void MoveData(OutputVariable<T>& destination, Variable<T>& source);
};

template<typename T>
void
TransformerCompressionToOutputVariable::MoveData(OutputVariable<T>& destination, Variable<T>& source)
{
    destination.name_ = std::move(source.name_);
    destination.global_domain_ = std::move(source.global_domain_);
    destination.missing_value_ = source.attributes_.missing_value_;
    destination.layout_ = source.attributes_.initial_data_layout_;
    destination.data_ = std::move(source.data_);
}


class TransformerCompressionToByteVariable
{
public:
    ByteVar operator()(Var& transform_to_byte_variable);

private:
    template<typename T> void MoveData(ByteVariable<T>& destination, Variable<T>& source);
};

template<typename T>
void
TransformerCompressionToByteVariable::MoveData(ByteVariable<T>& destination, Variable<T>& source)
{
    destination.name_ = std::move(source.name_);
    destination.global_domain_ = std::move(source.global_domain_);
    destination.attributes_ = std::move(source.attributes_);
    destination.utilities_ = std::move(source.utilities_);
    destination.initial_data_ = std::move(source.data_);
    destination.mesh_ = std::move(source.mesh_);

    //TODO: remove below
    #if 0
    cmc_debug_msg("Remove this file IO");
    FILE* file = fopen("compr_initial_data.bin", "wb");
    fwrite(destination.initial_data_.data(), sizeof(float), destination.initial_data_.size(), file);
    fclose(file);
    #endif
}

}

#endif /* !CMC_VARIABLE_TRANSFORMER_HXX */
