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
    if (source.is_initial_data_kept_)
    {
        cmc_debug_msg("initial data is transferred");
        cmc_debug_msg("Size of initial data in transformation: ", source.initial_data_.size() );
        destination.uncompressed_mesh_ = std::move(source.initial_mesh_);
        destination.uncompressed_data_ = std::move(source.initial_data_);
        destination.are_uncompressed_states_stored_ = true;
    }
}

}

#endif /* !CMC_VARIABLE_TRANSFORMER_HXX */
