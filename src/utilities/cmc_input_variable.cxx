#include "utilities/cmc_input_variable.hxx"

namespace cmc
{


static CmcInputVariable
SetupInputVar(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value, const DataLayout layout, const GeoDomain& domain)
{
    switch(type)
    {
        case CmcType::Double:
        {
            cmc_assert(std::holds_alternative<double>(missing_value));
            InputVariable<double> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<double>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Float:
        {
            cmc_assert(std::holds_alternative<float>(missing_value));
            InputVariable<float> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<float>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Int32_t:
        {
            cmc_assert(std::holds_alternative<int32_t>(missing_value));
            InputVariable<int32_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int32_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Uint32_t:
        {
            cmc_assert(std::holds_alternative<uint32_t>(missing_value));
            InputVariable<uint32_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint32_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Int64_t:
        {
            cmc_assert(std::holds_alternative<int64_t>(missing_value));
            InputVariable<int64_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int64_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Uint64_t:
        {
            cmc_assert(std::holds_alternative<uint64_t>(missing_value));
            InputVariable<uint64_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint64_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Int16_t:
        {
            cmc_assert(std::holds_alternative<int16_t>(missing_value));
            InputVariable<int16_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int16_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Uint16_t:
        {
            cmc_assert(std::holds_alternative<uint16_t>(missing_value));
            InputVariable<uint16_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint16_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Int8_t:
        {
            cmc_assert(std::holds_alternative<int8_t>(missing_value));
            InputVariable<int8_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int8_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        case CmcType::Uint8_t:
        {
            cmc_assert(std::holds_alternative<uint8_t>(missing_value));
            InputVariable<uint8_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint8_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return CmcInputVariable(std::move(var));
        }
        break;
        default:
            cmc_err_msg("The supplied data type (", type, ") is not supported");
            return CmcInputVariable(std::in_place_index<0>, "ErrorVariable", CMC_ERR, DataLayout::LayoutUndefined);
    }
}

InputVar::InputVar(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value, const DataLayout layout, const GeoDomain& domain)
{
    var_ = SetupInputVar(type, name, id, num_elements, missing_value, layout, domain);
}

struct ObtainCmcType
{
public:
    ObtainCmcType() = default;

    CmcType operator()([[maybe_unused]] const InputVariable<int8_t>& var) {
        return CmcType::Int8_t;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<char>& var) {
        return CmcType::Char; 
    }
    CmcType operator()([[maybe_unused]] const InputVariable<int16_t>& var) {
        return CmcType::Int16_t;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<int32_t>& var) {
        return CmcType::Int32_t;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<float>& var) {
        return CmcType::Float;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<double>& var) {
        return CmcType::Double;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<uint8_t>& var) {
        return CmcType::Uint8_t;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<uint16_t>& var) {
        return CmcType::Uint16_t;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<uint32_t>& var) {
        return CmcType::Uint32_t;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<int64_t>& var) {
        return CmcType::Int64_t;
    }
    CmcType operator()([[maybe_unused]] const InputVariable<uint64_t>& var) {
        return CmcType::Uint64_t;
    }
    template<typename T>
    CmcType operator()(const T& var) {
        return CmcType::TypeUndefined;
    }
};

InputVar MetaCopy(const InputVar& variable)
{
    return std::visit([&](auto&& var) -> InputVar{
        return InputVar(HollowCopy(var));
    }, variable.var_);

}

void
InputVar::SetUpFilledVariable(const size_t num_elements, const CmcUniversalType& fill_value)
{
    std::visit([&](auto&& var){
        var.SetUpFilledVariable(num_elements, fill_value);
    }, var_);
}

CmcType 
InputVar::GetType() const
{
    return std::visit(ObtainCmcType(), var_);
}


int
InputVar::GetID() const
{
    return std::visit([](auto&& var){
        return var.GetID();
    }, var_);
}

const std::string&
InputVar::GetName() const
{
    return std::visit([](auto&& var) -> const std::string& {
        return var.GetName();
    }, var_);   
}

void
InputVar::SetMissingValue(const CmcUniversalType& missing_value)
{
    std::visit([&](auto&& var){
        var.SetMissingValue(missing_value);
    }, var_);
}

void
InputVar::SetAddOffset(const CmcUniversalType& add_offset)
{
    std::visit([&](auto&& var){
        var.SetAddOffset(add_offset);
    }, var_);
}

void
InputVar::SetScaleFactor(const CmcUniversalType& scale_factor)
{
    std::visit([&](auto&& var){
        var.SetScaleFactor(scale_factor);
    }, var_);
}

void
InputVar::TransformCoordinatesToLinearIndices()
{
    std::visit([](auto&& var){
        var.TransformCoordinatesToMortonIndices();
    }, var_);
}

//TODO: The same funciton for the Variable message (for the parallel case)
void
InputVar::AssignDataAtLinearIndices(const InputVar& variable, const MortonIndex start_offset, const UpdateLinearIndices& update)
{
    std::visit([&](auto&& var){
        var.AssignDataAtLinearIndices(variable, start_offset, update);
    }, var_);
}

const GeoDomain&
InputVar::GetGlobalDomain() const
{
    return std::visit([](auto&& var) -> const GeoDomain& {
        return var.GetGlobalDomain();
    }, var_);
}

CmcUniversalType
InputVar::GetMissingValue() const
{
    return std::visit([](auto&& var){
        return CmcUniversalType{var.GetMissingValue()};
    }, var_);
}

CmcUniversalType 
InputVar::GetAddOffset() const
{
    return std::visit([](auto&& var){
        return CmcUniversalType{var.GetAddOffset()};
    }, var_);
}

CmcUniversalType
InputVar::GetScaleFactor() const
{
    return std::visit([](auto&& var){
        return CmcUniversalType{var.GetScaleFactor()};
    }, var_);
}

int
InputVar::GetGlobalContextInformation() const
{
    return std::visit([](auto&& var){
        return var.GetGlobalContextInformation();
    }, var_);
}

void
InputVar::SetGlobalContextInformation(const int global_context_information)
{
    std::visit([&](auto&& var){
        var.SetGlobalContextInformation(global_context_information);
    }, var_); 
}

DataLayout
InputVar::GetInitialDataLayout() const
{
    return std::visit([](auto&& var){
        return var.GetInitialDataLayout();
    }, var_);
}

const CmcInputVariable&
InputVar::GetInternalVariant() const
{
    return var_;
}


void
InputVar::ApplyScalingAndOffset()
{
    cmc_debug_msg("Vor apply offset ad scaling");
    std::visit([this](auto&& var) {
        const CmcDefaultDataType d_scale_factor = GetUniversalDataAs<CmcDefaultDataType>(var.GetScaleFactor());
        const CmcDefaultDataType d_add_offset = GetUniversalDataAs<CmcDefaultDataType>(var.GetAddOffset());
        
        if (!ApproxCompare(d_scale_factor, static_cast<CmcDefaultDataType>(1.0)) && !ApproxCompare(d_add_offset, static_cast<CmcDefaultDataType>(0.0)))
        {
            /* Scaling and offset need to be applied */
            ApplyAxpyScalingAndOffset(var);
        } else if (!ApproxCompare(d_scale_factor, static_cast<CmcDefaultDataType>(1.0)) && ApproxCompare(d_add_offset, static_cast<CmcDefaultDataType>(0.0)))
        {
            /* Only scaling needs to be applied */
            ApplyScaling(var);
        } else if (ApproxCompare(d_scale_factor, static_cast<CmcDefaultDataType>(1.0)) && !ApproxCompare(d_add_offset, static_cast<CmcDefaultDataType>(0.0)))
        {
            /* Only offsets need to be applied */
            ApplyOffset(var);
        }
    }, var_);
    //cmc_debug_msg("nach apply offset ad scaling");
};

struct SplitInputVariables
{
public:
    SplitInputVariables() = delete;
    SplitInputVariables(const Dimension dimension)
    : split_dimension{dimension}{};

    std::vector<InputVar> operator()(const InputVariable<int8_t>& var) {
        std::vector<InputVariable<int8_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<char>& var) {
        std::vector<InputVariable<char>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<int16_t>& var) {
        std::vector<InputVariable<int16_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<int32_t>& var) {
        std::vector<InputVariable<int32_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<float>& var) {
        std::vector<InputVariable<float>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<double>& var) {
        std::vector<InputVariable<double>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<uint8_t>& var) {
        std::vector<InputVariable<uint8_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<uint16_t>& var) {
        std::vector<InputVariable<uint16_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<uint32_t>& var) {
        std::vector<InputVariable<uint32_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<int64_t>& var) {
        std::vector<InputVariable<int64_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<InputVar> operator()(const InputVariable<uint64_t>& var) {
        std::vector<InputVariable<uint64_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<InputVar> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
private:
    const Dimension split_dimension;
};

std::vector<InputVar>
SplitIntoSubVariables(const InputVar& variable, const Dimension dimension)
{
    return std::visit(SplitInputVariables(dimension), variable.var_);
}

struct GatherSendData
{
public:
    GatherSendData() = delete;
    GatherSendData(const DataOffsets& mpi_offsets)
    : offsets{mpi_offsets}{};

    std::vector<VariableSendMessage> operator()(const InputVariable<int8_t>& var) {
        ReceiverMap<int8_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<int8_t>& var_message = send_data[sd_iter->first];
            VariableMessage<int8_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<char>& var) {
        ReceiverMap<char> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<char>& var_message = send_data[sd_iter->first];
            VariableMessage<char> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<int16_t>& var) {
        ReceiverMap<int16_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<int16_t>& var_message = send_data[sd_iter->first];
            VariableMessage<int16_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<int32_t>& var) {
        ReceiverMap<int32_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<int32_t>& var_message = send_data[sd_iter->first];
            VariableMessage<int32_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<float>& var) {
        ReceiverMap<float> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<float>& var_message = send_data[sd_iter->first];
            VariableMessage<float> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<double>& var) {
        ReceiverMap<double> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<double>& var_message = send_data[sd_iter->first];
            VariableMessage<double> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<uint8_t>& var) {
        ReceiverMap<uint8_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<uint8_t>& var_message = send_data[sd_iter->first];
            VariableMessage<uint8_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<uint16_t>& var) {
        ReceiverMap<uint16_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<uint16_t>& var_message = send_data[sd_iter->first];
            VariableMessage<uint16_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<uint32_t>& var) {
        ReceiverMap<uint32_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<uint32_t>& var_message = send_data[sd_iter->first];
            VariableMessage<uint32_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<int64_t>& var) {
        ReceiverMap<int64_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<int64_t>& var_message = send_data[sd_iter->first];
            VariableMessage<int64_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
    std::vector<VariableSendMessage> operator()(const InputVariable<uint64_t>& var) {
        ReceiverMap<uint64_t> send_data = GatherDataToBeDistributed(var, offsets);
        std::vector<VariableSendMessage> variables_send_data;
        variables_send_data.reserve(send_data.size());
        for (auto sd_iter = send_data.begin(); sd_iter != send_data.end(); ++sd_iter)
        {
            VariableMessage<uint64_t>& var_message = send_data[sd_iter->first];
            VariableMessage<uint64_t> send_var;
            DetachVarMessage(var_message, send_var);
            variables_send_data.emplace_back(std::move(send_var));
        }
        return variables_send_data;
    }
private:
    const DataOffsets& offsets;
};


std::vector<VariableSendMessage>
GatherDistributionData(const InputVar& variable, const DataOffsets& offsets)
{
    return std::visit(GatherSendData(offsets), variable.var_);
}




}
