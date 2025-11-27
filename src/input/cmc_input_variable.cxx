#include "input/cmc_input_variable.hxx"

namespace cmc::input
{

static GeneralVariable
SetupVar(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value, const DataLayout layout, const GeoDomain& domain)
{
    switch(type)
    {
        case CmcType::Double:
        {
            cmc_assert(std::holds_alternative<double>(missing_value));
            Variable<double> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<double>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Float:
        {
            cmc_assert(std::holds_alternative<float>(missing_value));
            Variable<float> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<float>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int32_t:
        {
            cmc_assert(std::holds_alternative<int32_t>(missing_value));
            Variable<int32_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int32_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint32_t:
        {
            cmc_assert(std::holds_alternative<uint32_t>(missing_value));
            Variable<uint32_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint32_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int64_t:
        {
            cmc_assert(std::holds_alternative<int64_t>(missing_value));
            Variable<int64_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int64_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint64_t:
        {
            cmc_assert(std::holds_alternative<uint64_t>(missing_value));
            Variable<uint64_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint64_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int16_t:
        {
            cmc_assert(std::holds_alternative<int16_t>(missing_value));
            Variable<int16_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int16_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint16_t:
        {
            cmc_assert(std::holds_alternative<uint16_t>(missing_value));
            Variable<uint16_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint16_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int8_t:
        {
            cmc_assert(std::holds_alternative<int8_t>(missing_value));
            Variable<int8_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<int8_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint8_t:
        {
            cmc_assert(std::holds_alternative<uint8_t>(missing_value));
            Variable<uint8_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, std::get<uint8_t>(missing_value));
            var.SetGlobalDomain(domain);
            var.SetMissingValue(missing_value);
            return GeneralVariable(std::move(var));
        }
        break;
        default:
            cmc_err_msg("The supplied data type (", type, ") is not supported");
            return GeneralVariable(std::in_place_index<0>, "ErrorVariable", CMC_ERR, DataLayout::LayoutUndefined);
    }
}


static GeneralVariable
SetupVar(const CmcType type, const std::string& name, const int id, const size_t num_elements, const DataLayout layout, const GeoDomain& domain)
{
    switch(type)
    {
        case CmcType::Double:
        {
            Variable<double> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, double{0.0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Float:
        {
            Variable<float> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, float{0.0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int32_t:
        {
            Variable<int32_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, int32_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint32_t:
        {
            Variable<uint32_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, uint32_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int64_t:
        {
            Variable<int64_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, int64_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint64_t:
        {
            Variable<uint64_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, uint64_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int16_t:
        {
            Variable<int16_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, int16_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint16_t:
        {
            Variable<uint16_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, uint16_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Int8_t:
        {
            Variable<int8_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, int8_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        case CmcType::Uint8_t:
        {
            Variable<uint8_t> var(name, id, layout);
            var.SetUpFilledVariable(num_elements, uint8_t{0});
            var.SetGlobalDomain(domain);
            return GeneralVariable(std::move(var));
        }
        break;
        default:
            cmc_err_msg("The supplied data type (", type, ") is not supported");
            return GeneralVariable(std::in_place_index<0>, "ErrorVariable", CMC_ERR, DataLayout::LayoutUndefined);
    }
}

Var::Var(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value, const DataLayout layout, const GeoDomain& domain)
{
    var_ = SetupVar(type, name, id, num_elements, missing_value, layout, domain);
}

Var::Var(const CmcType type, const std::string& name, const int id, const size_t num_elements, const DataLayout layout, const GeoDomain& domain)
{
    var_ = SetupVar(type, name, id, num_elements, layout, domain);
}

struct ObtainCmcType
{
public:
    ObtainCmcType() = default;

    CmcType operator()([[maybe_unused]] const Variable<int8_t>& var) {
        return CmcType::Int8_t;
    }
    CmcType operator()([[maybe_unused]] const Variable<char>& var) {
        return CmcType::Char; 
    }
    CmcType operator()([[maybe_unused]] const Variable<int16_t>& var) {
        return CmcType::Int16_t;
    }
    CmcType operator()([[maybe_unused]] const Variable<int32_t>& var) {
        return CmcType::Int32_t;
    }
    CmcType operator()([[maybe_unused]] const Variable<float>& var) {
        return CmcType::Float;
    }
    CmcType operator()([[maybe_unused]] const Variable<double>& var) {
        return CmcType::Double;
    }
    CmcType operator()([[maybe_unused]] const Variable<uint8_t>& var) {
        return CmcType::Uint8_t;
    }
    CmcType operator()([[maybe_unused]] const Variable<uint16_t>& var) {
        return CmcType::Uint16_t;
    }
    CmcType operator()([[maybe_unused]] const Variable<uint32_t>& var) {
        return CmcType::Uint32_t;
    }
    CmcType operator()([[maybe_unused]] const Variable<int64_t>& var) {
        return CmcType::Int64_t;
    }
    CmcType operator()([[maybe_unused]] const Variable<uint64_t>& var) {
        return CmcType::Uint64_t;
    }
    template<typename T>
    CmcType operator()(const T& var) {
        return CmcType::TypeUndefined;
    }
};

Var MetaCopy(const Var& variable)
{
    return std::visit([&](auto&& var) -> Var{
        return Var(HollowCopy(var));
    }, variable.var_);

}

void
Var::SetUpFilledVariable(const size_t num_elements, const CmcUniversalType& fill_value)
{
    std::visit([&](auto&& var){
        var.SetUpFilledVariable(num_elements, fill_value);
    }, var_);
}

CmcType 
Var::GetType() const
{
    return std::visit(ObtainCmcType(), var_);
}


int
Var::GetID() const
{
    return std::visit([](auto&& var){
        return var.GetID();
    }, var_);
}

const std::string&
Var::GetName() const
{
    return std::visit([](auto&& var) -> const std::string& {
        return var.GetName();
    }, var_);   
}

void
Var::SetMissingValue(const CmcUniversalType& missing_value)
{
    std::visit([&](auto&& var){
        var.SetMissingValue(missing_value);
    }, var_);
}

void
Var::SetAddOffset(const CmcUniversalType& add_offset)
{
    std::visit([&](auto&& var){
        var.SetAddOffset(add_offset);
    }, var_);
}

void
Var::SetScaleFactor(const CmcUniversalType& scale_factor)
{
    std::visit([&](auto&& var){
        var.SetScaleFactor(scale_factor);
    }, var_);
}

int
Var::GetInternID() const
{
    return std::visit([](auto&& var) -> int {
        return var.GetInternID();
    }, var_);
}

GeneralVariable&
Var::GetInternalVariant([[maybe_unused]] const AccessKey& key)
{
    return var_;
}

void
Var::SetInternID(const int id)
{
    std::visit([&](auto&& var){
        var.SetInternID(id);
    }, var_);
}

DataFormat
Var::GetActiveDataFormat() const
{
    return std::visit([](auto&& var) -> DataFormat {
        return var.GetActiveDataFormat();
    }, var_);
}
void
Var::TransformCoordinatesToLinearIndices()
{
    std::visit([](auto&& var){
        var.TransformCoordinatesToMortonIndices();
    }, var_);
}

//TODO: The same funciton for the Variable message (for the parallel case)
void
Var::AssignDataAtLinearIndices(const Var& variable, const UpdateLinearIndices& update)
{
    std::visit([&](auto&& var){
        var.AssignDataAtLinearIndices(variable, update);
    }, var_);
}

const GeoDomain&
Var::GetGlobalDomain() const
{
    return std::visit([](auto&& var) -> const GeoDomain& {
        return var.GetGlobalDomain();
    }, var_);
}

CmcUniversalType
Var::GetMissingValue() const
{
    return std::visit([](auto&& var){
        return CmcUniversalType{var.GetMissingValue()};
    }, var_);
}

CmcUniversalType 
Var::GetAddOffset() const
{
    return std::visit([](auto&& var){
        return CmcUniversalType{var.GetAddOffset()};
    }, var_);
}

CmcUniversalType
Var::GetScaleFactor() const
{
    return std::visit([](auto&& var){
        return CmcUniversalType{var.GetScaleFactor()};
    }, var_);
}

int
Var::GetGlobalContextInformation() const
{
    return std::visit([](auto&& var){
        return var.GetGlobalContextInformation();
    }, var_);
}

void
Var::SetGlobalContextInformation(const int global_context_information)
{
    std::visit([&](auto&& var){
        var.SetGlobalContextInformation(global_context_information);
    }, var_); 
}

DataLayout
Var::GetInitialDataLayout() const
{
    return std::visit([](auto&& var){
        return var.GetInitialDataLayout();
    }, var_);
}

const GeneralVariable&
Var::GetInternalVariant() const
{
    return var_;
}

void
Var::ApplyScalingAndOffset()
{
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
}

class SplitVariables
{
public:
    SplitVariables() = delete;
    SplitVariables(const Dimension dimension)
    : split_dimension{dimension}{};

    std::vector<Var> operator()(const Variable<int8_t>& var) {
        std::vector<Variable<int8_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<char>& var) {
        std::vector<Variable<char>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<int16_t>& var) {
        std::vector<Variable<int16_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<int32_t>& var) {
        std::vector<Variable<int32_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<float>& var) {
        std::vector<Variable<float>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<double>& var) {
        std::vector<Variable<double>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<uint8_t>& var) {
        std::vector<Variable<uint8_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<uint16_t>& var) {
        std::vector<Variable<uint16_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<uint32_t>& var) {
        std::vector<Variable<uint32_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<int64_t>& var) {
        std::vector<Variable<int64_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
    std::vector<Var> operator()(const Variable<uint64_t>& var) {
        std::vector<Variable<uint64_t>> split_variables = ExtractSubVariables(var, split_dimension);
        std::vector<Var> variables(std::make_move_iterator(split_variables.begin()), std::make_move_iterator(split_variables.end()));
        return variables;
    }
private:
    const Dimension split_dimension;
};

std::vector<Var>
SplitIntoSubVariables(const Var& variable, const Dimension dimension)
{
    return std::visit(SplitVariables(dimension), variable.var_);
}

CmcType
GetDataTypeFromVariableViaID(const std::vector<Var>& input_variables, const int variable_id)
{
    /* Find the variable with the corresponding id */
    auto var_iter = std::find_if(input_variables.begin(), input_variables.end(), [&](auto& var){
        return (var.GetID() == variable_id);
    });

    /* Check if the variable has been found and return its corresponding CmcType */
    if (var_iter != input_variables.end())
    {
        return var_iter->GetType();
    } else
    {
        /* In case there is no variable with the corresponding ID */
        return CmcType::TypeUndefined;
    }
}

CmcType
GetDataTypeFromVariableViaInternID(const std::vector<Var>& input_variables, const int intern_id)
{
    /* Find the variable with the corresponding id */
    auto var_iter = std::find_if(input_variables.begin(), input_variables.end(), [&](auto& var){
        return (var.GetInternID() == intern_id);
    });

    /* Check if the variable has been found and return its corresponding CmcType */
    if (var_iter != input_variables.end())
    {
        return var_iter->GetType();
    } else
    {
        /* In case there is no variable with the corresponding ID */
        return CmcType::TypeUndefined;
    }
}


#ifdef CMC_ENABLE_MPI

void
Var::SetMPIComm(const MPI_Comm comm)
{
    std::visit([&](auto&& var){
        var.SetMPIComm(comm);
    }, var_);
}

MPI_Comm
Var::GetMPIComm() const
{
    return std::visit([](auto&& var) -> MPI_Comm {
        return var.GetMPIComm();
    }, var_);
}

void
Var::AssignDataAtLinearIndices(const VariableRecvMessage& message, const UpdateLinearIndices& update)
{
    std::visit([&](auto&& var){
        var.AssignDataAtLinearIndices(message, update);
    }, var_);
}


void
Var::GatherDistributionData(const DataOffsets& offsets, std::vector<VariableSendMessage>& messages)
{
    std::visit(overloaded {
            [&offsets, &messages](Variable<int8_t>& var)
            {
                ReceiverMap<int8_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<char>& var)
            {
                ReceiverMap<char> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<int16_t>& var)
            {
                ReceiverMap<int16_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<int32_t>& var)
            {
                ReceiverMap<int32_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<float>& var)
            {
                ReceiverMap<float> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<double>& var)
            {
                ReceiverMap<double> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<uint8_t>& var)
            {
                ReceiverMap<uint8_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<uint16_t>& var)
            {
                ReceiverMap<uint16_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<uint32_t>& var)
            {
                ReceiverMap<uint32_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<int64_t>& var)
            {
                ReceiverMap<int64_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [&offsets, &messages](Variable<uint64_t>& var)
            {
                ReceiverMap<uint64_t> send_data = var.GatherDataToBeDistributed(offsets);
                AppendSendData(messages, std::move(send_data));
            },
            [](auto& var){cmc_err_msg("Type Error in GatherDistributionData");}
        }, var_);
}

#endif

}

