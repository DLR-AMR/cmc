#include "utilities/cmc_variable_transformer.hxx"

namespace cmc
{

Var
TransformerInputToCompressionVariable::operator()(InputVar& input_var)
{
    if (InputVariable<double>* vp = std::get_if<InputVariable<double>>(&input_var.var_))
    {
        Variable<double> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Double, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<float>* vp = std::get_if<InputVariable<float>>(&input_var.var_))
    {
        Variable<float> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Float, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<int32_t>* vp = std::get_if<InputVariable<int32_t>>(&input_var.var_))
    {
        Variable<int32_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Int32_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<int64_t>* vp = std::get_if<InputVariable<int64_t>>(&input_var.var_))
    {
        Variable<int64_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Int64_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<uint32_t>* vp = std::get_if<InputVariable<uint32_t>>(&input_var.var_))
    {
        Variable<uint32_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Uint32_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<uint64_t>* vp = std::get_if<InputVariable<uint64_t>>(&input_var.var_))
    {
        Variable<uint64_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Uint64_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<int16_t>* vp = std::get_if<InputVariable<int16_t>>(&input_var.var_))
    {
        Variable<int16_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Int16_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<uint16_t>* vp = std::get_if<InputVariable<uint16_t>>(&input_var.var_))
    {
        Variable<uint16_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Uint16_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<int8_t>* vp = std::get_if<InputVariable<int8_t>>(&input_var.var_))
    {
        Variable<int8_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Int8_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<uint8_t>* vp = std::get_if<InputVariable<uint8_t>>(&input_var.var_))
    {
        Variable<uint8_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Uint8_t, CmcVariable(std::move(transformed_variable)));
    } else if (InputVariable<char>* vp = std::get_if<InputVariable<char>>(&input_var.var_))
    {
        Variable<char> transformed_variable;
        MoveData(transformed_variable, *vp);
        return Var(input_var.GetID(), input_var.GetInternID(), CmcType::Char, CmcVariable(std::move(transformed_variable)));
    } else
    {
        cmc_err_msg("The input variable could not be transformed to a compression variable.");
        return Var(CMC_ERR, CmcType::TypeUndefined, CmcVariable(Variable<int32_t>()));
    }
}


OutputVar
TransformerCompressionToOutputVariable::operator()(Var& compression_var)
{
    if (Variable<double>* vp = std::get_if<Variable<double>>(&compression_var.var_))
    {
        OutputVariable<double> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<float>* vp = std::get_if<Variable<float>>(&compression_var.var_))
    {
        OutputVariable<float> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<int32_t>* vp = std::get_if<Variable<int32_t>>(&compression_var.var_))
    {
        OutputVariable<int32_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<int64_t>* vp = std::get_if<Variable<int64_t>>(&compression_var.var_))
    {
        OutputVariable<int64_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<uint32_t>* vp = std::get_if<Variable<uint32_t>>(&compression_var.var_))
    {
        OutputVariable<uint32_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<uint64_t>* vp = std::get_if<Variable<uint64_t>>(&compression_var.var_))
    {
        OutputVariable<uint64_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<int16_t>* vp = std::get_if<Variable<int16_t>>(&compression_var.var_))
    {
        OutputVariable<int16_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<uint16_t>* vp = std::get_if<Variable<uint16_t>>(&compression_var.var_))
    {
        OutputVariable<uint16_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<int8_t>* vp = std::get_if<Variable<int8_t>>(&compression_var.var_))
    {
        OutputVariable<int8_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<uint8_t>* vp = std::get_if<Variable<uint8_t>>(&compression_var.var_))
    {
        OutputVariable<uint8_t> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else if (Variable<char>* vp = std::get_if<Variable<char>>(&compression_var.var_))
    {
        OutputVariable<char> transformed_variable;
        MoveData(transformed_variable, *vp);
        return OutputVar(compression_var.GetID(), CmcOutputVariable(std::move(transformed_variable)));
    } else 
    {
        cmc_err_msg("The input variable could not be transformed to a compression variable.");
        return OutputVar(CMC_ERR, CmcOutputVariable(OutputVariable<int32_t>()));
    }
}


ByteVar
TransformerCompressionToByteVariable::operator()(Var& compression_var)
{
    if (Variable<double>* vp = std::get_if<Variable<double>>(&compression_var.var_))
    {
        ByteVariable<double> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Double, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<float>* vp = std::get_if<Variable<float>>(&compression_var.var_))
    {
        ByteVariable<float> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Float, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<int32_t>* vp = std::get_if<Variable<int32_t>>(&compression_var.var_))
    {
        ByteVariable<int32_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Int32_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<uint32_t>* vp = std::get_if<Variable<uint32_t>>(&compression_var.var_))
    {
        ByteVariable<uint32_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Uint32_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<int64_t>* vp = std::get_if<Variable<int64_t>>(&compression_var.var_))
    {
        ByteVariable<int64_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Int64_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<uint64_t>* vp = std::get_if<Variable<uint64_t>>(&compression_var.var_))
    {
        ByteVariable<uint64_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Uint64_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<int16_t>* vp = std::get_if<Variable<int16_t>>(&compression_var.var_))
    {
        ByteVariable<int16_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Int16_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<uint16_t>* vp = std::get_if<Variable<uint16_t>>(&compression_var.var_))
    {
        ByteVariable<uint16_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Uint16_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<int8_t>* vp = std::get_if<Variable<int8_t>>(&compression_var.var_))
    {
        ByteVariable<int8_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Int8_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<uint8_t>* vp = std::get_if<Variable<uint8_t>>(&compression_var.var_))
    {
        ByteVariable<uint8_t> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Uint8_t, CmcByteVariable(std::move(transformed_variable)));
    } else if (Variable<char>* vp = std::get_if<Variable<char>>(&compression_var.var_))
    {
        ByteVariable<char> transformed_variable(vp->GetDataAsCompressionValues());
        MoveData(transformed_variable, *vp);
        return ByteVar(compression_var.GetID(), CmcType::Char, CmcByteVariable(std::move(transformed_variable)));
    } else
    {
        cmc_err_msg("The input variable could not be transformed to a compression variable.");
        return ByteVar(CMC_ERR, CmcType::TypeUndefined, CmcByteVariable(ByteVariable<int32_t>()));
    }
}

}
