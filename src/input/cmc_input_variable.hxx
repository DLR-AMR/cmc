#ifndef INPUT_CMC_INPUT_VARIABLE_HXX
#define INPUT_CMC_INPUT_VARIABLE_HXX

#include "input/cmc_forward_decl_input_variable.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_dimension_interval.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_cart_coordinate.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_geo_utilities.hxx"
#include "utilities/cmc_morton.hxx"
#include "input/cmc_binary_reader_forward.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
//#include "mpi/cmc_mpi_data.hxx"
#endif

//#include "t8code/cmc_t8_mpi.hxx"


#include <vector>

namespace cmc::input
{

constexpr int kNoGlobalContext = -1;
constexpr int kNoInternalIDSet = -1;

template<class T>
class Variable
{
public:
    Variable() = default;
    Variable(const std::string& name, const int id, const DataLayout layout)
    : name_{name}, id_{id}, initial_layout_{layout}{}
    Variable(std::string&& name, const int id, const DataLayout layout)
    : name_{std::move(name)}, id_{id}, initial_layout_{layout}{}
    Variable(const DataLayout layout, const DataFormat format)
    : initial_layout_{layout}, active_format_{format}{}
    ~Variable() = default;

    Variable(const Variable& other) = default;
    Variable& operator=(const Variable& other) = default;
    Variable(Variable&& other) = default;
    Variable& operator=(Variable&& other) = default;

    template<typename U> auto operator[](U index) const -> std::enable_if_t<std::is_integral_v<U>, T> {return data_[index];}
    template<typename U> auto operator[](U index) -> std::enable_if_t<std::is_integral_v<U>, T&> {return data_[index];}

    int GetID() const;
    CmcType GetType() const;
    const std::string& GetName() const;

    void SetUpFilledVariable(const size_t num_elements, const T fill_value);
    void SetUpFilledVariable(const size_t num_elements, const CmcUniversalType& fill_value);

    void PushBack(const T value, const LinearIndex index);

    void PushBack(const T value, const CartesianCoordinate& coordinate);
    void PushBack(const T value, CartesianCoordinate&& coordinate);
    void PushBack(const std::vector<T>& values, const Hyperslab& hyperslab);
    void PushBack(const std::vector<T>& values, Hyperslab&& hyperslab);
    void SetDataAndCoordinates(std::vector<T>&& values, std::vector<Hyperslab>&& hyperslabs);
    size_t GetNumberCoordinates() const;

    LinearIndex GetLinearIndexCoordinate(const size_t position) const;
    CartesianCoordinate GetCartesianCoordinate(const size_t position) const;
    Hyperslab GetHyperslab(const size_t position) const;

    DataLayout GetInitialDataLayout() const;
    DataLayout GetPreCompressionDataLayout() const;
    GeoDomain& GetGlobalDomain();
    const GeoDomain& GetGlobalDomain() const;
    void SetGlobalDomain(const GeoDomain& domain);
    void SetGlobalDomain(GeoDomain&& domain);

    int GetInternID() const;
    void SetInternID(const int id);
    DataFormat GetActiveDataFormat() const;
    CmcUniversalType GetAddOffset() const;
    CmcUniversalType GetScaleFactor() const;
    T GetMissingValue() const;
    void SetMissingValue(const CmcUniversalType& missing_value);
    void SetMissingValue(const T missing_value);
    void SetScaleFactor(const CmcUniversalType& scale_factor);
    void SetScaleFactor(const T scale_factor);
    void SetAddOffset(const CmcUniversalType& add_offset);
    void SetAddOffset(const T add_offset);
#ifdef CMC_ENABLE_MPI
    void SetMPIComm(const MPI_Comm comm) {comm_ = comm;}
    MPI_Comm GetMPIComm() const {return comm_;}
    ReceiverMap<T> GatherDataToBeDistributed(const DataOffsets& offsets);
    void AssignDataAtLinearIndices(const VariableRecvMessage& message, const UpdateLinearIndices& update);
#endif
    int GetGlobalContextInformation() const;
    void SetGlobalContextInformation(const int global_context_information);

    DataFormat GetActiveFormat() const;
    void TransformCoordinatesToMortonIndices();

    void Clear();
    bool IsValid() const;

    //TODO: Make domain index
    size_t GetIndexWithinDimension(const Dimension dimension, const size_t coordinate_position) const;

    std::vector<T>&& DetachData();
    void ClearData();

    std::vector<LinearIndex>&& DetachMortonIndices();
    void ClearMortonIndices();

    GeoDomain&& DetachGlobalDomain();

    void AssignDataAtLinearIndices(const Var& var, const UpdateLinearIndices& update);
    void AssignDataAtLinearIndices(const Variable& source_var, const UpdateLinearIndices& update);

    struct ExtractedVar
    {
        std::vector<T> data;
        std::vector<MortonIndex> morton_indices;
    };

    std::vector<ExtractedVar>
    ExtractAllSubVariablesByDimensionFromHyperslabs(const Dimension dimension) const;

    std::vector<ExtractedVar>
    ExtractAllSubVariablesByDimension(const Dimension dimension) const;

    template<typename U> Variable<U> ObtainMetaCopy() const;
    
    friend std::vector<Variable> ExtractSubVariables <> (const Variable& variable, const Dimension split_dimension);
    friend Variable HollowCopy <> (const Variable& variable);

    template <class U> friend class Variable;
    friend class TransformerInputToCompressionVariable;

    std::vector<T> GetDataFromHyperslab(const Hyperslab& hyperslab) const;

    /* Those functions are only accessible for certain functions */
    const std::vector<T>& GetDataForReading() const {return data_;}
    const std::vector<T>& GetData([[maybe_unused]] const AccessKey& key) const {return data_;}
    void SetData([[maybe_unused]] const AccessKey& key, std::vector<T>&& data) {data_ = std::move(data);}
    void MoveDataInto(std::vector<T>& vec_to_move_data_into);

private:
    std::string name_;
    int id_;
    CmcType type_{ConvertToCmcType<T>()};

    T missing_value_{std::numeric_limits<T>::lowest()};
    CmcUniversalType add_offset_{static_cast<double>(0.0)};
    CmcUniversalType scale_factor_{static_cast<double>(1.0)};

    std::vector<T> data_;

    std::vector<CartesianCoordinate> cartesian_coordinates_;
    std::vector<LinearIndex> linear_indices_;
    std::vector<Hyperslab> hyperslabs_;

    DataLayout initial_layout_{DataLayout::LayoutUndefined};
    DataFormat active_format_{DataFormat::FormatUndefined};

    GeoDomain global_domain_;

    int intern_id_{kNoInternalIDSet};
    int global_context_information_{kNoGlobalContext};
    Dimension has_split_dimension_{Dimension::DimensionUndefined};
    bool _has_been_invalidated_by_moving_{false};
    DataLayout pre_compression_layout_{DataLayout::LayoutUndefined};

#ifdef CMC_ENABLE_MPI
    MPI_Comm comm_{MPI_COMM_NULL};
    void CreateInternIDForRedistribution();
#endif
};


class Var
{
public:
    Var() = default;
    Var(const CmcType type, const std::string& name, const int id, const size_t num_elements, const DataLayout layout, const GeoDomain& domain);
    Var(const CmcType type, const std::string& name, const int id, const size_t num_elements, const CmcUniversalType missing_value, const DataLayout layout, const GeoDomain& domain);

    template<typename T> Var(const Variable<T>& var);

    template<typename T> Var(Variable<T>&& var);

    ~Var() = default;

    Var(const Var& other) = default;
    Var& operator=(const Var& other) = default;
    Var(Var&& other) = default;
    Var& operator=(Var&& other) = default;

    int GetID() const;

    const std::string& GetName() const;

    CmcType GetType() const;

    void TransformCoordinatesToLinearIndices();
    
    void AssignDataAtLinearIndices(const Var& variable, const UpdateLinearIndices& update);

    const GeoDomain& GetGlobalDomain() const;

    CmcUniversalType GetMissingValue() const;

    DataLayout GetInitialDataLayout() const;

    const GeneralVariable& GetInternalVariant() const;
    GeneralVariable& GetInternalVariant([[maybe_unused]] const AccessKey& key);

    void ApplyScalingAndOffset();

    CmcUniversalType GetAddOffset() const;
    CmcUniversalType GetScaleFactor() const;

    void SetMissingValue(const CmcUniversalType& missing_value);
    void SetAddOffset(const CmcUniversalType& add_offset);
    void SetScaleFactor(const CmcUniversalType& scale_factor);
    int GetGlobalContextInformation() const;
    void SetGlobalContextInformation(const int global_context_information);
    int GetInternID() const;
    void SetInternID(const int id);
#ifdef CMC_ENABLE_MPI
    void SetMPIComm(const MPI_Comm comm);
    MPI_Comm GetMPIComm() const;
    void AssignDataAtLinearIndices(const VariableRecvMessage& message, const UpdateLinearIndices& update);
    void GatherDistributionData(const DataOffsets& offsets, std::vector<VariableSendMessage>& messages);
#endif
    DataFormat GetActiveDataFormat() const;
    void SetUpFilledVariable(const size_t num_elements, const CmcUniversalType& fill_value);

    template<typename T> bool holds_alternative() const;

    friend std::vector<Var> SplitIntoSubVariables(const Var& variable, const Dimension dimension);

    friend class TransformerInputToCompressionVariable;

    template <typename T> friend class Variable;

    friend Var MetaCopy(const Var& variable);

    template<typename T> void SetData([[maybe_unused]] const AccessKey& key, std::vector<T>&& data) {std::visit([&](auto& var){var.SetData(key, std::move(data));}, var_);}
    friend class AccessKey;
private:
    template<class T> void ApplyAxpyScalingAndOffset(const Variable<T>&);
    template<class T> void ApplyScaling(const Variable<T>&);
    template<class T> void ApplyOffset(const Variable<T>&);

    GeneralVariable var_;
};


struct IndexReduction
{
    IndexReduction(const MortonIndex index, MortonIndex skipped_elements)
    : start_index{index}, indices_to_subtract{skipped_elements}{}

    MortonIndex start_index{0};
    MortonIndex indices_to_subtract{0};
};

class UpdateLinearIndices
{
public:
    UpdateLinearIndices() = delete;
    UpdateLinearIndices(std::vector<IndexReduction>&& correction)
    : correction_{correction}{}
    UpdateLinearIndices(const std::vector<IndexReduction>& correction)
    : correction_{correction}{}

    template<typename T>
    auto operator()(T index_to_be_updated) const
     -> std::enable_if_t<std::is_same_v<std::decay_t<T>, MortonIndex>, MortonIndex>;

private:
    std::vector<IndexReduction> correction_;
};

class AccessKey
{
private:
    AccessKey(){}
    AccessKey([[maybe_unused]] const AccessKey& other){}

    template<class T> friend void Var::ApplyAxpyScalingAndOffset(const Variable<T>&);
    template<class T> friend void Var::ApplyScaling(const Variable<T>&);
    template<class T> friend void Var::ApplyOffset(const Variable<T>&);

    friend class input::binary::Reader;
};

CmcType
GetDataTypeFromVariableViaID(const std::vector<Var>& input_variables, const int variable_id);

CmcType
GetDataTypeFromVariableViaInternID(const std::vector<Var>& input_variables, const int intern_id);


/** VARIABLE<T> MEMBER FUNCITONS **/

template<typename T>
template<typename U>
Variable<U>
Variable<T>::ObtainMetaCopy() const
{
    Variable<U> hollow_variable(name_, id_, initial_layout_);
    hollow_variable.missing_value_ = static_cast<U>(missing_value_);
    hollow_variable.add_offset_ = add_offset_;
    hollow_variable.scale_factor_ = scale_factor_;
    hollow_variable.active_format_ = active_format_;
    hollow_variable.global_domain_ = global_domain_;
    hollow_variable.global_context_information_ = global_context_information_;
    hollow_variable.has_split_dimension_ = has_split_dimension_;
    hollow_variable._has_been_invalidated_by_moving_ = _has_been_invalidated_by_moving_;
    hollow_variable.pre_compression_layout_ = pre_compression_layout_;
#ifdef CMC_ENABLE_MPI
    hollow_variable.comm_ = comm_;
#endif
    hollow_variable.intern_id_ = intern_id_;
    hollow_variable.cartesian_coordinates_ = cartesian_coordinates_;
    hollow_variable.linear_indices_ = linear_indices_;
    hollow_variable.hyperslabs_ = hyperslabs_;
    hollow_variable.type_ = ConvertToCmcType<U>();
    return hollow_variable;
}

template<class T>
void Variable<T>::SetUpFilledVariable(const size_t num_elements, const T fill_value)
{
    std::vector<T> filled_vector(num_elements, fill_value);
    std::swap(filled_vector, data_); 
}

template<class T>
void Variable<T>::SetUpFilledVariable(const size_t num_elements, const CmcUniversalType& fill_value)
{
    cmc_assert(std::holds_alternative<T>(fill_value));
    std::vector<T> filled_vector(num_elements, std::get<T>(fill_value));
    std::swap(filled_vector, data_); 
}

template<class T>
void Variable<T>::PushBack(const T value, const LinearIndex index)
{
    data_.push_back(value);
    linear_indices_.push_back(index);
}

template<class T>
void Variable<T>::PushBack(const T value, const CartesianCoordinate& coordinate)
{
    data_.push_back(value);
    cartesian_coordinates_.push_back(coordinate);
}

template<class T>
void Variable<T>::PushBack(const T value, CartesianCoordinate&& coordinate)
{
    data_.push_back(value);
    cartesian_coordinates_.push_back(std::move(coordinate));
}

template<class T>
void Variable<T>::PushBack(const std::vector<T>& values, const Hyperslab& hyperslab)
{
    std::copy(values.begin(), values.end(), std::back_inserter(data_));
    hyperslabs_.push_back(hyperslab);
}

template<class T>
void Variable<T>::PushBack(const std::vector<T>& values, Hyperslab&& hyperslab)
{
    std::copy(values.begin(), values.end(), std::back_inserter(data_));
    hyperslabs_.push_back(std::move(hyperslab));
}

template<class T>
void Variable<T>::SetDataAndCoordinates(std::vector<T>&& values, std::vector<Hyperslab>&& hyperslabs)
{
    data_ = std::move(values);
    hyperslabs_ = std::move(hyperslabs);
    active_format_ = DataFormat::HyperslabFormat;
}

template<class T>
LinearIndex Variable<T>::GetLinearIndexCoordinate(const size_t position) const
{
    return linear_indices_[position];
}

template<class T>
CartesianCoordinate Variable<T>::GetCartesianCoordinate(const size_t position) const
{
    return cartesian_coordinates_[position];
}

template<class T>
Hyperslab Variable<T>::GetHyperslab(const size_t position) const
{
    return hyperslabs_[position];
}

template<class T>
DataLayout Variable<T>::GetInitialDataLayout() const
{
    return initial_layout_;
}

template<class T>
DataLayout Variable<T>::GetPreCompressionDataLayout() const
{
    return pre_compression_layout_;
}

template<class T>
GeoDomain& Variable<T>::GetGlobalDomain()
{
    return global_domain_;
}

template<class T>
const GeoDomain& 
Variable<T>::GetGlobalDomain() const
{
    return global_domain_;
}

template<class T>
void Variable<T>::SetGlobalDomain(const GeoDomain& domain)
{
    global_domain_ = domain;
}

template<class T>
void Variable<T>::SetGlobalDomain(GeoDomain&& domain)
{
    global_domain_ = std::move(domain);
}

template<class T>
void Variable<T>::MoveDataInto(std::vector<T>& vec_to_move_data_into)
{
    vec_to_move_data_into = std::move(data_);
    data_ = std::vector<T>();
}

template<class T>
void Variable<T>::Clear()
{
    data_.clear();
    cartesian_coordinates_.clear();
    linear_indices_.clear();
    hyperslabs_.clear();
    active_format_ = DataFormat::FormatUndefined;
}

template<class T>
void
Variable<T>::SetScaleFactor(const CmcUniversalType& scale_factor)
{
    scale_factor_ = scale_factor;
}

template<class T>
void
Variable<T>::SetScaleFactor(const T scale_factor)
{
    scale_factor_ = scale_factor;
}

template<class T>
void
Variable<T>::SetAddOffset(const CmcUniversalType& add_offset)
{
    add_offset_ = add_offset;
}

template<class T>
void
Variable<T>::SetAddOffset(const T add_offset)
{
    add_offset_ = add_offset;
}

template<class T>
DataFormat
Variable<T>::GetActiveFormat() const
{
    return active_format_;
}

template<class T>
size_t
Variable<T>::GetNumberCoordinates() const
{
    cmc_debug_msg("In Number coords: Active format is: ", active_format_);
    switch (active_format_)
    {
        case DataFormat::LinearFormat:
            return linear_indices_.size();
        break;
        case DataFormat::CartesianFormat:
            return cartesian_coordinates_.size();
        break;
        case DataFormat::HyperslabFormat:
        {
            HyperslabIndex count = 1;
            for (auto hs_iter = hyperslabs_.begin(); hs_iter != hyperslabs_.end(); ++hs_iter)
            {
                count *= hs_iter->GetNumberCoordinates();
            }
            return count;
        }
        break;
        default:
        cmc_debug_msg("Is default case here");
            return 0;
    }
}

template<class T>
void
Variable<T>::TransformCoordinatesToMortonIndices()
{
    cmc_debug_msg("In TransformCoordinatesToMortonIndices: active_format: ", GetActiveDataFormat());

    switch (GetActiveDataFormat())
    {
        case DataFormat::LinearFormat:
            return;
        break;
        case DataFormat::CartesianFormat:
            //TODO: maybe this is suited for a parallel search indentifying elements corresponding to the points or Morton
            // indices calculation for integer coordinates
            cmc_err_msg("Currently, it is not possible to transform coordinates from the Cartesian format to Morton indices.");
            //TransformCartesianCoordinatesToMortonIndex(const std::vector<CartesianCoordinate>& coordiantes, const GeoDomain& global_domain);
        break;
        case DataFormat::HyperslabFormat:
            linear_indices_ = TransformHyperslabCoordinatesToMortonIndices(hyperslabs_, initial_layout_, global_domain_);
            hyperslabs_.clear();
        break;
        default:
            cmc_err_msg("The variable's coordinates cannot to transformed to Morton indices.");
    }

    active_format_ = DataFormat::LinearFormat;
}

template<class T>
bool
Variable<T>::IsValid() const
{
    if ((linear_indices_.size() > 0 && cartesian_coordinates_.size() > 0) ||
        (linear_indices_.size() > 0 && hyperslabs_.size() > 0) ||
        (cartesian_coordinates_.size() > 0 && hyperslabs_.size() > 0))
    {
        return false;
    }
    if (data_.size() != GetNumberCoordinates())
    {
        return false;
    }

    return true;
}

template<class T>
size_t
Variable<T>::GetIndexWithinDimension(const Dimension dimension, const size_t coordinate_position) const
{
    switch (GetActiveDataFormat())
    {
        case DataFormat::LinearFormat:
            cmc_err_msg("Currently not implemented.");
            return 0;
        break;
        case DataFormat::CartesianFormat:
            return cartesian_coordinates_[coordinate_position].GetDimensionCoordinate(dimension);
        break;
        case DataFormat::HyperslabFormat:
        {
            size_t count = 0;
            for (auto hs_iter = hyperslabs_.begin(); hs_iter != hyperslabs_.end(); ++hs_iter)
            {
                size_t current_count = 1;
                for (auto c_iter = hs_iter->CountIndicesBegin(); c_iter != hs_iter->CountIndicesEnd(); ++c_iter)
                {
                    current_count *= *c_iter;
                }
                /* Check if the coordinate position lies within the current hyperslab */
                if (coordinate_position < count + current_count)
                {
                    /* Get the hyperslab local index */
                    size_t pos = coordinate_position - count;
                    
                    const UpdateHyperslabCoordinateFn HsIndexIterationFuntion = GetHyperslabCoordinatesIterationFunction(initial_layout_);
                    std::vector<HyperslabIndex> reference_coordinates(GetDimensionalityOfDataLayout(initial_layout_));
                    HsIndexIterationFuntion(*hs_iter, reference_coordinates, pos);
                    const DimensionValueExtractionFn HsDimensionValueExtraction = GetDimensionValueFunctionForReferenceCoords(initial_layout_);
                    return HsDimensionValueExtraction(reference_coordinates, dimension);
                } else
                {
                    count += current_count;
                }
            }
            return CMC_ERR;
        }
        break;
        default:
            return 0;
    }
}

template<class T>
std::vector<T>&&
Variable<T>::DetachData()
{
    _has_been_invalidated_by_moving_ = true;
    return std::move(data_);
}

template<class T>
void
Variable<T>::ClearData()
{
    data_.clear();
}

template<class T>
std::vector<LinearIndex>&&
Variable<T>::DetachMortonIndices()
{
    _has_been_invalidated_by_moving_ = true;
    return std::move(linear_indices_);
}

template<class T>
void
Variable<T>::ClearMortonIndices()
{
    linear_indices_.clear();
}

template<class T>
T
Variable<T>::GetMissingValue() const
{
    return missing_value_;
}

template<class T>
void
Variable<T>::SetMissingValue(const CmcUniversalType& missing_value)
{
    cmc_assert(std::holds_alternative<T>(missing_value));
    missing_value_ = std::get<T>(missing_value);
}

template<class T>
void
Variable<T>::SetMissingValue(const T missing_value)
{
    missing_value_ = missing_value;
}

template<class T>
CmcUniversalType
Variable<T>::GetAddOffset() const
{
    return add_offset_;
}

template<class T>
CmcUniversalType
Variable<T>::GetScaleFactor() const
{
    return scale_factor_;
}

template<class T>
int
Variable<T>::GetID() const
{
    return id_;
}

template<class T>
const std::string&
Variable<T>::GetName() const
{
    return name_;
}

template<class T>
int
Variable<T>::GetGlobalContextInformation() const
{
    return global_context_information_;
}

template<class T>
void
Variable<T>::SetGlobalContextInformation(const int global_context_information)
{
    global_context_information_ = global_context_information;  
}

template<class T>
GeoDomain&&
Variable<T>::DetachGlobalDomain()
{
    return std::move(global_domain_);
}

template<class T>
CmcType
Variable<T>::GetType() const
{
    return type_;
}

template<class T>
std::vector<typename Variable<T>::ExtractedVar>
Variable<T>::ExtractAllSubVariablesByDimensionFromHyperslabs(const Dimension dimension) const
{
    const size_t num_extracted_variables = global_domain_.GetDimensionLength(dimension);
    
    std::vector<ExtractedVar> extracted_vars;
    extracted_vars.reserve(num_extracted_variables);

    /* Get the number of cooridnates for a single subvariable */
    HyperslabIndex _num_coordinates = 0;
    for (auto hs_iter = hyperslabs_.begin(); hs_iter != hyperslabs_.end(); ++hs_iter)
    {
        _num_coordinates += hs_iter->GetNumberCoordinatesWithoutCertainDimension(dimension);
    }

    const HyperslabIndex& num_coordinates = _num_coordinates;

    /* Reserve memory for each variable which will be extracted */
    for (size_t iter = 0; iter < num_extracted_variables; ++iter)
    {
        extracted_vars.emplace_back(ExtractedVar());
        extracted_vars.back().data.reserve(num_coordinates);
        extracted_vars.back().morton_indices.reserve(num_coordinates);
    }

    const int initial_dimensionality = GetDimensionalityOfDataLayout(initial_layout_);
    std::vector<HyperslabIndex> reference_coordinates(initial_dimensionality);
    const UpdateHyperslabCoordinateFn HsUpdateFn = GetHyperslabCoordinatesIterationFunction(initial_layout_);
    const DimensionValueExtractionFn ExtractVariableIDFn = GetDimensionValueFunctionForReferenceCoords(initial_layout_);
    const DataLayout resulting_layout = GetDataLayoutAfterRemoval(GetInitialDataLayout(), dimension);
    const int resulting_dimensionality  = GetDimensionalityOfDataLayout(resulting_layout);
    const TrimRefCoordVectorFn ReceiveTrimmedReferenceCoordsFn = GetReferenceCoordTrimmingFunction(resulting_layout);

    int data_accessor = 0;
    /* Iterate over all hyperslabs and assign each datum and coordinate to the correct sub variable */
    for (auto hs_iter = hyperslabs_.begin(); hs_iter != hyperslabs_.end(); ++hs_iter)
    {
        const HyperslabIndex hs_num_coordinates = hs_iter->GetNumberCoordinates();
        for (HyperslabIndex iter = 0; iter < hs_num_coordinates; ++iter)
        {
            /* The vector holding the refernce coordinates is updated to the current coordinate from the hyperslab */
            HsUpdateFn(*hs_iter, reference_coordinates, iter);
            /* The variable id is extracted from the reference coordinates */
            const HyperslabIndex extracted_var_id = ExtractVariableIDFn(reference_coordinates, dimension);
            /* Copy the corresponding data value to the new variable */
            extracted_vars[extracted_var_id].data.push_back(data_[data_accessor]);
            /* Transform the coordinate to a Morton index */
            extracted_vars[extracted_var_id].morton_indices.push_back(GetMortonIndex(ReceiveTrimmedReferenceCoordsFn(reference_coordinates), resulting_dimensionality));
            /* Update the counter to access the data vector of the input variable */
            ++data_accessor;
        }
    }

    return extracted_vars;
}

template<class T>
std::vector<typename Variable<T>::ExtractedVar>
Variable<T>::ExtractAllSubVariablesByDimension(const Dimension dimension) const
{
    switch (GetActiveFormat())
    {
        case DataFormat::LinearFormat:
            cmc_err_msg("Extraction for Linear indices is currently not implemented.");
            return std::vector<ExtractedVar>();
        break;
        case DataFormat::CartesianFormat:
            cmc_err_msg("Extraction for Cartesian Coordinates is currently not implemented.");
            return std::vector<ExtractedVar>();
        break;
        case DataFormat::HyperslabFormat:
            return ExtractAllSubVariablesByDimensionFromHyperslabs(dimension);
        break;
        default:
            cmc_err_msg("A not recognized data format is present. Therefore, the variable's data cannot be extracted.");
            return std::vector<ExtractedVar>();
    }
}

template<class T>
Variable<T>
HollowCopy(const Variable<T>& variable)
{
    Variable<T> hollow_variable(variable.name_, variable.id_, variable.initial_layout_);
    hollow_variable.missing_value_ = variable.missing_value_;
    hollow_variable.add_offset_ = variable.add_offset_;
    hollow_variable.scale_factor_ = variable.scale_factor_;
    hollow_variable.active_format_ = variable.active_format_;
    hollow_variable.global_domain_ = variable.global_domain_;
    hollow_variable.global_context_information_ = variable.global_context_information_;
    hollow_variable.has_split_dimension_ = variable.has_split_dimension_;
    hollow_variable._has_been_invalidated_by_moving_ = variable._has_been_invalidated_by_moving_;
    hollow_variable.pre_compression_layout_ = variable.pre_compression_layout_;
#ifdef CMC_ENABLE_MPI
    hollow_variable.comm_ = variable.comm_;
#endif
    hollow_variable.intern_id_ = variable.intern_id_;
    hollow_variable.cartesian_coordinates_ = variable.cartesian_coordinates_;
    hollow_variable.linear_indices_ = variable.linear_indices_;
    hollow_variable.hyperslabs_ = variable.hyperslabs_;
    hollow_variable.type_ = variable.type_;
    return hollow_variable;
}

template<class T>
std::vector<Variable<T>>
ExtractSubVariables(const Variable<T>& variable, const Dimension split_dimension)
{
    cmc_assert(variable.IsValid());

    const size_t num_extracted_variables = variable.global_domain_.GetDimensionLength(split_dimension);

    std::vector<Variable<T>> extracted_variables(num_extracted_variables, HollowCopy(variable));

    std::vector<typename Variable<T>::ExtractedVar> extracted_variable_data = variable.ExtractAllSubVariablesByDimension(split_dimension);

    for (size_t iter = 0; iter < num_extracted_variables; ++iter)
    {
        extracted_variables[iter].global_context_information_ = iter;

        extracted_variables[iter].data_ = std::move(extracted_variable_data[iter].data);

        extracted_variables[iter].linear_indices_ = std::move(extracted_variable_data[iter].morton_indices);

        extracted_variables[iter].initial_layout_ = GetDataLayoutAfterRemoval(variable.initial_layout_, split_dimension);
        
        extracted_variables[iter].pre_compression_layout_ = variable.initial_layout_;

        extracted_variables[iter].active_format_ = DataFormat::LinearFormat;
        
        extracted_variables[iter].global_domain_.ClearDimension(split_dimension);

        extracted_variables[iter].has_split_dimension_ = split_dimension;
#ifdef CMC_ENABLE_MPI
        extracted_variables[iter].comm_ = variable.GetMPIComm();
#endif
    }

    return extracted_variables;
}

template<class T>
int
Variable<T>::GetInternID() const
{
    return intern_id_;
}

template<class T>
void
Variable<T>::SetInternID(const int id)
{
    intern_id_ = id;
}

template<class T>
DataFormat
Variable<T>::GetActiveDataFormat() const
{
    return active_format_;
}

/** INPUT VAR MEMBER FUNCTIONS **/
template<typename T>
Var::Var(const Variable<T>& var)
: var_{var}{}

template<typename T>
Var::Var(Variable<T>&& var)
: var_{std::move(var)}{}

template<typename T>
bool
Var::holds_alternative() const
{
    return std::holds_alternative<T>(var_);
}

template<class T>
void 
Var::ApplyAxpyScalingAndOffset(const Variable<T>& variable)
{
    if (std::holds_alternative<T>(variable.GetScaleFactor()) &&
        std::holds_alternative<T>(variable.GetAddOffset()))
    {
        const T scale_factor = std::get<T>(variable.GetScaleFactor());
        const T add_offset = std::get<T>(variable.GetAddOffset());
        const T missing_value = variable.GetMissingValue();

        /* The data type remains unchanged */
        Variable<T> transformed_variable = HollowCopy(variable);

        /* Get the current data of te variable (as read only) */
        const std::vector<T>& source_data = variable.GetData(AccessKey());

        /* Create new data vector which will hold the transformed values */
        std::vector<T> destination_data;
        destination_data.reserve(source_data.size());

        for(auto iter = source_data.begin(); iter != source_data.end(); ++iter)
        { 
            if (!ApproxCompare(*iter, missing_value))
            {
                destination_data.push_back(scale_factor * (*iter) + add_offset);
            } else
            {
                destination_data.push_back(missing_value);
            }
        }

        /* Set the transformed data of the variable within the new variable */
        transformed_variable.SetData(AccessKey(), std::move(destination_data));

        var_ = GeneralVariable(std::move(transformed_variable));
    } else
    {
        /* The data types do not coincide and therefore they are transformed to the default data type */
        const CmcDefaultDataType scale_factor = GetUniversalDataAs<CmcDefaultDataType>(variable.GetScaleFactor());
        const CmcDefaultDataType add_offset = GetUniversalDataAs<CmcDefaultDataType>(variable.GetAddOffset());

        const T initial_missing_value = variable.GetMissingValue();
        const CmcDefaultDataType missing_value = static_cast<CmcDefaultDataType>(initial_missing_value);

        Variable<CmcDefaultDataType> transformed_variable = variable.template ObtainMetaCopy<CmcDefaultDataType>();

        /* Get the current data of the variable (as read only) */
        const std::vector<T>& source_data = variable.GetData(AccessKey());

        /* Create new data vector which will hold the transformed values */
        std::vector<CmcDefaultDataType> destination_data;
        destination_data.reserve(source_data.size());

        for (auto iter = source_data.begin(); iter != source_data.end(); ++iter)
        {
            if (!ApproxCompare(*iter, initial_missing_value))
            {
                destination_data.push_back(scale_factor * static_cast<CmcDefaultDataType>(*iter) + add_offset);
            } else
            {
                destination_data.push_back(missing_value);
            }
        }
        
        /* Set the transformed data of the variable within the new variable */
        transformed_variable.SetData(AccessKey(), std::move(destination_data));

        var_ = GeneralVariable(std::move(transformed_variable));
    }
}

template<class T>
void 
Var::ApplyScaling(const Variable<T>& variable)
{
    if (std::holds_alternative<T>(variable.GetScaleFactor()))
    {
        //TODO: Apply transformation in palce and do not create a new variable
        /* The data type remains unchanged */
        const T scale_factor = std::get<T>(variable.GetScaleFactor());
        const T missing_value = variable.GetMissingValue();

        Variable<T> transformed_variable = HollowCopy(variable);

        /* Get the current data of te variable (as read only) */
        const std::vector<T>& source_data = variable.GetData(AccessKey());

        /* Create new data vector which will hold the transformed values */
        std::vector<T> destination_data;
        destination_data.reserve(source_data.size());

        for(auto iter = source_data.begin(); iter != source_data.end(); ++iter)
        { 
            if (!ApproxCompare(*iter, missing_value))
            {
                destination_data.push_back(scale_factor * (*iter));
            } else
            {
                destination_data.push_back(missing_value);
            }
        }

        /* Set the transformed data of the variable within the new variable */
        transformed_variable.SetData(AccessKey(), std::move(destination_data));

        var_ = GeneralVariable(std::move(transformed_variable));
    } else
    {
        /* The data types do not coincide and therefore they are transformed to the default data type */
        const CmcDefaultDataType scale_factor = GetUniversalDataAs<CmcDefaultDataType>(variable.GetScaleFactor());

        const T initial_missing_value = variable.GetMissingValue();
        const CmcDefaultDataType missing_value = static_cast<CmcDefaultDataType>(initial_missing_value);

        Variable<CmcDefaultDataType> transformed_variable = variable.template ObtainMetaCopy<CmcDefaultDataType>();

        /* Get the current data of te variable (as read only) */
        const std::vector<T>& source_data = variable.GetData(AccessKey());

        /* Create new data vector which will hold the transformed values */
        std::vector<CmcDefaultDataType> destination_data;
        destination_data.reserve(source_data.size());

        for (auto iter = source_data.begin(); iter != source_data.end(); ++iter)
        {
            if (!ApproxCompare(*iter, initial_missing_value))
            {
                destination_data.push_back(scale_factor * static_cast<CmcDefaultDataType>(*iter));
            } else
            {
                destination_data.push_back(missing_value);
            }
        }

        /* Set the transformed data of the variable within the new variable */
        transformed_variable.SetData(AccessKey(), std::move(destination_data));
        
        var_ = GeneralVariable(std::move(transformed_variable));
    }
}

template<class T>
void 
Var::ApplyOffset(const Variable<T>& variable)
{
    if (std::holds_alternative<T>(variable.GetAddOffset()))
    {
        //TODO: Apply transformation in palce and do not create a new variable
        /* The data type remains unchanged */
        const T add_offset = std::get<T>(variable.GetAddOffset());
        const T missing_value = variable.GetMissingValue();

        /* The data type remains unchanged */
        Variable<T> transformed_variable = HollowCopy(variable);

        /* Get the current data of te variable (as read only) */
        const std::vector<T>& source_data = variable.GetData(AccessKey());

        /* Create new data vector which will hold the transformed values */
        std::vector<T> destination_data;
        destination_data.reserve(source_data.size());

        for(auto iter = source_data.begin(); iter != source_data.end(); ++iter)
        {
            if (!ApproxCompare(*iter, missing_value))
            {
                destination_data.push_back((*iter) + add_offset);
            } else
            {
                destination_data.push_back(missing_value);
            }
        }

        /* Set the transformed data of the variable within the new variable */
        transformed_variable.SetData(AccessKey(), std::move(destination_data));

        var_ = GeneralVariable(std::move(transformed_variable));
    } else
    {
        /* The data types do not coincide and therefore they are transformed to the default data type */
        const CmcDefaultDataType add_offset = GetUniversalDataAs<CmcDefaultDataType>(variable.GetAddOffset());

        const T initial_missing_value = variable.GetMissingValue();
        const CmcDefaultDataType missing_value = static_cast<CmcDefaultDataType>(initial_missing_value);

        Variable<CmcDefaultDataType> transformed_variable = variable.template ObtainMetaCopy<CmcDefaultDataType>();

        /* Get the current data of te variable (as read only) */
        const std::vector<T>& source_data = variable.GetData(AccessKey());

        /* Create new data vector which will hold the transformed values */
        std::vector<CmcDefaultDataType> destination_data;
        destination_data.reserve(source_data.size());

        for (auto iter = source_data.begin(); iter != source_data.end(); ++iter)
        {
            if (!ApproxCompare(*iter, initial_missing_value))
            {
                destination_data.push_back((static_cast<CmcDefaultDataType>(*iter) + add_offset));
            } else
            {
                destination_data.push_back(missing_value);
            }
        }

        /* Set the transformed data of the variable within the new variable */
        transformed_variable.SetData(AccessKey(), std::move(destination_data));

        var_ = GeneralVariable(std::move(transformed_variable));
    }
}

template<typename T>
void Variable<T>::AssignDataAtLinearIndices(const Var& var, const UpdateLinearIndices& update)
{
    cmc_assert(var.holds_alternative<Variable<T>>());

    const Variable<T>& source_var = std::get<Variable<T>>(var.GetInternalVariant());

    size_t index = 0;
    for (auto index_iter = source_var.linear_indices_.begin(); index_iter != source_var.linear_indices_.end(); ++index_iter, ++index)
    {
        data_[update(*index_iter)] = source_var.data_[index];
    }
}

template<typename T>
void Variable<T>::AssignDataAtLinearIndices(const Variable<T>& source_var, const UpdateLinearIndices& update)
{
    size_t index = 0;
    for (auto index_iter = source_var.linear_indices_.begin(); index_iter != source_var.linear_indices_.end(); ++index_iter, ++index)
    {
        data_[update(*index_iter)] = source_var.data_[index];
    }
}

template<typename T>
std::vector<T>
Variable<T>::GetDataFromHyperslab(const Hyperslab& hyperslab) const
{
    /* Get the global domain as hyperslab */
    const Hyperslab global_domain = TransformGeoDomainToHyperslab(global_domain_);

    /* Get the correct function to extract the linear indices for this hyperslab */
    LinearIndicesExtractionFn indices_extraction_fn = GetIndicesExtractionFunction(initial_layout_);

    /* Get the indices to extract the correct data */
    std::vector<HyperslabIndex> indices = indices_extraction_fn(global_domain, hyperslab);

    /* Allocate a vector which will hold the data; potentially padded with missing values */
    std::vector<T> data;
    data.reserve(indices.size());
    //cmc_debug_msg("\n\nIn Getdata from hyperslab");
    //cmc_debug_msg("Size of data in GetDatafromHyperslab: ", data_.size());
    /* Iterate over all indices and extract the data */
    for (auto index_iter = indices.begin(); index_iter != indices.end(); ++index_iter)
    {
        //cmc_debug_msg("IndexIterationVal: ", *index_iter);
        if (*index_iter != kOutsideOfHyperslabDomain)
        {
            data.push_back(data_[*index_iter]);
        } else
        {
            /* If the index is outside the global domain, we will assign a missing value */
            data.push_back(missing_value_);
        }
    }
    
    return data;
}

/** UPDATE LINEAR INDICES MEMBER FUNCTIONS **/
template<typename T>
auto UpdateLinearIndices::operator()(T index_to_be_updated) const
 -> std::enable_if_t<std::is_same_v<std::decay_t<T>, MortonIndex>, MortonIndex>
{
    auto previous_index_correction = std::upper_bound(correction_.begin(), correction_.end(), index_to_be_updated,
                                                      [](const MortonIndex& val1, const IndexReduction& val2) -> bool
                                                      {
                                                        return (val1 < val2.start_index);
                                                      });
    
    cmc_assert(previous_index_correction != correction_.begin());

    /* The binary search yields the upper bound of the Morton Indices, therefore we need to take a step back
     * in order to obtain the previous uniform indices which were skipped (due to the possible coarse element within the mesh) */
    const MortonIndex correction = std::prev(previous_index_correction)->indices_to_subtract;

    return index_to_be_updated - correction;
}


#ifdef CMC_ENABLE_MPI

template<class T>
void
Variable<T>::CreateInternIDForRedistribution()
{
    if (GetGlobalContextInformation() != kNoGlobalContext)
    {
        /* If the variable is part of a higher dimensional one and resembles a subband of the variable */
        const int intern_id = CreateSubbandID(GetID(), GetGlobalContextInformation());
        SetInternID(intern_id);
    } else
    {
        /* If the variable has not been split */
        const int intern_id = GetID();
        SetInternID(intern_id);
    }
}

template <typename T>
void
AppendSendData(std::vector<VariableSendMessage>& messages, ReceiverMap<T>&& send_data)
{
    for (auto sd_iter = send_data.begin(); sd_iter != send_data.end();)
    {
        messages.push_back(std::move(sd_iter->second));
        sd_iter = send_data.erase(sd_iter);
    }
}

template<class T>
ReceiverMap<T>
Variable<T>::GatherDataToBeDistributed(const DataOffsets& offsets)
{
    CreateInternIDForRedistribution();
    
    //TransformCoordinatesToMortonIndices();
    
    cmc_assert(GetActiveDataFormat() == DataFormat::LinearFormat);

    ReceiverMap<T> send_messages;

    const size_t num_coordinates = GetNumberCoordinates();
    cmc_debug_msg("Num coordinates in gather data to be distributed: ", num_coordinates);

    LinearIndex previous_lower_bound = 0;
    LinearIndex previous_upper_bound = 0;
    int owner_rank = -1;

    for (size_t index = 0; index < num_coordinates; ++index)
    {
        const LinearIndex current_linear_index = GetLinearIndexCoordinate(index);

        /* Check if the current linear index belongs to the same rank as the previous one, if not we need to find the new receiving rank */
        if (!(previous_lower_bound <= current_linear_index && current_linear_index < previous_upper_bound))
        {
            /* Find an iterator to the rank which is ought to hold the value corresponding to this linear index */
            auto owner_rank_iter = std::upper_bound(offsets.Begin(), offsets.End(), current_linear_index);

            /* Get the integer number of the corresponding rank */
            owner_rank = std::distance(offsets.Begin(), owner_rank_iter) - 1;
            
            previous_lower_bound = offsets[owner_rank];
            previous_upper_bound = offsets[owner_rank + 1];

            /* Check if the rank is already a receiving rank */
            if (auto search_rk = send_messages.find(owner_rank); search_rk != send_messages.end())
            {
                /* If the rank is already listed in the ReceiverMap */
                search_rk->second.data_.push_back(data_[index]);
                search_rk->second.morton_indices_.push_back(current_linear_index);
            } else
            {
                cmc_assert(GetInternID() != kNoInternalIDSet);
                /* If the receiving rank is not yet listed within the ReceiverMap */
                send_messages[owner_rank] = VariableMessage<T>(owner_rank, GetInternID());

                const size_t estimate_receiving_data_points = 2 * (num_coordinates / offsets.size());

                /* Reserve an estimate of data */
                send_messages[owner_rank].data_.reserve(estimate_receiving_data_points);
                send_messages[owner_rank].morton_indices_.reserve(estimate_receiving_data_points);
                
                /* Save the value and index */
                send_messages[owner_rank].data_.push_back(data_[index]);
                send_messages[owner_rank].morton_indices_.push_back(current_linear_index);
            }
        } else
        { 
            /* If the current value belongs to same rank the previous value already belongs to */
            send_messages[owner_rank].data_.push_back(data_[index]);
            send_messages[owner_rank].morton_indices_.push_back(current_linear_index);
        }
    }

    return send_messages;
}


template<typename T>
void Variable<T>::AssignDataAtLinearIndices(const VariableRecvMessage& message, const UpdateLinearIndices& update)
{
    cmc_assert(std::holds_alternative<VariableMessage<T>>(message.GetInternalVariant()));

    /* Get the message holding the actual data and their Morton idnices from the message */
    const VariableMessage<T>& msg = std::get<VariableMessage<T>>(message.GetInternalVariant());

    /* Iterate throught the indices and the data and assign it accordingly */
    auto data_iter = msg.DataBegin();
    for (auto morton_idx_iter = msg.MortonIndicesBegin(); morton_idx_iter != msg.MortonIndicesEnd(); ++morton_idx_iter, ++data_iter)
    {
        //cmc_debug_msg("We will write at pos: ", update(*morton_idx_iter), " the val: ", *data_iter);
        data_[update(*morton_idx_iter)] = *data_iter;
    }
}

#endif

}



#endif /* !INPUT_CMC_INPUT_VARIABLE_HXX */
