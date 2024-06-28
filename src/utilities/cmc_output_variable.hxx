#ifndef CMC_OUTPUT_VARIABLE_HXX
#define CMC_OUTPUT_VARIABLE_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_output_variable_forward.hxx"
#include "utilities/cmc_variable_transformer_forward.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_dimension_interval.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "t8code/cmc_t8_morton.hxx"
#include "utilities/cmc_cart_coordinate.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "t8code/cmc_t8_mpi.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy_forward.hxx"
#include "t8code/cmc_t8_interpolation.hxx"

#include <iterator>
#include <type_traits>

namespace cmc
{

template<class T>
class OutputVariable
{
public:
    OutputVariable() = default;
    OutputVariable(const std::string& name, const DataLayout layout, const GeoDomain& domain)
    : name_{name}, layout_{layout}, global_domain_{domain}{};

    void AllocateData(const int size_hint);

    void SetMissingValue(const CmcUniversalType missing_value);

    DataLayout GetDataLayout() const;

    void PushBackHyperslab(Hyperslab&& hyperslab);

    void PushBackHyperslab(const Hyperslab& hyperslab);

    void AssignDataToHyperslab(const CmcUniversalType& value, const Hyperslab& hyperslab, const GeoDomain& resulting_domain);

    std::vector<T> GetData() const;

    std::vector<T>&& GetData();

    friend class TransformerCompressionToOutputVariable;
private:
    std::string name_;
    std::vector<T> data_;
    T missing_value_;
    DataLayout layout_;
    std::vector<Hyperslab> hyperslabs_;
    GeoDomain global_domain_;
};

class OutputVar
{
public:
    OutputVar() = default;
    OutputVar(const int id)
    : id_{id}{};
    OutputVar(const int id, CmcOutputVariable&& variable)
    : id_{id}, var_{std::move(variable)}{};
    OutputVar(const int id, const CmcType type, const int size_hint, const std::string& name, const DataLayout decompressed_layout, const GeoDomain& global_domain,const CmcUniversalType missing_value)
    : id_{id} {
        SetupOutputVar(type, size_hint, name, decompressed_layout, global_domain, missing_value);
    };

    void AssignDataToHyperslab(const CmcUniversalType& value, const Hyperslab& hyperslab, const GeoDomain& resulting_domain);

    void PushBackHyperslab(Hyperslab&& hyperslab);

    void PushBackHyperslab(const Hyperslab& hyperslab);

    template<typename T> 
    OutputVariable<T>&& SeizeOutputVariable();

private:
    void SetupOutputVar(const CmcType type, const int size_hint, const std::string& name, const DataLayout decompressed_layout, const GeoDomain& global_domain,const CmcUniversalType missing_value);
    
    int id_;
    CmcOutputVariable var_;
};


/** OUTPUT VARIABLE<T> MEMBER FUNCTIONS **/
template<class T>
void
OutputVariable<T>::AllocateData(const int size_hint)
{
    data_ = std::vector<T>(size_hint, 0);
}

template<class T>
void
OutputVariable<T>::SetMissingValue(const CmcUniversalType missing_value)
{
    cmc_assert(std::holds_alternative<T>(missing_value));
    missing_value_ = std::get<T>(missing_value);
}

template<class T>
DataLayout
OutputVariable<T>::GetDataLayout() const 
{
    return layout_;
}

template<class T>
void
OutputVariable<T>::PushBackHyperslab(Hyperslab&& hyperslab)
{
    hyperslabs_.push_back(std::move(hyperslab));
}

template<class T>
void
OutputVariable<T>::PushBackHyperslab(const Hyperslab& hyperslab)
{
    hyperslabs_.push_back(hyperslab);
}

template<class T>
void
OutputVariable<T>::AssignDataToHyperslab(const CmcUniversalType& value, const Hyperslab& hyperslab, const GeoDomain& resulting_domain)
{
    cmc_assert(std::holds_alternative<T>(value));
    const T val = std::get<T>(value);
    const UpdateHyperslabCoordinateFn hs_iteration_fn = GetHyperslabCoordinatesIterationFunction(layout_);
    const HyperslabIndex num_coordinates = hyperslab.GetNumberCoordinates();
    const int initial_dimensionality = GetDimensionalityOfDataLayout(layout_);
    std::vector<HyperslabIndex> linear_indices(initial_dimensionality, 0);
    TransformCartesianToIndexAccessorFn data_accessor_fn = GetCartesianCoordsToLinearIndexFunction(layout_);
    for (HyperslabIndex iter = 0; iter < num_coordinates; ++iter)
    {
        /* Get the current coordinates of the hyperslab */
        hs_iteration_fn(hyperslab, linear_indices, iter);
        /* Get an index regarding the data layout of the data in order to assign the value */
        data_[data_accessor_fn(linear_indices, resulting_domain)] = val;
    }
}

template<class T>
std::vector<T>
OutputVariable<T>::GetData() const
{
    return data_;
}

template<class T>
std::vector<T>&&
OutputVariable<T>::GetData()
{
    return std::move(data_);
}

/** OUTPUT VAR MEMBER FUNTIONS **/
inline void 
OutputVar::AssignDataToHyperslab(const CmcUniversalType& value, const Hyperslab& hyperslab, const GeoDomain& resulting_domain)
{
    std::visit([&](auto&& var){
        var.AssignDataToHyperslab(value, hyperslab, resulting_domain);
    }, var_);
}

inline void
OutputVar::PushBackHyperslab(Hyperslab&& hyperslab)
{
    std::visit([&](auto&& var){
        var.PushBackHyperslab(std::forward<Hyperslab>(hyperslab));
    }, var_);
}

inline void
OutputVar::PushBackHyperslab(const Hyperslab& hyperslab)
{
    std::visit([&](auto&& var){
        var.PushBackHyperslab(hyperslab);
    }, var_);
}

template<typename T> 
OutputVariable<T>&&
OutputVar::SeizeOutputVariable()
{
    if (!std::holds_alternative<OutputVariable<T>>(var_))
    {
        cmc_err_msg("The variable is of a different data type.");
    }
    return std::move(std::get<OutputVariable<T>>(var_));
}

}

#endif /* !CMC_OUTPUT_VARIABLE_HXX */
