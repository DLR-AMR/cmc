#ifndef CMC_T8_DATA_VARIABLES_HXX
#define CMC_T8_DATA_VARIABLES_HXX
/**
 * @file cmc_t8_data_variables.hxx
 */
#include "utilities/cmc_utilities.hxx"

#include <cstdlib>
#include <memory>

namespace cmc {

class DataArray
{
public:
    DataArray() = delete;
    constexpr DataArray(const CmcType type_of_datum, const size_t number_of_data_points)
    : data_type{type_of_datum}, num_elements{number_of_data_points}
    {
        data = std::malloc(CmcTypeToBytes(type_of_datum) * number_of_data_points);
    };

    CmcType GetDataType() const;
    size_t GetNumberElements() const;

    void Axpy(const CmcUniversalType sclaing_factor, const CmcUniversalType offset);
    void Scale(const CmcUniversalType scaling_factor);
    void Add(const CmcUniversalType constant_summand);

    void Resize(const size_t number_of_data_points);

private:
    const CmcType data_type{CmcType::TypeUndefined};
    const size_t num_elements{0};
    void* data{nullptr};
};

/* The global domain cooridnates are saved, therefore, there is no need to store the local offset */
class VariableDomain
{
public:
    enum DataRepresentation {RepresentationUndefined, SpaceFillingCurve, CartesianCoordinates, Hyperslab};

    VariableDomain() = delete;
    VariableDomain(DataRepresentation data_representation)
    : current_domain_representation{data_representation}{};

    void TransformToMortonIndices();
    void TransformToCartesianCoordinates();

    //Not sure if iterators are needed
    auto GetLinearIndicesIteratorBegin() const;
    auto GetLinearIndicesIteratorEnd() const;

    auto GetCartesianCoordinatesIteratorBegin() const;
    auto GetCartesianCoordinatesIteratorEnd() const;

    auto GetHyperslabIteratorBegin() const;
    auto GetHyperslabIteratorEnd() const;

    void PushBackCoordinate(const CartesianCoordinate& cartesian_coordinate);
    void PushBackCoordinate(CartesianCoordinate&& cartesian_coordinate);
    void PushBackCoordinate(int64_t linear_index);
    void PushBackCoordinate(const Hyperslab& hyperslab);
    void PushBackCoordinate(Hyperslab&& hyperslab);

    void ClearDomain();
    
private:
    DataRepresentation current_domain_representation{RepresentationUndefined};
    std::vector<CartesianCoordinate> cartesian_coordinates;
    std::vector<int64_t> linear_indices;
    std::vector<Hyperslab> hyperslab_coordinates;
};

class VariableAttributes
{
public:

private:
    CmcUniversalType missing_value;
    CmcUniversalType add_offset{static_cast<double>(0.0)};
    CmcUniversalType scale_factor{static_cast<double>(1.0)};
    bool is_scaling_and_offset_applied{true};
    int global_context_information{0}; //!< A variable for describing an optional additional relation (the meaning of this varibale may vary)
};


class Variable
{
public:
    Variable(){};
    ~Variable(){};

    std::string GetName() const;
    bool isValidForCompression() const;

private:

    std::string name{}; //!< The name of the variable
    std::unique_ptr<DataArray> data; //!< The actual data of the variable
    AmrMesh mesh; //!< The mesh on which the variable is defined
    VariableAttributes attributes; //!< Application specific attributes for the variable
    VariableDomain local_domain; //!< The domain on which the data is defined
    GlobalDomain global_domain; //!< The specifications about the global domain
};


}

#endif /* !CMC_T8_DATA_VARIABLES_HXX */
