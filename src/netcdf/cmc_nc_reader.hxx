#ifndef CMC_NC_READER_HXX
#define CMC_NC_READER_HXX

#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_netcdf.hxx"
#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#endif
#include <vector>
#include <cstring>

namespace cmc::nc
{

class Reader
{
public:
    Reader() = delete;
    ~Reader() = default;

    Reader(const Reader& other) = default;
    Reader& operator=(const Reader& other) = default;
    Reader(Reader&& other) = default;
    Reader& operator=(Reader&& other) = default;

    Reader(const std::string& file_name)
    : file_name_{file_name} {};
#ifdef CMC_ENABLE_MPI
    Reader(const std::string& file_name, const MPI_Comm comm)
    : file_name_{file_name}, comm_{comm} {};
#endif
    /* Read a single variable from the file */
    Variable ReadVariable(const std::string& variable_name);

    /* Read and return the data of a single specified variable */
    template<typename T> std::vector<T> ReadVariableData(const std::string& variable_name);
    template<typename T> std::vector<T> ReadVariableData(const std::string& variable_name, const GeneralHyperslab& hyperslab);

    /* Read only the attributes of single variable */
    std::vector<Attribute> ReadVariableAttributes(const std::string& variable_name);

    /* Read and return all global attributes */
    std::vector<Attribute> ReadGlobalAttrtibutes();

    /* Save variable names and hyperslabs which will be read together from the file on calling ReadVariables() */
    void StashVariableForReading(const std::string& variable_name, const std::vector<GeneralHyperslab>& hyperslabs);
    void StashVariableForReading(const std::string& variable_name, std::vector<GeneralHyperslab>&& hyperslabs);
    void StashVariableForReading(const std::string& variable_name, const GeneralHyperslab& hyperslab);
    void StashVariableForReading(const std::string& variable_name, GeneralHyperslab&& hyperslab);

    /* Remove all stored information from the stash */
    void ClearStashedVariables();

    /* Read and return only the stashed variables and the global attributes */
    std::pair<std::vector<Variable>, std::vector<Attribute>> ReadVariables();

    /* Read all the meta data of all variables and return those "variable hulls" */
    std::vector<Variable> ReadVariableMetaData();

    /* Read all the meta data of all variables and return those "variable hulls". Additionally, read and return all global attributes  */
    std::pair<std::vector<Variable>, std::vector<Attribute>> ReadVariableMetaDataAndGlobalAttributes();

    /* Get the dimensions of a single variable */
    std::vector<Dimension> GetVariableDimensions(const std::string& variable_name);
    
    /* Get the data type of a single variable */
    CmcType GetTypeOfVariable(const std::string& variable_name);

    /* INquire some basic information about the file which can be obtained by the Getter-functions afterwards */
    void InquireGeneralFileInformation();

    const std::string& GetFileName() const;
    int GetNetcdfFormat() const;
    int GetNumberOfDimensions() const;
    int GetNumberOfVariables() const;
    int GetNumberOfGlobalAttributes() const;
    int GetNumberOfUnlimitedDimensions() const;
    std::vector<int> GetUnlimitedDimensionIDs() const;

private:
    struct StashedVariable;

    int Open();
    void Close(const int ncid);
    int FindVariableID(const int ncid, const std::string& variable_name);
    GeneralHyperslab GetDataDomainAsGeneralHyperslab(const std::string& variable_name);
    GeneralHyperslab GetDataDomainAsGeneralHyperslab(const int ncid, const std::string& variable_name);
    void InquireGeneralFileInformation(const int ncid);
    std::vector<Attribute> InquireAttributes(const int ncid, const int var_id);
    std::vector<Variable> InquireVariableMetaData(const int ncid);
    std::vector<Dimension> ConvertDimensionIDs(const int ncid, const std::vector<int>& dim_ids);
    void ReadVariableDataFromFile(const int ncid, const nc_type var_type, const std::string& var_name, const int var_id, const std::vector<GeneralHyperslab>& hyperslabs, Variable& variable);
    
    const std::string file_name_;

    int netcdf_format_{-1};
    int num_dimensions_{0};
    int num_variables_{0};
    int num_global_attributes_{0};
    std::vector<int> unlimited_dimension_ids_;

    std::vector<Dimension> dimensions_;
    std::vector<Variable> variables_;
    std::vector<Attribute> global_attributes_;

    std::vector<StashedVariable> variable_stash_;
    bool is_file_opened_{false};
    bool has_general_information_been_inquired_{false};

#ifdef CMC_ENABLE_MPI
    const MPI_Comm comm_{MPI_COMM_SELF};
#endif
};

struct Reader::StashedVariable
{
    StashedVariable() = delete;

    StashedVariable(const std::string& variable_name, const std::vector<GeneralHyperslab>& variable_hyperslabs)
    : name{variable_name}, hyperslabs{variable_hyperslabs} {};
    StashedVariable(const std::string& variable_name, std::vector<GeneralHyperslab>&& variable_hyperslabs)
    : name{variable_name}, hyperslabs{std::move(variable_hyperslabs)} {};
    StashedVariable(const std::string& variable_name, const GeneralHyperslab& variable_hyperslab)
    : name{variable_name}, hyperslabs{variable_hyperslab} {};
    StashedVariable(const std::string& variable_name, GeneralHyperslab&& variable_hyperslab)
    : name{variable_name}, hyperslabs{std::move(variable_hyperslab)} {};


    const std::string name;
    const std::vector<GeneralHyperslab> hyperslabs;
};

template<typename T>
std::vector<T>
Reader::ReadVariableData(const std::string& variable_name, const GeneralHyperslab& hyperslab)
{
    /* Open the file to be read */
    const int ncid = Open();

    InquireGeneralFileInformation(ncid);

    /* Get the corresponding ID to the supplied variable name */
    const int var_id = FindVariableID(ncid, variable_name);

    /* We have found the correct variable with the given name.
     * Now, we are able to inquire the data type of the variable */
    nc_type var_type{NC_NAT};
    const int type_err = nc_inq_vartype(ncid, var_id, &var_type);
    CheckError(type_err);

    /* Get the Cmc Type of the variable */
    const CmcType type = ConvertNcTypeToCmcType(var_type);

    /* Compare whether the CmcType complies with the template type */
    if (type != ConvertToCmcType<T>())
    {
        cmc_err_msg("The variable ", variable_name, " is of a different type (CmcType: ", type, ") than the template type.");
    }

    /* Get the number of data values which will be read for this hyperslab */
    const size_t size = hyperslab.GetNumberOfCoveredCoordinates();

    /* Allocate a vector of the given size */
    std::vector<T> data(size);

    /* Read in the data */
    const int data_err = nc_get_vara(ncid, var_id, hyperslab.start_values.data(), hyperslab.count_values.data(), static_cast<void*>(data.data()));
    CheckError(data_err);

    /* Close the file after the reading process is finished */
    Close(ncid);

    return data;
}

template<typename T>
std::vector<T>
Reader::ReadVariableData(const std::string& variable_name)
{
    return ReadVariableData<T>(variable_name, GetDataDomainAsGeneralHyperslab(variable_name));
}

}

#endif /* !CMC_NC_READER_HXX */
