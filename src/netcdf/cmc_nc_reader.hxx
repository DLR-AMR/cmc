#ifndef CMC_NC_READER_HXX
#define CMC_NC_READER_HXX

#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "mpi/cmc_mpi.hxx"

#include <vector>
#include <cstring>

namespace cmc
{

class NcReader
{
public:
    NcReader() = delete;
    ~NcReader() = default;

    NcReader(const NcReader& other) = default;
    NcReader& operator=(const NcReader& other) = default;
    NcReader(NcReader&& other) = default;
    NcReader& operator=(NcReader&& other) = default;

    NcReader(const std::string& file_name, const MPI_Comm comm = MPI_COMM_SELF)
    : file_name_{file_name}, comm_{comm} {};

    /* Read and return the data of a single specified variable */
    template<typename T> std::vector<T> ReadVariableData(const std::string& variable_name);
    template<typename T> std::vector<T> ReadVariableData(const std::string& variable_name, const GeneralHyperslab& hyperslab);

    //TODO: Implement
    NcVariable ReadVariable(const std::string& variable_name);
    std::vector<NcAttribute> ReadVariableAttributes(const std::string& variable_name);


    void StashVariableForReading(const std::string& variable_name, const std::vector<GeneralHyperslab>& hyperslabs);
    void StashVariableForReading(const std::string& variable_name, std::vector<GeneralHyperslab>&& hyperslabs);
    void StashVariableForReading(const std::string& variable_name, const GeneralHyperslab& hyperslab);
    void StashVariableForReading(const std::string& variable_name, GeneralHyperslab&& hyperslab);

    void ClearStashedVariables();

    /* Read and return only the stashed variables and the global attributes */
    std::pair<std::vector<NcVariable>, std::vector<NcAttribute>> ReadVariables();

    /* Read and return all global attributes */
    std::vector<NcAttribute> ReadGlobalAttrtibutes();

    /* Read all the meta data of all variables and return those "variable hulls" */
    std::vector<NcVariable> ReadVariableMetaData();

    /* Read all the meta data of all variables and return those "variable hulls". Additionally, read and return all global attributes  */
    std::pair<std::vector<NcVariable>, std::vector<NcAttribute>> ReadVariableMetaDataAndGlobalAttributes();

    const std::string& GetFileName() const;
    int GetNetcdfFormat() const;
    int GetNumberOfDimensions() const;
    int GetNumberOfVariables() const;
    int GetNumberOfGlobalAttributes() const;
    int GetNumberOfUnlimitedDimensions() const;
    std::vector<int> GetUnlimitedDimensionIDs() const;

    NcVariable&& GetVariable(const std::string& variable_name);

    std::vector<NcDimension> GetVariableDimensions(const std::string& variable_name);

    CmcType GetTypeOfVariable(const std::string& variable_name);

private:
    struct StashedVariable;

    int NcOpen();
    void NcClose(const int ncid);
    int FindVariableID(const int ncid, const std::string& variable_name);
    void InquireGeneralFileInformation(const int ncid);
    std::vector<NcAttribute> InquireAttributes(const int ncid, const int var_id);
    std::vector<NcVariable> InquireVariableMetaData(const int ncid);
    std::vector<NcDimension> ConvertDimensionIDs(const int ncid, const std::vector<int>& dim_ids);
    void ReadVariableDataFromFile(const int ncid, const nc_type var_type, const std::string& var_name, const int var_id, const std::vector<GeneralHyperslab>& hyperslabs, NcVariable& variable);
    
    const std::string file_name_;
    const MPI_Comm comm_;

    int netcdf_format_{-1};
    int num_dimensions_{0};
    int num_variables_{0};
    int num_global_attributes_{0};
    std::vector<int> unlimited_dimension_ids_;

    std::vector<NcDimension> dimensions_;
    std::vector<NcVariable> variables_;
    std::vector<NcAttribute> global_attributes_;

    std::vector<StashedVariable> variable_stash_;
    bool is_file_opened_{false};
    bool has_general_information_been_inquired_{false};
};

struct NcReader::StashedVariable
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

}

#endif /* !CMC_NC_READER_HXX */
