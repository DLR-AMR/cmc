#ifndef CMC_NC_WRITER_HXX
#define CMC_NC_WRITER_HXX

#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "mpi/cmc_mpi.hxx"

#include <vector>

namespace cmc::nc
{

class Writer
{
public:
    Writer() = delete;
    Writer(const std::string& file_name, const int netcdf_format, const MPI_Comm comm = MPI_COMM_SELF)
    : file_name_{file_name}, netcdf_format_{netcdf_format}, comm_{comm} {};
    ~Writer() = default;

    Writer(const Writer& other) = default;
    Writer& operator=(const Writer& other) = default;
    Writer(Writer&& other) = default;
    Writer& operator=(Writer&& other) = default;

    void ReserveVariables(const size_t num_variables);
    void AddGlobalAttribute(const Attribute& attribute);
    void AddGlobalAttribute(Attribute&& attribute);

    void AddVariable(const Variable& variable);
    void AddVariable(Variable&& variable);
    template<typename T> void AddVariable(const SpecificVariable<T>& variable, const std::vector<Attribute>& attributes);
    template<typename T> void AddVariable(const SpecificVariable<T>& variable, std::vector<Attribute>&&  attributes);
    template<typename T> void AddVariable(SpecificVariable<T>&& variable, const std::vector<Attribute>& attributes);
    template<typename T> void AddVariable(SpecificVariable<T>&& variable, std::vector<Attribute>&& attributes);

    void ClearVariablesAndGlobalAttributes();

    void Write();

private:
    int Open();
    int Create();
    void Close(const int ncid) const;
    std::vector<int> DefineVariableDimensions(const int ncid, const std::vector<Dimension>& dimensions);
    void DefineAttributes(const int ncid, const int var_id, const std::vector<Attribute>& attributes);
    std::vector<int> DefineVariables(const int ncid);
    void DefineGlobalAttributes(const int ncid);
    void WriteData(const int ncid, const std::vector<int>& var_ids);

    const std::string file_name_;
    const int netcdf_format_;
    const MPI_Comm comm_;
    std::vector<Variable> variables_;
    std::vector<Attribute> global_attributes_;
    bool file_has_been_created_{false};
};

template<typename T>
void
Writer::AddVariable(const SpecificVariable<T>& variable, const std::vector<Attribute>& attributes)
{
    AddVariable(Variable(variable, attributes));
}

template<typename T>
void
Writer::AddVariable(const SpecificVariable<T>& variable, std::vector<Attribute>&&  attributes)
{
    AddVariable(Variable(variable, std::move(attributes)));
}

template<typename T>
void
Writer::AddVariable(SpecificVariable<T>&& variable, const std::vector<Attribute>& attributes)
{
    AddVariable(Variable(std::move(variable), attributes));
}

template<typename T>
void
Writer::AddVariable(SpecificVariable<T>&& variable, std::vector<Attribute>&& attributes)
{
    AddVariable(Variable(std::move(variable), std::move(attributes)));
}

}

#endif /* !CMC_NC_WRITER_HXX */
