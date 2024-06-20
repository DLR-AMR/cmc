#ifndef CMC_NC_WRITER_HXX
#define CMC_NC_WRITER_HXX

#include "cmc.h"
#include "utilities/cmc_utilities.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "mpi/cmc_mpi.hxx"

#include <vector>

namespace cmc
{

class NcWriter
{
public:
    NcWriter() = delete;
    NcWriter(const std::string& file_name, const int netcdf_format, const MPI_Comm comm = MPI_COMM_SELF)
    : file_name_{file_name}, netcdf_format_{netcdf_format}, comm_{comm} {};
    ~NcWriter() = default;

    NcWriter(const NcWriter& other) = default;
    NcWriter& operator=(const NcWriter& other) = default;
    NcWriter(NcWriter&& other) = default;
    NcWriter& operator=(NcWriter&& other) = default;

    void AddGlobalAttribute(const NcAttribute& attribute);
    void AddGlobalAttribute(NcAttribute&& attribute);

    void AddVariable(const NcVariable& variable);
    void AddVariable(NcVariable&& variable);
    template<typename T> void AddVariable(const NcSpecificVariable<T>& variable, const std::vector<NcAttribute>& attributes);
    template<typename T> void AddVariable(const NcSpecificVariable<T>& variable, std::vector<NcAttribute>&&  attributes);
    template<typename T> void AddVariable(NcSpecificVariable<T>&& variable, const std::vector<NcAttribute>& attributes);
    template<typename T> void AddVariable(NcSpecificVariable<T>&& variable, std::vector<NcAttribute>&& attributes);

    void Write();

private:
    int NcOpen();
    int NcCreate();
    void NcClose(const int ncid) const;
    std::vector<int> DefineVariableDimensions(const int ncid, const std::vector<NcDimension>& dimensions);
    void DefineAttributes(const int ncid, const int var_id, const std::vector<NcAttribute>& attributes);
    std::vector<int> DefineVariables(const int ncid);
    void DefineGlobalAttributes(const int ncid);
    void WriteData(const int ncid, const std::vector<int>& var_ids);

    const std::string file_name_;
    const int netcdf_format_;
    const MPI_Comm comm_;
    std::vector<NcVariable> variables_;
    std::vector<NcAttribute> global_attributes_;
    bool file_has_been_created_{false};
};

template<typename T>
void
NcWriter::AddVariable(const NcSpecificVariable<T>& variable, const std::vector<NcAttribute>& attributes)
{
    AddVariable(NcVariable(variable, attributes));
}

template<typename T>
void
NcWriter::AddVariable(const NcSpecificVariable<T>& variable, std::vector<NcAttribute>&&  attributes)
{
    AddVariable(NcVariable(variable, std::move(attributes)));
}

template<typename T>
void
NcWriter::AddVariable(NcSpecificVariable<T>&& variable, const std::vector<NcAttribute>& attributes)
{
    AddVariable(NcVariable(std::move(variable), attributes));
}

template<typename T>
void
NcWriter::AddVariable(NcSpecificVariable<T>&& variable, std::vector<NcAttribute>&& attributes)
{
    AddVariable(NcVariable(std::move(variable), std::move(attributes)));
}

}

#endif /* !CMC_NC_WRITER_HXX */
