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
    NcWriter(const std::string& file_name, const MPI_Comm comm = MPI_COMM_SELF)
    : file_name_{file_name}, comm_{comm} {};
    ~NcWriter() = default;

    NcWriter(const NcWriter& other) = default;
    NcWriter& operator=(const NcWriter& other) = default;
    NcWriter(NcWriter&& other) = default;
    NcWriter& operator=(NcWriter&& other) = default;

    void AddGlobalAttribute(const NcAttribute& attribute);
    void AddGlobalAttribute(NcAttribute&& attribute);

    void AddVariable(const NcVariable& variable);
    void AddVariable(NcVariable&& variable);

    void Write();

private:
    int NcOpen() const;
    void NcClose(const int ncid) const;
    std::vector<int> DefineVariableDimensions(const int ncid, const std::vector<NcDimension>& dimensions);
    void DefineAttributes(const int ncid, const int var_id, const std::vector<NcAttribute>& attributes);
    void DefineVariables(const int ncid);
    void DefineGlobalAttributes(const int ncid);
    void WriteData(const int ncid);

    const std::string file_name_;
    const MPI_Comm comm_;
    std::vector<NcVariable> variables_;
    std::vector<NcAttribute> global_attributes_;
};

}

#endif /* !CMC_NC_WRITER_HXX */
