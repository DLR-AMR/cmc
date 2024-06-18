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

    void AddGlobalAttribute(const NcAttribute& attribute);
    void AddGlobalAttribute(NcAttribute&& attribute);

    void AddVariable(const NcVariable& variable);
    void AddVariable(NcVariable&& variable);

    void Write();

private:
    int NcOpen() const;
    void NcClose(const int ncid) const;

    const std::string file_name_;
    const MPI_Comm comm_;
    std::vector<NcVariable> variables_;
    std::vector<NcAttribute> global_attributes_;
};


int
NcWriter::NcOpen() const
{
    int ncid{-1};

    if (comm_ != MPI_COMM_SELF)
    {
        /* Open for parallel access */
        const MPI_Info info = MPI_INFO_NULL;
        const int err = nc_open_par(file_name_.c_str(), NC_WRITE, comm_, info, &ncid);
        NcCheckError(err);
    } else
    {
        /* Otherwise open for serial access */
        const int err = nc__open(file_name_.c_str(), NC_WRITE, NULL, &ncid);
        NcCheckError(err);
    }

    return ncid;
}

void
NcWriter::NcClose(const int ncid) const
{
    const int err = nc_close(ncid);
    NcCheckError(err);
}

void
NcWriter::Write()
{
    const int ncid = NcOpen();



    NcClose(ncid);
}

}

#endif /* !CMC_NC_WRITER_HXX */
