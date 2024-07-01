#ifndef CMC_NETCDF_HXX
#define CMC_NETCDF_HXX
/**
 * @file cmc_netcdf.hxx
 * @brief Via the 'include' of @file cmc_netcdf.h, this file collects all functions used for accessing netCDF files and storing the data of the netCDF variables as well as the geo-spatial domain on which the variables are defined.
 * Additionally, this file supplies C++ only functions for inquiring variables from a netCDF file
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_coordinate_array.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_input_variable.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf.h"
#endif
#ifdef CMC_WITH_NETCDF_PAR
#include "netcdf_par.h"
#endif

#include <vector>
#include <string>

namespace cmc
{

constexpr int kDimensionNotConsiderdered = -1;
constexpr int kVariableNotConsiderdered = -1;

/* Enum describing an access mode for netCDF files */
enum NcOpeningMode {Serial = 0, Parallel};

[[noreturn]] void
NcExit(const int _err_code, const char* _location);

#define NcCheckError(err) ((err) == NC_NOERR ? (void) 0 : NcExit(err, CMC_FILE_LOCATION))

class NcData
{
public:
    NcData() = delete;
    NcData(const std::string& path_to_file, const NcOpeningMode mode, const MPI_Comm comm = MPI_COMM_WORLD)
    {
        NcOpen(path_to_file, mode, comm);
    };
    ~NcData()
    {
        if (!_file_has_been_closed_)
        {
            /* Close the still open file */
            CloseFileHandle();
        }
    };

    void CloseFileHandle()
    {
        int err = nc_close(ncid_);
        NcCheckError(err);
        _file_has_been_closed_ = true;
    }

    template<typename... Ts>
    void InquireVariables(const Hyperslab& hyperslab, Ts&&... variable_names);

    void SetHintLongitudeDimension(const int longitude_dimension_id);

    void SetHintLatitudeDimension(const int latitude_dimension_id);

    void SetHintHeightDimension(const int height_dimension_id);

    void SetHintTimeDimension(const int time_dimension_id);

    [[nodiscard]] std::vector<InputVar>&& TransferData();

private:
    void NcOpen(const std::string& path_to_file, const NcOpeningMode mode, const MPI_Comm comm);
    void InquireCoordinates();
    void InquireCoordinateDimensions();
    void InquireVariableData();
    InputVar SetupVariableData(int, int, std::string&&, DataLayout, DomainIndex, std::vector<Hyperslab>&&, GeoDomain&&, int, const std::array<int, NC_MAX_VAR_DIMS>&);
    InputVar InquireVariable(const Hyperslab&, std::string&&);
    template<typename... Ts> void InquireAllVariables(const Hyperslab&, Ts&&...);

    int ncid_;

    int num_dimensions_{0};
    int num_global_attributes_{0};
    int id_unlimited_dimension_{-1};

    CoordinateArray<int> coordinate_dimension_ids_{kDimensionNotConsiderdered};
    CoordinateArray<int> coordinate_variable_ids_{kVariableNotConsiderdered};
    CoordinateArray<DomainIndex> coordinate_lengths_;

    std::vector<size_t> dimension_lengths_;
    std::vector<std::string> dimension_names_;

    std::vector<InputVar> variables_;

    bool _file_has_been_closed_{false};
    bool _data_has_been_transfered_{false};
};

template<typename... Ts>
void
NcData::InquireAllVariables(const Hyperslab& hyperslab, Ts&&... var_names)
{
    #ifdef CMC_WITH_NETCDF
    (variables_.push_back(InquireVariable(hyperslab, std::forward<Ts>(var_names))), ...);
    #endif
};


template<typename... Ts>
void
NcData::InquireVariables(const Hyperslab& hyperslab, Ts&&... variable_names)
{
    #ifdef CMC_WITH_NETCDF
    /* Inquire information about dimensions and read coordinate variables */
    InquireCoordinates();

    /* Reserve memory for the variables */
    variables_.reserve(sizeof...(Ts));

    InquireAllVariables(hyperslab, std::forward<Ts>(variable_names)...);

    #else
    cmc_err_msg("CMC is not compiled with netCDF, please reconfigure with netCDF linkage in order to use this function.");
    #endif
};

int NcOpenSerial(const char* path_to_file);

int
NcOpenParallel(const char* path_to_file, MPI_Comm comm);


}

#endif /* CMC_NETCDF_HXX */
