#ifndef INPUT_CMC_NETCDF_HXX
#define INPUT_CMC_NETCDF_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_coordinate_array.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "input/cmc_input_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"

namespace cmc::input::netcdf
{

constexpr int kDimensionNotConsiderdered = -1;
constexpr int kVariableNotConsiderdered = -1;


class Data
{
public:
    Data() = delete;
    Data(const std::string& path_to_file, const cmc::nc::OpeningMode mode, const MPI_Comm comm = MPI_COMM_WORLD)
    : comm_{comm}
    {
        Open(path_to_file, mode, comm);
    };
    ~Data()
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
        nc::CheckError(err);
        _file_has_been_closed_ = true;
    }

    template<typename... Ts>
    void InquireVariables(const Hyperslab& hyperslab, Ts&&... variable_names);

    void SetHintLongitudeDimension(const int longitude_dimension_id);
    void SetHintLatitudeDimension(const int latitude_dimension_id);
    void SetHintHeightDimension(const int height_dimension_id);
    void SetHintTimeDimension(const int time_dimension_id);

    void SetHintLongitudeDimension(const std::string& longitude_dimension_name);
    void SetHintLatitudeDimension(const std::string& latitude_dimension_name);
    void SetHintHeightDimension(const std::string& height_dimension_name);
    void SetHintTimeDimension(const std::string& time_dimension_name);

    [[nodiscard]] std::vector<input::Var>&& TransferData();

private:
    void Open(const std::string& path_to_file, const nc::OpeningMode mode, const MPI_Comm comm);
    void InquireCoordinates();
    void InquireCoordinateDimensions();
    input::Var SetupVariableData(int, int, std::string&&, DataLayout, DomainIndex, std::vector<Hyperslab>&&, GeoDomain&&, int, const std::array<int, NC_MAX_VAR_DIMS>&);
    input::Var InquireVariable(const Hyperslab&, std::string&&);
    template<typename... Ts> void InquireAllVariables(const Hyperslab&, Ts&&...);

    int ncid_;
    MPI_Comm comm_{MPI_COMM_WORLD};
    
    int num_dimensions_{0};
    int num_global_attributes_{0};
    int id_unlimited_dimension_{-1};

    CoordinateArray<int> coordinate_dimension_ids_{kDimensionNotConsiderdered};
    CoordinateArray<int> coordinate_variable_ids_{kVariableNotConsiderdered};
    CoordinateArray<DomainIndex> coordinate_lengths_;

    std::vector<size_t> dimension_lengths_;
    std::vector<std::string> dimension_names_;

    std::vector<input::Var> variables_;

    bool _file_has_been_closed_{false};
    bool _data_has_been_transfered_{false};
};

template<typename... Ts>
void
Data::InquireAllVariables(const Hyperslab& hyperslab, Ts&&... var_names)
{
    #ifdef CMC_WITH_NETCDF
    (variables_.push_back(InquireVariable(hyperslab, std::forward<Ts>(var_names))), ...);
    #endif
}


template<typename... Ts>
void
Data::InquireVariables(const Hyperslab& hyperslab, Ts&&... variable_names)
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
}

}

#endif /* !INPUT_CMC_NETCDF_HXX */
