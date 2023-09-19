#include "utilities/cmc_geo_util.h"
#include "utilities/cmc_log_functions.h"

void
cmc_var::switch_data()
{
    if (data != nullptr)
    {
        delete data;
    }
    data = data_new;
    data_new = nullptr;
}

void
cmc_var::get_data_layout_from_axis_ordering()
{
    /* Only 2D and 3D geo-spatial data variables are supported for the data layout */
    cmc_assert(axis_ordering.size() >= 2);

    if (axis_ordering.size() == 2)
    {
        switch (axis_ordering[0])
        {
            case CMC_COORD_IDS::CMC_LON:
                if (axis_ordering[1] == CMC_COORD_IDS::CMC_LAT)
                {
                    data_layout = DATA_LAYOUT::CMC_2D_LON_LAT;
                } else
                {
                    data_layout = DATA_LAYOUT::CMC_2D_LON_LEV;
                }
            break;
            case CMC_COORD_IDS::CMC_LAT:
                if (axis_ordering[1] == CMC_COORD_IDS::CMC_LON)
                {
                    data_layout = DATA_LAYOUT::CMC_2D_LAT_LON;
                } else
                {
                    data_layout = DATA_LAYOUT::CMC_2D_LAT_LEV;
                }
            break;
            case CMC_COORD_IDS::CMC_LEV:
                if (axis_ordering[1] == CMC_COORD_IDS::CMC_LAT)
                {
                    data_layout = DATA_LAYOUT::CMC_2D_LEV_LAT;
                } else
                {
                    data_layout = DATA_LAYOUT::CMC_2D_LEV_LON;
                }
            break;
            default:
                cmc_err_msg("No geo-spatial coordinate was found in the axis ordering's first position.");
        }
    }
    else if (axis_ordering.size() == 3)
    {
        switch (axis_ordering[0])
        {
            case CMC_COORD_IDS::CMC_LON:
                if (axis_ordering[1] == CMC_COORD_IDS::CMC_LAT)
                {
                    data_layout = DATA_LAYOUT::CMC_3D_LON_LAT_LEV;
                } else {
                    data_layout = DATA_LAYOUT::CMC_3D_LON_LEV_LAT;
                }
            break;
            case CMC_COORD_IDS::CMC_LAT:
                if (axis_ordering[1] == CMC_COORD_IDS::CMC_LON)
                {
                    data_layout = DATA_LAYOUT::CMC_3D_LAT_LON_LEV;
                } else {
                    data_layout = DATA_LAYOUT::CMC_3D_LAT_LEV_LON;
                }
            break;
            case CMC_COORD_IDS::CMC_LEV:
                if (axis_ordering[1] == CMC_COORD_IDS::CMC_LON)
                {
                    data_layout = DATA_LAYOUT::CMC_3D_LEV_LON_LAT;
                } else {
                    data_layout = DATA_LAYOUT::CMC_3D_LEV_LAT_LON;
                }
            break;
            default:
                cmc_err_msg("No geo-spatial coordinate was found in the axis ordering's first position.");
        }
    }
}


std::string
get_coord_name(const CMC_COORD_IDS coord_id)
{
    switch (coord_id)
    {
        case CMC_LAT:
            return std::string{"latitude", 8};
            break;
        case CMC_LON:
            return std::string{"longitude", 9};
            break;
        case CMC_LEV:
            return std::string{"leverage", 8};
            break;
        case CMC_TIME:
            return std::string{"time", 4};
            break;
        default:
            return std::string{};
    }
}

void
cmc_var::assign_an_arbitrary_missing_value()
{
    switch(var_type)
    {
        case cmc_type::CMC_INT8_T:
            missing_value = std::numeric_limits<int8_t>::lowest();
            break;
        case cmc_type::CMC_CHAR:
            missing_value = std::numeric_limits<char>::lowest();
            break;
        case cmc_type::CMC_INT16_T:
            missing_value = std::numeric_limits<int16_t>::lowest();
            break;
        case cmc_type::CMC_INT32_T:
            missing_value = std::numeric_limits<int32_t>::lowest();
            break;
        case cmc_type::CMC_FLOAT:
            missing_value = std::numeric_limits<float>::lowest();
            break;
        case cmc_type::CMC_DOUBLE:
            //Apparently, choosing 'lowest' as missing value, leads to invalid vtk representations in t8code. Therefore, another value is chosen */
            //missing_value = std::numeric_limits<double>::lowest();
            missing_value = static_cast<double>(-32768.0);
            break;
        case cmc_type::CMC_UINT8_T:
            missing_value = std::numeric_limits<uint8_t>::lowest();
            break;
        case cmc_type::CMC_UINT16_T:
            missing_value = std::numeric_limits<uint16_t>::lowest();
            break;
        case cmc_type::CMC_UINT32_T:
            missing_value = std::numeric_limits<uint32_t>::lowest();
            break;
        case cmc_type::CMC_INT64_T:
            missing_value = std::numeric_limits<int64_t>::lowest();
            break;
        case cmc_type::CMC_UINT64_T:
            missing_value = std::numeric_limits<uint64_t>::lowest();
            break;
        case cmc_type::CMC_BYTE:
            missing_value = std::numeric_limits<std::byte>::lowest();
            break;
        default:
            cmc_err_msg("A missing value cannot be supplied for the given cmc type.");
            missing_value = static_cast<int32_t>(0);
    }
}

bool cmc_var::is_equal_to_missing_value(const size_t index) const
{
    if (missing_value_present == false)
    {
        return false;
    }
    cmc_assert(index < data->size());

    /* It is assumed that the missing value is of the same type as the vairbale's data */
    switch (data->get_data_type())
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->get_initial_data_ptr())};
            const int32_t reference_val{data_ptr[index]};
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->get_initial_data_ptr())};
            const float reference_val{data_ptr[index]};
            const float missing_val{std::get<float>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->get_initial_data_ptr())};
            const double reference_val{data_ptr[index]};
            const double missing_val{std::get<double>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->get_initial_data_ptr())};
            const int16_t reference_val{data_ptr[index]};
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->get_initial_data_ptr())};
            const int64_t reference_val{data_ptr[index]};
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->get_initial_data_ptr())};
            const uint64_t reference_val{data_ptr[index]};
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->get_initial_data_ptr())};
            const uint32_t reference_val{data_ptr[index]};
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->get_initial_data_ptr())};
            const int8_t reference_val{data_ptr[index]};
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->get_initial_data_ptr())};
            const uint8_t reference_val{data_ptr[index]};
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->get_initial_data_ptr())};
            const uint16_t reference_val{data_ptr[index]};
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        case CMC_BYTE:
        {
            std::byte* data_ptr{static_cast<std::byte*>(data->get_initial_data_ptr())};
            const std::byte reference_val{data_ptr[index]};
            const std::byte missing_val{std::get<std::byte>(missing_value)};
            return (reference_val == missing_val);
        }
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->get_initial_data_ptr())};
            const char reference_val{data_ptr[index]};
            const char missing_val{std::get<char>(missing_value)};
            return cmc_approx_cmp(reference_val, missing_val);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
            return false;
    }
}

bool
cmc_value_equal_to_missing_value(const cmc_universal_type_t& value, const cmc_universal_type_t& missing_value)
{
    //TODO: not implemented yet
    return false;
}

bool
cmc_value_equal_to_zero(const cmc_universal_type_t& value)
{
    /* Get the current type of the universal_type */
    switch (static_cast<cmc_type>(value.index()))
    {
        case CMC_INT32_T:
            if(std::get<int32_t>(value) == static_cast<int32_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_FLOAT:
            if(std::get<float>(value) == static_cast<float>(0.0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_DOUBLE:
            if(std::get<double>(value) == static_cast<double>(0.0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_INT16_T:
            if(std::get<int16_t>(value) == static_cast<int16_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_INT64_T:
            if(std::get<int64_t>(value) == static_cast<int64_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_UINT64_T:
            if(std::get<uint64_t>(value) == static_cast<uint64_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_UINT32_T:
            if(std::get<uint32_t>(value) == static_cast<uint32_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_INT8_T:
            if(std::get<int8_t>(value) == static_cast<int8_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_UINT8_T:
            if(std::get<uint8_t>(value) == static_cast<uint8_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_UINT16_T:
            if(std::get<uint16_t>(value) == static_cast<uint16_t>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_BYTE:
            if(std::get<std::byte>(value) == static_cast<std::byte>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        case CMC_CHAR:
            if(std::get<char>(value) == static_cast<char>(0))
            {
                return true;
            } else
            {
                return false;
            }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

/**
 * @brief Calculate an upper bound for a coarsening step. This functions estimates the deviation from the potential next coarsening @var current_mean 
 *        to the initial data given the parameters below. 
 * 
 * @param previous_max_deviation The previous maximum deviation from the former coarsening step
 * @param previous_value The value of one of the elements whose value will be (eventually) repalced by the next coarsening by the value @var current_mean 
 * @param new_value The calculated value which will replace @var previous_mean and it's siblings' values
 * @return double The estimate of the relative deviation 
 */
double
calculate_two_step_relative_max_deviation(const double previous_max_deviation, const cmc_universal_type_t& previous_value, const cmc_universal_type_t& new_value)
{
    cmc_assert(previous_value.index() == new_value.index());
    cmc_assert(previous_max_deviation >= 0.0);

    /* Convert the new value, which is ought be a replacement for the previous element value, to double */
    const double dnew_value = std::visit([](auto& val) -> double {return static_cast<double>(val);}, new_value);

    /* Convert the previous element value to double */
    double dprevious_value = std::abs(std::visit([](auto& val) -> double {return static_cast<double>(val);}, previous_value));

    /* Check if the new_value is unequal to zero */
    if (!cmc_approx_cmp(dprevious_value, static_cast<double>(0.0)))
    {
        return std::max(std::abs(dprevious_value + previous_max_deviation - dnew_value) / std::abs(dprevious_value + previous_max_deviation),
                        std::abs(dprevious_value - previous_max_deviation - dnew_value) / std::abs(dprevious_value - previous_max_deviation));
    } else
    {
        /* Check if the previous value is zero */
        if (cmc_approx_cmp(dnew_value, static_cast<double>(0.0)))
        {
            /* If so, there is no deviation */
            return static_cast<double>(0.0);
        } else
        {
            /* If the new value is unequal to zero, we use a measure of relative difference */
            return static_cast<double>(2.0);
        }
    }
}

double
calculate_two_step_absolute_max_deviation(const double previous_max_deviation, const cmc_universal_type_t& previous_value, const cmc_universal_type_t& new_value)
{
    cmc_assert(previous_value.index() == new_value.index());
    cmc_assert(previous_max_deviation >= 0.0);

    /* Convert the new value, which is ought be a replacement for the previous element value, to double */
    const double dnew_value = std::visit([](auto& val) -> double {return static_cast<double>(val);}, new_value);

    /* Convert the previous element value to double */
    double dprevious_value = std::visit([](auto& val) -> double {return static_cast<double>(val);}, previous_value);

    return std::abs(dnew_value - dprevious_value) + previous_max_deviation;
}

CMC_COORD_IDS
cmc_get_split_dim_id(const DATA_LAYOUT initial_layout, const DATA_LAYOUT preferred_layout)
{
    cmc_assert(initial_layout > preferred_layout);
    
    switch(static_cast<int>(DATA_LAYOUT::_INTERN_ID_END_2D_START_3D) - static_cast<int>(preferred_layout))
    {
        case 1:
            [[fallthrough]];
        case 2:
            return CMC_COORD_IDS::CMC_LAT;
        break;
        case 3:
            [[fallthrough]];
        case 4:
            return CMC_COORD_IDS::CMC_LON;
        break;
        case 5:
            [[fallthrough]];
        case 6:
            return CMC_COORD_IDS::CMC_LEV;
        break;
        default:
            cmc_err_msg("The split dimension spuld not be determined.");
            return CMC_LAT; //This return value has to be set only because there is no error value at the moment
    }
}

/**
 * @brief This functions returns the axis ordering given a certain @var layout
 * 
 * @param layout The DATA_LAYOUT which should be transformed to a vector
 * @return std::vector<int> A vector containing the @enum's CMC_COORD_IDS ordered as given by the @var layout
 */
std::vector<int>
cmc_get_axis_ordering_from_layout(const DATA_LAYOUT layout)
{
    switch(layout)
    {
        case DATA_LAYOUT::CMC_2D_LAT_LON:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LAT), static_cast<int>(CMC_COORD_IDS::CMC_LON)};
        break;
        case DATA_LAYOUT::CMC_2D_LON_LAT:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LON), static_cast<int>(CMC_COORD_IDS::CMC_LAT)};
        break;
        case DATA_LAYOUT::CMC_2D_LAT_LEV:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LAT), static_cast<int>(CMC_COORD_IDS::CMC_LEV)};
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LAT:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LEV), static_cast<int>(CMC_COORD_IDS::CMC_LAT)};
        break;
        case DATA_LAYOUT::CMC_2D_LON_LEV:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LON), static_cast<int>(CMC_COORD_IDS::CMC_LEV)};
        break;
        case DATA_LAYOUT::CMC_2D_LEV_LON:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LEV), static_cast<int>(CMC_COORD_IDS::CMC_LON)};
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LON_LEV:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LAT), static_cast<int>(CMC_COORD_IDS::CMC_LON), static_cast<int>(CMC_COORD_IDS::CMC_LEV)};
        break;
        case DATA_LAYOUT::CMC_3D_LAT_LEV_LON:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LAT), static_cast<int>(CMC_COORD_IDS::CMC_LEV), static_cast<int>(CMC_COORD_IDS::CMC_LON)};
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LAT_LON:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LEV), static_cast<int>(CMC_COORD_IDS::CMC_LAT), static_cast<int>(CMC_COORD_IDS::CMC_LON)};
        break;
        case DATA_LAYOUT::CMC_3D_LEV_LON_LAT:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LEV), static_cast<int>(CMC_COORD_IDS::CMC_LON), static_cast<int>(CMC_COORD_IDS::CMC_LAT)};
        break;
        case DATA_LAYOUT::CMC_3D_LON_LEV_LAT:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LON), static_cast<int>(CMC_COORD_IDS::CMC_LEV), static_cast<int>(CMC_COORD_IDS::CMC_LAT)};
        break;
        case DATA_LAYOUT::CMC_3D_LON_LAT_LEV:
            return std::vector<int>{static_cast<int>(CMC_COORD_IDS::CMC_LAT), static_cast<int>(CMC_COORD_IDS::CMC_LON), static_cast<int>(CMC_COORD_IDS::CMC_LEV)};
        break;
        default:
            cmc_err_msg("An unknown data layout has been supplied.");
            return std::vector<int>{static_cast<int>(CMC_ERR)};
    }
}

inline
uint64_t
cmc_coordinate::cmc_get_cart_coord_by_id(const cmc_coordinate_id _id) const
{
    switch(_id)
    {
        case cmc_coordinate_id::_CMC_LON:
            return std::get<0>(cartesian_coordinate);
        break;
        case cmc_coordinate_id::_CMC_LAT:
            return std::get<1>(cartesian_coordinate);
        break;
        case cmc_coordinate_id::_CMC_LEV:
            return std::get<2>(cartesian_coordinate);
        break;
        default:
            cmc_err_msg("Unknown coordinate ID.");
            return 0;
    }
}

inline
uint64_t
cmc_coordinate::cmc_get_ref_coord(const cmc_coordinate_id _coord_id) const
{
    switch(repr)
    {
        case current_representation::_CMC_CART_COORD:
            return cmc_get_cart_coord_by_id(_coord_id);
        break;
        case current_representation::_CMC_MORTON_IDX:
            cmc_err_msg("not yet implemented");
            return 0;
        break;
        default:
            cmc_err_msg("Cannot obtain a reference coordinate from the current representation.");
    }
}

uint64_t
cmc_coordinate::get_longitude_ref_coordinate() const
{
    return cmc_get_ref_coord(cmc_coordinate_id::_CMC_LON);
}

uint64_t
cmc_coordinate::get_latitude_ref_coordinate() const
{
    return cmc_get_ref_coord(cmc_coordinate_id::_CMC_LAT);
}

uint64_t
cmc_coordinate::get_elevation_ref_coordinate() const
{
    return cmc_get_ref_coord(cmc_coordinate_id::_CMC_LEV);
}

void
cmc_ref_coordinates::ref()
{
    /* Increment the reference count */
    ++reference_count;
}

void
cmc_ref_coordinates::unref()
{
    /* Decrement the reference count */
    if(--reference_count < 0)
    {
        //Nothing to be deallocated by now
    }
}

cmc_universal_type_t
cmc_global_coordinate_system::get_coordinate_value_at_dim(const CMC_COORD_IDS coord_id, const size_t reference_position)
{

    cmc_universal_type_t ret_val{static_cast<int>(0)};

    switch (coord_id)
    {
        case CMC_COORD_IDS::CMC_LON:
            cmc_assert(reference_position < ((coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LON))).size()));
            ret_val = (coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LON))).operator[](reference_position);
        break;
        case CMC_COORD_IDS::CMC_LAT:
            cmc_assert(reference_position < (coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LAT)).size()));
            ret_val = (coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LAT))).operator[](reference_position);
        break;
        case CMC_COORD_IDS::CMC_LEV:
            cmc_assert(reference_position < (coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LEV)).size()));
            ret_val = (coords.operator[](static_cast<size_t>(CMC_COORD_IDS::CMC_LEV))).operator[](reference_position);
        break;
        case CMC_COORD_IDS::CMC_TIME:
            cmc_err_msg("Currently, time series data is not supported.");
        break;
        default:
            cmc_err_msg("An unknown coordinate dimension was supplied.");
    }

    return ret_val;
}

cmc_universal_type_t
cmc_global_coordinate_system::get_coordinate_value_at_longitude(const size_t reference_position)
{
    /* Return the longitude value at the reference position */
    return get_coordinate_value_at_dim(CMC_COORD_IDS::CMC_LON, reference_position);
}

cmc_universal_type_t
cmc_global_coordinate_system::get_coordinate_value_at_latitude(const size_t reference_position)
{
    /* Return the latitude value at the reference position */
    return get_coordinate_value_at_dim(CMC_COORD_IDS::CMC_LAT, reference_position);
}

cmc_universal_type_t
cmc_global_coordinate_system::get_coordinate_value_at_elevation(const size_t reference_position)
{
    /* Return the elevation value at the reference position */
    return get_coordinate_value_at_dim(CMC_COORD_IDS::CMC_LEV, reference_position);
}

void
cmc_global_coordinate_system::ref()
{
    /* Increase the reference count */
    ++(reference_count);
}

void
cmc_global_coordinate_system::unref()
{
    /* Decrease the reference count */
    if((--reference_count) < 0)
    {
        /* Deallocate all global coordinate data */
        //delete coords;
        //coords = nullptr;
    }
}

void
cmc_global_coordinate_system::reset_reference_count()
{
    /* Reset the reference count to zero */
    reference_count = 0;
}
