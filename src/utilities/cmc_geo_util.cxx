#include "utilities/cmc_geo_util.h"
#include "utilities/cmc_log_functions.h"

////new try

void
cmc_var::switch_data()
{
    if (data != nullptr)
    {
        delete data;
        data = data_new;
        data_new = nullptr;
    }
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
