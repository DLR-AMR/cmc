#include "utilities/cmc_container.h"

struct var_array
{
private:
public:
    var_array(const size_t _num_elems, const cmc_type _type)
    : num_elements{_num_elems}, type{_type}{
        initial_data_ptr = malloc(_num_elems * cmc_type_to_bytes[_type]);
    };
    ~var_array(){
        if (initial_data_ptr != nullptr)
        {
            free(initial_data_ptr);
            initial_data_ptr = nullptr;
        }
    };

    const size_t num_elements;
    const cmc_type type;
    void* initial_data_ptr{nullptr};
};

void
var_array_t::_init_array(const size_t num_elements, const cmc_type _type)
{
    data = new var_array(num_elements, _type);
}

void
var_array_t::_destroy_array()
{
    if (data != nullptr)
    {
        delete data;
        data = nullptr;
    }
}

var_array_t::var_array_t(const var_array_t& array)
{
    data = new var_array(array.data->num_elements, array.data->type);
    if (array.data->initial_data_ptr != nullptr)
    {
        memcpy(data->initial_data_ptr, array.data->initial_data_ptr, array.data->num_elements * cmc_type_to_bytes[array.data->type]);
    }
}

var_array_t&
var_array_t::operator=(const var_array_t& array)
{
    if (this == &array)
    {
        return *this;
    }
    if (data != nullptr)
    {
        delete data;
        data = nullptr;
    }
    data = new var_array(array.data->num_elements, array.data->type);
    if (array.data->initial_data_ptr != nullptr)
    {
        memcpy(data->initial_data_ptr, array.data->initial_data_ptr, array.data->num_elements * cmc_type_to_bytes[array.data->type]);
    }
    return (*this);
}


#if 0
switch (data->type)
{
    case CMC_INT32_T:

    break;
    case CMC_FLOAT:

    break;
    case CMC_DOUBLE:

    break;
    case CMC_INT16_T:

    break;
    case CMC_INT64_T:

    break;
    case CMC_UINT64_T:

    break;
    case CMC_UINT32_T:

    break;
    case CMC_INT8_T:

    break;
    case CMC_UINT8_T:

    break;
    case CMC_UINT16_T:

    break;
    case CMC_BYTE:

    break;
    case CMC_CHAR:

    break;
    default:
        cmc_err_msg("An unknown cmc data type has been supplied.");
}
#endif

//TODO: complete these functions
template<typename T>
static auto
get_cmc_type()
    -> std::enable_if_t<std::is_same_v<T, int32_t>, cmc_type>
{
    return cmc_type::CMC_INT32_T;
}
template<typename T>
static auto
get_cmc_type()
    -> std::enable_if_t<std::is_same_v<T, double>, cmc_type>
{
    return cmc_type::CMC_DOUBLE;
}
template<typename T>
static auto
get_cmc_type()
    -> std::enable_if_t<std::is_same_v<T, float>, cmc_type>
{
    return cmc_type::CMC_FLOAT;
}
template<typename T>
static auto
get_cmc_type()
    -> std::enable_if_t<std::is_same_v<T, int16_t>, cmc_type>
{
    return cmc_type::CMC_INT16_T;
}

void
var_array_t::axpy_scalar_with_missing_vals(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& add_offset, const cmc_universal_type_t& missing_value)
{
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            if (const int32_t* missing_val = std::get_if<int32_t>(&missing_value))
            {
                if(std::holds_alternative<int32_t>(scale_factor) && std::holds_alternative<int32_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    int32_t* dest_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                    const int32_t scaling{std::get<int32_t>(scale_factor)};
                    const int32_t offset{std::get<int32_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int32_t* src_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int32_t missing value present.");
            }
        }
        break;
        case CMC_FLOAT:
        {
            if (const float* missing_val = std::get_if<float>(&missing_value))
            {
                if(std::holds_alternative<float>(scale_factor) && std::holds_alternative<float>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    float* dest_ptr{static_cast<float*>(data->initial_data_ptr)};
                    const float scaling{std::get<float>(scale_factor)};
                    const float offset{std::get<float>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    float* src_ptr{static_cast<float*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no float missing value present.");
            }
        }
        break;
        case CMC_DOUBLE:
        {
            if (const double* missing_val = std::get_if<double>(&missing_value))
            {
                if(std::holds_alternative<double>(scale_factor) && std::holds_alternative<double>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    double* dest_ptr{static_cast<double*>(data->initial_data_ptr)};
                    const double scaling{std::get<double>(scale_factor)};
                    const double offset{std::get<double>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    double* src_ptr{static_cast<double*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no double missing value present.");
            }
        }
        break;
        case CMC_INT16_T:
        {
            if (const int16_t* missing_val = std::get_if<int16_t>(&missing_value))
            {
                if(std::holds_alternative<int16_t>(scale_factor) && std::holds_alternative<int16_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    int16_t* dest_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                    const int16_t scaling{std::get<int16_t>(scale_factor)};
                    const int16_t offset{std::get<int16_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int16_t* src_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int16_t missing value present.");
            }
        }
        break;
        case CMC_INT64_T:
        {
            if (const int64_t* missing_val = std::get_if<int64_t>(&missing_value))
            {
                if(std::holds_alternative<int64_t>(scale_factor) && std::holds_alternative<int64_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    int64_t* dest_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                    const int64_t scaling{std::get<int64_t>(scale_factor)};
                    const int64_t offset{std::get<int64_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int64_t* src_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int64_t missing value present.");
            }
        }
        break;
        case CMC_UINT64_T:
        {
            if (const uint64_t* missing_val = std::get_if<uint64_t>(&missing_value))
            {
                if(std::holds_alternative<uint64_t>(scale_factor) && std::holds_alternative<uint64_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    uint64_t* dest_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                    const uint64_t scaling{std::get<uint64_t>(scale_factor)};
                    const uint64_t offset{std::get<uint64_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint64_t* src_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint64_t missing value present.");
            }
        }
        break;
        case CMC_UINT32_T:
        {
            if (const uint32_t* missing_val = std::get_if<uint32_t>(&missing_value))
            {
                if(std::holds_alternative<uint32_t>(scale_factor) && std::holds_alternative<uint32_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    uint32_t* dest_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                    const uint32_t scaling{std::get<uint32_t>(scale_factor)};
                    const uint32_t offset{std::get<uint32_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint32_t* src_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint32_t missing value present.");
            }
        }
        break;
        case CMC_INT8_T:
        {
            if (const int8_t* missing_val = std::get_if<int8_t>(&missing_value))
            {
                if(std::holds_alternative<int8_t>(scale_factor) && std::holds_alternative<int8_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    int8_t* dest_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                    const int8_t scaling{std::get<int8_t>(scale_factor)};
                    const int8_t offset{std::get<int8_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int8_t* src_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int8_t missing value present.");
            }
        }
        break;
        case CMC_UINT8_T:
        {
            if (const uint8_t* missing_val = std::get_if<uint8_t>(&missing_value))
            {
                if(std::holds_alternative<uint8_t>(scale_factor) && std::holds_alternative<uint8_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    uint8_t* dest_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                    const uint8_t scaling{std::get<uint8_t>(scale_factor)};
                    const uint8_t offset{std::get<uint8_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint8_t* src_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint8_t missing value present.");
            }
        }
        break;
        case CMC_UINT16_T:
        {
            if (const uint16_t* missing_val = std::get_if<uint16_t>(&missing_value))
            {
                if(std::holds_alternative<uint16_t>(scale_factor) && std::holds_alternative<uint16_t>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    uint16_t* dest_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                    const uint16_t scaling{std::get<uint16_t>(scale_factor)};
                    const uint16_t offset{std::get<uint16_t>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint16_t* src_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint16_t missing value present.");
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add to data to type byte.");
        break;
        case CMC_CHAR:
        {
            if (const char* missing_val = std::get_if<char>(&missing_value))
            {
                if(std::holds_alternative<char>(scale_factor) && std::holds_alternative<char>(add_offset))
                {
                    /* If the scale factor is of the same type as the data */
                    char* dest_ptr{static_cast<char*>(data->initial_data_ptr)};
                    const char scaling{std::get<char>(scale_factor)};
                    const char offset{std::get<char>(add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * dest_ptr[i] + offset;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    char* src_ptr{static_cast<char*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no char missing value present.");
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void
var_array_t::add_const_with_missing_vals(const cmc_universal_type_t& offset, const cmc_universal_type_t& missing_value)
{
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            if (const int32_t* missing_val = std::get_if<int32_t>(&missing_value))
            {
                if(std::holds_alternative<int32_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    int32_t* dest_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                    const int32_t offset_val{std::get<int32_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int32_t* src_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int32_t missing value present.");
            }
        }
        break;
        case CMC_FLOAT:
        {
            if (const float* missing_val = std::get_if<float>(&missing_value))
            {
                if(std::holds_alternative<float>(offset))
                {
                    /* If the offset is of the same type as the data */
                    float* dest_ptr{static_cast<float*>(data->initial_data_ptr)};
                    const float offset_val{std::get<float>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    float* src_ptr{static_cast<float*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no float missing value present.");
            }
        }
        break;
        case CMC_DOUBLE:
        {
            if (const double* missing_val = std::get_if<double>(&missing_value))
            {
                if(std::holds_alternative<double>(offset))
                {
                    /* If the offset is of the same type as the data */
                    double* dest_ptr{static_cast<double*>(data->initial_data_ptr)};
                    const double offset_val{std::get<double>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    double* src_ptr{static_cast<double*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no double missing value present.");
            }
        }
        break;
        case CMC_INT16_T:
        {
            if (const int16_t* missing_val = std::get_if<int16_t>(&missing_value))
            {
                if(std::holds_alternative<int16_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    int16_t* dest_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                    const int16_t offset_val{std::get<int16_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int16_t* src_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int16_t missing value present.");
            }
        }
        break;
        case CMC_INT64_T:
        {
            if (const int64_t* missing_val = std::get_if<int64_t>(&missing_value))
            {
                if(std::holds_alternative<int64_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    int64_t* dest_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                    const int64_t offset_val{std::get<int64_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int64_t* src_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int64_t missing value present.");
            }
        }
        break;
        case CMC_UINT64_T:
        {
            if (const uint64_t* missing_val = std::get_if<uint64_t>(&missing_value))
            {
                if(std::holds_alternative<uint64_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    uint64_t* dest_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                    const uint64_t offset_val{std::get<uint64_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint64_t* src_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint64_t missing value present.");
            }
        }
        break;
        case CMC_UINT32_T:
        {
            if (const uint32_t* missing_val = std::get_if<uint32_t>(&missing_value))
            {
                if(std::holds_alternative<uint32_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    uint32_t* dest_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                    const uint32_t offset_val{std::get<uint32_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint32_t* src_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint32_t missing value present.");
            }
        }
        break;
        case CMC_INT8_T:
        {
            if (const int8_t* missing_val = std::get_if<int8_t>(&missing_value))
            {
                if(std::holds_alternative<int8_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    int8_t* dest_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                    const int8_t offset_val{std::get<int8_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int8_t* src_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int8_t missing value present.");
            }
        }
        break;
        case CMC_UINT8_T:
        {
            if (const uint8_t* missing_val = std::get_if<uint8_t>(&missing_value))
            {
                if(std::holds_alternative<uint8_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    uint8_t* dest_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                    const uint8_t offset_val{std::get<uint8_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint8_t* src_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint8_t missing value present.");
            }
        }
        break;
        case CMC_UINT16_T:
        {
            if (const uint16_t* missing_val = std::get_if<uint16_t>(&missing_value))
            {
                if(std::holds_alternative<uint16_t>(offset))
                {
                    /* If the offset is of the same type as the data */
                    uint16_t* dest_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                    const uint16_t offset_val{std::get<uint16_t>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint16_t* src_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint16_t missing value present.");
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add to dat of type byte.");
        break;
        case CMC_CHAR:
        {
            if (const char* missing_val = std::get_if<char>(&missing_value))
            {
                if(std::holds_alternative<char>(offset))
                {
                    /* If the offset is of the same type as the data */
                    char* dest_ptr{static_cast<char*>(data->initial_data_ptr)};
                    const char offset_val{std::get<char>(offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] += offset_val;
                        }
                    }
                } else
                {
                    /* If the offset and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array(data->num_elements, get_cmc_type<cmc_standard_type>());
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    char* src_ptr{static_cast<char*>(data->initial_data_ptr)};
                    const cmc_standard_type offset_val{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, offset)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(offset_val) + static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the offset data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no char missing value present.");
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

///////hier gehts weiter; fuer alle faelle und noch fuer Funtkionen (add_const, axpy_scalar) fertig machen
void
var_array_t::scale_with_missing_vals(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& missing_value)
{
    switch(data->type)
    {
        case CMC_INT32_T:
        {
            if (const int32_t* missing_val = std::get_if<int32_t>(&missing_value))
            {
                if(std::holds_alternative<int32_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    int32_t* dest_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                    const int32_t scaling{std::get<int32_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int32_t* src_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int32_t missing value present.");
            }
        }
        break;
        case CMC_FLOAT:
        {
            if (const float* missing_val = std::get_if<float>(&missing_value))
            {
                if(std::holds_alternative<float>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    float* dest_ptr{static_cast<float*>(data->initial_data_ptr)};
                    const float scaling{std::get<float>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    float* src_ptr{static_cast<float*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no float missing value present.");
            }
        }
        break;
        case CMC_DOUBLE:
        {
            if (const double* missing_val = std::get_if<double>(&missing_value))
            {
                if(std::holds_alternative<double>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    double* dest_ptr{static_cast<double*>(data->initial_data_ptr)};
                    const double scaling{std::get<double>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    double* src_ptr{static_cast<double*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no double missing value present.");
            }
        }
        break;
        case CMC_INT16_T:
        {
            if (const int16_t* missing_val = std::get_if<int16_t>(&missing_value))
            {
                if(std::holds_alternative<int16_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    int16_t* dest_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                    const int16_t scaling{std::get<int16_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int16_t* src_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int16_t missing value present.");
            }
        }
        break;
        case CMC_INT64_T:
        {
            if (const int64_t* missing_val = std::get_if<int64_t>(&missing_value))
            {
                if(std::holds_alternative<int64_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    int64_t* dest_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                    const int64_t scaling{std::get<int64_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int64_t* src_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int64_t missing value present.");
            }
        }
        break;
        case CMC_UINT64_T:
        {
            if (const uint64_t* missing_val = std::get_if<uint64_t>(&missing_value))
            {
                if(std::holds_alternative<uint64_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    uint64_t* dest_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                    const uint64_t scaling{std::get<uint64_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint64_t* src_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint64_t missing value present.");
            }
        }
        break;
        case CMC_UINT32_T:
        {
            if (const uint32_t* missing_val = std::get_if<uint32_t>(&missing_value))
            {
                if(std::holds_alternative<uint32_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    uint32_t* dest_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                    const uint32_t scaling{std::get<uint32_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint32_t* src_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint32_t missing value present.");
            }
        }
        break;
        case CMC_INT8_T:
        {
            if (const int8_t* missing_val = std::get_if<int8_t>(&missing_value))
            {
                if(std::holds_alternative<int8_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    int8_t* dest_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                    const int8_t scaling{std::get<int8_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    int8_t* src_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no int8_t missing value present.");
            }
        }
        break;
        case CMC_UINT8_T:
        {
            if (const uint8_t* missing_val = std::get_if<uint8_t>(&missing_value))
            {
                if(std::holds_alternative<uint8_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    uint8_t* dest_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                    const uint8_t scaling{std::get<uint8_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint8_t* src_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint8_t missing value present.");
            }
        }
        break;
        case CMC_UINT16_T:
        {
            if (const uint16_t* missing_val = std::get_if<uint16_t>(&missing_value))
            {
                if(std::holds_alternative<uint16_t>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    uint16_t* dest_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                    const uint16_t scaling{std::get<uint16_t>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    uint16_t* src_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no uint16_t missing value present.");
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot apply scaling to data of type byte.");
        break;
        case CMC_CHAR:
        {
            if (const char* missing_val = std::get_if<char>(&missing_value))
            {
                if(std::holds_alternative<char>(scale_factor))
                {
                    /* If the scale factor is of the same type as the data */
                    char* dest_ptr{static_cast<char*>(data->initial_data_ptr)};
                    const char scaling{std::get<char>(scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(dest_ptr[i], *missing_val))
                        {
                            dest_ptr[i] *= scaling;
                        }
                    }
                } else
                {
                    /* If the scale factor and the data are NOT of the same data type */
                    /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                    var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                    cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                    char* src_ptr{static_cast<char*>(data->initial_data_ptr)};
                    const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                    for(size_t i{0}; i < data->num_elements; ++i)
                    {
                        /* If the data array holds a missing value, it will be skipped */
                        if (!cmc_approx_cmp(src_ptr[i], *missing_val))
                        {
                            dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                        }
                    }
                    /* Deallocate the old data array */
                    delete data;
                    /* Assign the scaled data array */
                    data = data_new;
                }
            } else
            {
                cmc_err_msg("There is no char missing value present.");
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}


void
var_array_t::add(var_array_t& array)
{
    assert(array.data->type == data->type);
    assert(data->num_elements == array.data->num_elements);

    switch(data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* dest_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            int32_t* src_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_FLOAT:
        {
            float* dest_ptr{static_cast<float*>(data->initial_data_ptr)};
            float* src_ptr{static_cast<float*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_DOUBLE:
        {
            double* dest_ptr{static_cast<double*>(data->initial_data_ptr)};
            double* src_ptr{static_cast<double*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* dest_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            int16_t* src_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* dest_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            int64_t* src_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* dest_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            uint64_t* src_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* dest_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            uint32_t* src_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* dest_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            int8_t* src_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* dest_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            uint8_t* src_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* dest_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            uint16_t* src_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* dest_ptr{static_cast<char*>(data->initial_data_ptr)};
            char* src_ptr{static_cast<char*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                dest_ptr[i] += src_ptr[i];
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

size_t
var_array_t::size() const
{
    return data->num_elements;
}

cmc_universal_type_t
var_array_t::sum_over_range(const size_t start_index, const size_t end_index) const
{
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            int32_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            float sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            double sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            int16_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            int64_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            uint64_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            uint32_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            int8_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            uint8_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            uint16_t sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            char sum{0};
            for(size_t i{start_index}; i <= end_index; ++i)
            {
                sum += data_ptr[i];
            }
            return cmc_universal_type_t{sum};
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
    return cmc_universal_type_t{0};
}

void
var_array_t::scale(const cmc_universal_type_t& scale_factor)
{
    switch(data->type)
    {
        case CMC_INT32_T:
        {
            if(std::holds_alternative<int32_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                int32_t* dest_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                const int32_t scaling{std::get<int32_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int32_t* src_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_FLOAT:
        {
            if(std::holds_alternative<float>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                float* dest_ptr{static_cast<float*>(data->initial_data_ptr)};
                const float scaling{std::get<float>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                float* src_ptr{static_cast<float*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_DOUBLE:
        {
            if(std::holds_alternative<double>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                double* dest_ptr{static_cast<double*>(data->initial_data_ptr)};
                const double scaling{std::get<double>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                double* src_ptr{static_cast<double*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT16_T:
        {
            if(std::holds_alternative<int16_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                int16_t* dest_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                const int16_t scaling{std::get<int16_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int16_t* src_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT64_T:
        {
            if(std::holds_alternative<int64_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                int64_t* dest_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                const int64_t scaling{std::get<int64_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int64_t* src_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT64_T:
        {
            if(std::holds_alternative<uint64_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                uint64_t* dest_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                const uint64_t scaling{std::get<uint64_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint64_t* src_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT32_T:
        {
            if(std::holds_alternative<uint32_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                uint32_t* dest_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                const uint32_t scaling{std::get<uint32_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint32_t* src_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT8_T:
        {
            if(std::holds_alternative<int8_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                int8_t* dest_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                const int8_t scaling{std::get<int8_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int8_t* src_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT8_T:
        {
            if(std::holds_alternative<uint8_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                uint8_t* dest_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                const uint8_t scaling{std::get<uint8_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint8_t* src_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT16_T:
        {
            if(std::holds_alternative<uint16_t>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                uint16_t* dest_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                const uint16_t scaling{std::get<uint16_t>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint16_t* src_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot apply scaling to data of type byte.");
        break;
        case CMC_CHAR:
        {
            if(std::holds_alternative<char>(scale_factor))
            {
                /* If the scale factor is of the same type as the data */
                char* dest_ptr{static_cast<char*>(data->initial_data_ptr)};
                const char scaling{std::get<char>(scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] *= scaling;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                char* src_ptr{static_cast<char*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(scaling) * static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void
var_array_t::add_const(const cmc_universal_type_t& constant_summand)
{
    switch(data->type)
    {
        case CMC_INT32_T:
        {
            if(std::holds_alternative<int32_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                int32_t* dest_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                const int32_t offset{std::get<int32_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int32_t* src_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_FLOAT:
        {
            if(std::holds_alternative<float>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                float* dest_ptr{static_cast<float*>(data->initial_data_ptr)};
                const float offset{std::get<float>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                float* src_ptr{static_cast<float*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_DOUBLE:
        {
            if(std::holds_alternative<double>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                double* dest_ptr{static_cast<double*>(data->initial_data_ptr)};
                const double offset{std::get<double>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                double* src_ptr{static_cast<double*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT16_T:
        {
            if(std::holds_alternative<int16_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                int16_t* dest_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                const int16_t offset{std::get<int16_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int16_t* src_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT64_T:
        {
            if(std::holds_alternative<int64_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                int64_t* dest_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                const int64_t offset{std::get<int64_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int64_t* src_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT64_T:
        {
            if(std::holds_alternative<uint64_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                uint64_t* dest_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                const uint64_t offset{std::get<uint64_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint64_t* src_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT32_T:
        {
            if(std::holds_alternative<uint32_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                uint32_t* dest_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                const uint32_t offset{std::get<uint32_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint32_t* src_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT8_T:
        {
            if(std::holds_alternative<int8_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                int8_t* dest_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                const int8_t offset{std::get<int8_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int8_t* src_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT8_T:
        {
            if(std::holds_alternative<uint8_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                uint8_t* dest_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                const uint8_t offset{std::get<uint8_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint8_t* src_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT16_T:
        {
            if(std::holds_alternative<uint16_t>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                uint16_t* dest_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                const uint16_t offset{std::get<uint16_t>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint16_t* src_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot apply scaling to data of type byte.");
        break;
        case CMC_CHAR:
        {
            if(std::holds_alternative<char>(constant_summand))
            {
                /* If the scale factor is of the same type as the data */
                char* dest_ptr{static_cast<char*>(data->initial_data_ptr)};
                const char offset{std::get<char>(constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] += offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                char* src_ptr{static_cast<char*>(data->initial_data_ptr)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, constant_summand)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = static_cast<cmc_standard_type>(offset) + static_cast<cmc_standard_type>(src_ptr[i]);
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void
var_array_t::axpy_scalar(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& add_offset)
{
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            if(std::holds_alternative<int32_t>(scale_factor) && std::holds_alternative<int32_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                int32_t* dest_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                const int32_t scaling{std::get<int32_t>(scale_factor)};
                const int32_t offset{std::get<int32_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int32_t* src_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_FLOAT:
        {
            if(std::holds_alternative<float>(scale_factor) && std::holds_alternative<float>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                float* dest_ptr{static_cast<float*>(data->initial_data_ptr)};
                const float scaling{std::get<float>(scale_factor)};
                const float offset{std::get<float>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                float* src_ptr{static_cast<float*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_DOUBLE:
        {
            if(std::holds_alternative<double>(scale_factor) && std::holds_alternative<double>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                double* dest_ptr{static_cast<double*>(data->initial_data_ptr)};
                const double scaling{std::get<double>(scale_factor)};
                const double offset{std::get<double>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                double* src_ptr{static_cast<double*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT16_T:
        {
            if(std::holds_alternative<int16_t>(scale_factor) && std::holds_alternative<int16_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                int16_t* dest_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                const int16_t scaling{std::get<int16_t>(scale_factor)};
                const int16_t offset{std::get<int16_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int16_t* src_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT64_T:
        {
            if(std::holds_alternative<int64_t>(scale_factor) && std::holds_alternative<int64_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                int64_t* dest_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                const int64_t scaling{std::get<int64_t>(scale_factor)};
                const int64_t offset{std::get<int64_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int64_t* src_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT64_T:
        {
            if(std::holds_alternative<uint64_t>(scale_factor) && std::holds_alternative<uint64_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                uint64_t* dest_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                const uint64_t scaling{std::get<uint64_t>(scale_factor)};
                const uint64_t offset{std::get<uint64_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint64_t* src_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT32_T:
        {
            if(std::holds_alternative<uint32_t>(scale_factor) && std::holds_alternative<uint32_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                uint32_t* dest_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                const uint32_t scaling{std::get<uint32_t>(scale_factor)};
                const uint32_t offset{std::get<uint32_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint32_t* src_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_INT8_T:
        {
            if(std::holds_alternative<int8_t>(scale_factor) && std::holds_alternative<int8_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                int8_t* dest_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                const int8_t scaling{std::get<int8_t>(scale_factor)};
                const int8_t offset{std::get<int8_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                int8_t* src_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT8_T:
        {
            if(std::holds_alternative<uint8_t>(scale_factor) && std::holds_alternative<uint8_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                uint8_t* dest_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                const uint8_t scaling{std::get<uint8_t>(scale_factor)};
                const uint8_t offset{std::get<uint8_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint8_t* src_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_UINT16_T:
        {
            if(std::holds_alternative<uint16_t>(scale_factor) && std::holds_alternative<uint16_t>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                uint16_t* dest_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                const uint16_t scaling{std::get<uint16_t>(scale_factor)};
                const uint16_t offset{std::get<uint16_t>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                uint16_t* src_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add to data to type byte.");
        break;
        case CMC_CHAR:
        {
            if(std::holds_alternative<char>(scale_factor) && std::holds_alternative<char>(add_offset))
            {
                /* If the scale factor is of the same type as the data */
                char* dest_ptr{static_cast<char*>(data->initial_data_ptr)};
                const char scaling{std::get<char>(scale_factor)};
                const char offset{std::get<char>(add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * dest_ptr[i] + offset;
                }
            } else
            {
                /* If the scale factor and the data are NOT of the same data type */
                /* All data will be casted to the standard type (in the hope the data will represented correct. However, this cannot be ensured) */
                var_array* data_new = new var_array{data->num_elements, get_cmc_type<cmc_standard_type>()};
                cmc_standard_type* dest_ptr{static_cast<cmc_standard_type*>(data_new->initial_data_ptr)};
                char* src_ptr{static_cast<char*>(data->initial_data_ptr)};
                const cmc_standard_type scaling{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, scale_factor)};
                const cmc_standard_type offset{std::visit([](auto& value)->cmc_standard_type {return static_cast<cmc_standard_type>(value);}, add_offset)};
                for(size_t i{0}; i < data->num_elements; ++i)
                {
                    dest_ptr[i] = scaling * static_cast<cmc_standard_type>(src_ptr[i]) + offset;
                }
                /* Deallocate the old data array */
                delete data;
                /* Assign the scaled data array */
                data = data_new;
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void
var_array_t::assign(const size_t index, const cmc_universal_type_t& _value)
{
    cmc_assert(index >= 0 && index < data->num_elements);
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> int32_t {return static_cast<int32_t>(value);}, _value);
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> float {return static_cast<float>(value);}, _value);
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> double {return static_cast<double>(value);}, _value);
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> int16_t {return static_cast<int16_t>(value);}, _value);
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> int64_t {return static_cast<int64_t>(value);}, _value);
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> uint64_t {return static_cast<uint64_t>(value);}, _value);
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> uint32_t {return static_cast<uint32_t>(value);}, _value);
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> int8_t {return static_cast<int8_t>(value);}, _value);
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> uint8_t {return static_cast<uint8_t>(value);}, _value);
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> uint16_t {return static_cast<uint16_t>(value);}, _value);
        }
        break;
        case CMC_BYTE:
        {
            std::byte* data_ptr{static_cast<std::byte*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> std::byte {return static_cast<std::byte>(value);}, _value);
        }
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            data_ptr[index] = std::visit([](auto&& value) -> char {return static_cast<char>(value);}, _value);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void
var_array_t::assign_value(const size_t index, const cmc_universal_type_t& value)
{
    cmc_assert(index >= 0 && index < data->num_elements);
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<int32_t>(value);
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<float>(value);
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<double>(value);
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<int16_t>(value);
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<int64_t>(value);
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<uint64_t>(value);
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<uint32_t>(value);
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<int8_t>(value);
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<uint8_t>(value);
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<uint16_t>(value);
        }
        break;
        case CMC_BYTE:
        {
            std::byte* data_ptr{static_cast<std::byte*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<std::byte>(value);
        }
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            data_ptr[index] = std::get<char>(value);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

cmc_universal_type_t 
var_array_t::mean_value(const size_t start_index, const size_t end_index) const
{
    cmc_assert(end_index > start_index && end_index < data->num_elements);
    switch (data->type)
    {
        case CMC_INT32_T:
            {
                return cmc_universal_type_t(std::get<int32_t>(sum_over_range(start_index, end_index)) / static_cast<int32_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_FLOAT:
            {
                return cmc_universal_type_t(std::get<float>(sum_over_range(start_index, end_index)) / static_cast<float>(end_index + 1 - start_index));
            }
        break;
        case CMC_DOUBLE:
            {
                return cmc_universal_type_t(std::get<double>(sum_over_range(start_index, end_index)) / static_cast<double>(end_index + 1 - start_index));
            }
        break;
        case CMC_INT16_T:
            {
                return cmc_universal_type_t(std::get<int16_t>(sum_over_range(start_index, end_index)) / static_cast<int16_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_INT64_T:
            {
                return cmc_universal_type_t(std::get<int64_t>(sum_over_range(start_index, end_index)) / static_cast<int64_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_UINT64_T:
            {
                return cmc_universal_type_t(std::get<uint64_t>(sum_over_range(start_index, end_index)) / static_cast<uint64_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_UINT32_T:
            {
                return cmc_universal_type_t(std::get<uint32_t>(sum_over_range(start_index, end_index)) / static_cast<uint32_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_INT8_T:
            {
                return cmc_universal_type_t(std::get<int8_t>(sum_over_range(start_index, end_index)) / static_cast<int8_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_UINT8_T:
            {
                return cmc_universal_type_t(std::get<uint8_t>(sum_over_range(start_index, end_index)) / static_cast<uint8_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_UINT16_T:
            {
                return cmc_universal_type_t(std::get<uint16_t>(sum_over_range(start_index, end_index)) / static_cast<uint16_t>(end_index + 1 - start_index));
            }
        break;
        case CMC_BYTE:
            cmc_err_msg("It is not possible to obtain a mean value of bytes.");
        break;
        case CMC_CHAR:
            {
                return cmc_universal_type_t(std::get<char>(sum_over_range(start_index, end_index)) / static_cast<char>(end_index + 1 - start_index));
            }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
            return cmc_universal_type_t(static_cast<int32_t>(0));
    }
}

/**
 * @brief Find the maximum within the given interval while considering missing values within the data
 */
cmc_universal_type_t
var_array_t::maximum_w_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const
{
    cmc_assert(end_index >= start_index && end_index < data->num_elements);
    /* It is assumed that the missing value is of the same data type as the data variable itself */
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            int32_t current_maximum{std::numeric_limits<int32_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_FLOAT:
        {
            const float missing_val{std::get<float>(missing_value)};
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            float current_maximum{std::numeric_limits<float>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_DOUBLE:
        {
            const double missing_val{std::get<double>(missing_value)};
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            double current_maximum{std::numeric_limits<double>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_INT16_T:
        {
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            int16_t current_maximum{std::numeric_limits<int16_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_INT64_T:
        {
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            int64_t current_maximum{std::numeric_limits<int64_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT64_T:
        {
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            uint64_t current_maximum{std::numeric_limits<uint64_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT32_T:
        {
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            uint32_t current_maximum{std::numeric_limits<uint32_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_INT8_T:
        {
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            int8_t current_maximum{std::numeric_limits<int8_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT8_T:
        {
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            uint8_t current_maximum{std::numeric_limits<uint8_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT16_T:
        {
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            uint16_t current_maximum{std::numeric_limits<uint16_t>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_BYTE:
        {
            const std::byte missing_val{std::get<std::byte>(missing_value)};
            std::byte* data_ptr{static_cast<std::byte*>(data->initial_data_ptr)};
            std::byte current_maximum{std::numeric_limits<std::byte>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (data_ptr[i] != missing_val)
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_CHAR:
        {
            const char missing_val{std::get<char>(missing_value)};
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            char current_maximum{std::numeric_limits<char>::min()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a greater value than currently saved */
                    if (current_maximum < data_ptr[i])
                    {
                        current_maximum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }   
}

/**
 * @brief Find the maximum within the given interval while considering missing values within the data
 */
cmc_universal_type_t
var_array_t::maximum(const size_t start_index, const size_t end_index) const
{
    cmc_assert(end_index >= start_index && end_index < data->num_elements);
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            int32_t current_maximum{std::numeric_limits<int32_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            float current_maximum{std::numeric_limits<float>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            double current_maximum{std::numeric_limits<double>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            int16_t current_maximum{std::numeric_limits<int16_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            int64_t current_maximum{std::numeric_limits<int64_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            uint64_t current_maximum{std::numeric_limits<uint64_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            uint32_t current_maximum{std::numeric_limits<uint32_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            int8_t current_maximum{std::numeric_limits<int8_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            uint8_t current_maximum{std::numeric_limits<uint8_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            uint16_t current_maximum{std::numeric_limits<uint16_t>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_BYTE:
        {
            std::byte* data_ptr{static_cast<std::byte*>(data->initial_data_ptr)};
            std::byte current_maximum{std::numeric_limits<std::byte>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            char current_maximum{std::numeric_limits<char>::min()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a greater value than currently saved */
                if (current_maximum < data_ptr[i])
                {
                    current_maximum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_maximum);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
            return cmc_universal_type_t(static_cast<int>(CMC_ERR));
    }   
}


/**
 * @brief Find the minimum within the given interval while considering missing values within the data
 */
cmc_universal_type_t
var_array_t::minimum_w_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const
{
    cmc_assert(end_index >= start_index && end_index < data->num_elements);
    /* It is assumed that the missing value is of the same data type as the data variable itself */
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            int32_t current_minimum{std::numeric_limits<int32_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_FLOAT:
        {
            const float missing_val{std::get<float>(missing_value)};
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            float current_minimum{std::numeric_limits<float>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_DOUBLE:
        {
            const double missing_val{std::get<double>(missing_value)};
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            double current_minimum{std::numeric_limits<double>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_INT16_T:
        {
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            int16_t current_minimum{std::numeric_limits<int16_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_INT64_T:
        {
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            int64_t current_minimum{std::numeric_limits<int64_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT64_T:
        {
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            uint64_t current_minimum{std::numeric_limits<uint64_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT32_T:
        {
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            uint32_t current_minimum{std::numeric_limits<uint32_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_INT8_T:
        {
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            int8_t current_minimum{std::numeric_limits<int8_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT8_T:
        {
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            uint8_t current_minimum{std::numeric_limits<uint8_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT16_T:
        {
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            uint16_t current_minimum{std::numeric_limits<uint16_t>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_BYTE:
        {
            const std::byte missing_val{std::get<std::byte>(missing_value)};
            std::byte* data_ptr{static_cast<std::byte*>(data->initial_data_ptr)};
            std::byte current_minimum{std::numeric_limits<std::byte>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (data_ptr[i] != missing_val)
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_CHAR:
        {
            const char missing_val{std::get<char>(missing_value)};
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            char current_minimum{std::numeric_limits<char>::max()};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Check for a smaller value than currently saved */
                    if (current_minimum > data_ptr[i])
                    {
                        current_minimum = data_ptr[i];
                    }
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
            return cmc_universal_type_t(static_cast<int>(CMC_ERR));
    }   
}

/**
 * @brief Find the minimum within the given interval while considering missing values within the data
 */
cmc_universal_type_t
var_array_t::minimum(const size_t start_index, const size_t end_index) const
{
    cmc_assert(end_index >= start_index && end_index < data->num_elements);
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            int32_t current_minimum{std::numeric_limits<int32_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            float current_minimum{std::numeric_limits<float>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            double current_minimum{std::numeric_limits<double>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            int16_t current_minimum{std::numeric_limits<int16_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            int64_t current_minimum{std::numeric_limits<int64_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            uint64_t current_minimum{std::numeric_limits<uint64_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            uint32_t current_minimum{std::numeric_limits<uint32_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            int8_t current_minimum{std::numeric_limits<int8_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            uint8_t current_minimum{std::numeric_limits<uint8_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            uint16_t current_minimum{std::numeric_limits<uint16_t>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_BYTE:
        {
            std::byte* data_ptr{static_cast<std::byte*>(data->initial_data_ptr)};
            std::byte current_minimum{std::numeric_limits<std::byte>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            char current_minimum{std::numeric_limits<char>::max()};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Check for a smaller value than currently saved */
                if (current_minimum > data_ptr[i])
                {
                    current_minimum = data_ptr[i];
                }
            }
            return cmc_universal_type_t(current_minimum);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
            return cmc_universal_type_t(static_cast<int>(CMC_ERR));
    }   
}

/**
 * @brief Calculate the mean value of the given interval considering the presence of missing values
 */
cmc_universal_type_t
var_array_t::mean_value_same_type_wo_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const
{
    cmc_assert(end_index >= start_index && end_index < data->num_elements);
    /* It is assumed that the missing value is of the same data type as the data variable itself */
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            std::pair<int32_t, int32_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<int32_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<int32_t>(0));
            }
        }
        break;
        case CMC_FLOAT:
        {
            const float missing_val{std::get<float>(missing_value)};
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            std::pair<float, float> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<float>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<float>(0));
            }
        }
        break;
        case CMC_DOUBLE:
        {
            const double missing_val{std::get<double>(missing_value)};
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            std::pair<double, double> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<double>(0));
            }
        }
        break;
        case CMC_INT16_T:
        {
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            std::pair<int16_t, int16_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<int16_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<int16_t>(0));
            }
        }
        break;
        case CMC_INT64_T:
        {
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            std::pair<int64_t, int64_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<int64_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<int64_t>(0));
            }
        }
        break;
        case CMC_UINT64_T:
        {
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            std::pair<uint64_t, uint64_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<uint64_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<uint64_t>(0));
            }
        }
        break;
        case CMC_UINT32_T:
        {
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            std::pair<uint32_t, uint32_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<uint32_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<uint32_t>(0));
            }
        }
        break;
        case CMC_INT8_T:
        {
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            std::pair<int8_t, int8_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<int8_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<int8_t>(0));
            }
        }
        break;
        case CMC_UINT8_T:
        {
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            std::pair<uint8_t, uint8_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<uint8_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<uint8_t>(0));
            }
        }
        break;
        case CMC_UINT16_T:
        {
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            std::pair<uint16_t, uint16_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<uint16_t>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<uint16_t>(0));
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            const char missing_val{std::get<char>(missing_value)};
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            std::pair<char, char> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return cmc_universal_type_t(static_cast<char>(mean_no_missing_vals.first / mean_no_missing_vals.second));
            } else
            {
                return cmc_universal_type_t(static_cast<char>(0));
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }   
}


double
var_array_t::mean_value_wo_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const
{
    cmc_assert(end_index > start_index && end_index < data->num_elements);
    /* It is assumed that the missing value is of the same data type as the data variable itself */
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            std::pair<int32_t, int32_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_FLOAT:
        {
            const float missing_val{std::get<float>(missing_value)};
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            std::pair<float, float> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_DOUBLE:
        {
            const double missing_val{std::get<double>(missing_value)};
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            std::pair<double, double> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_INT16_T:
        {
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            std::pair<int16_t, int16_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_INT64_T:
        {
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            std::pair<int64_t, int64_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_UINT64_T:
        {
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            std::pair<uint64_t, uint64_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_UINT32_T:
        {
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            std::pair<uint32_t, uint32_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_INT8_T:
        {
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            std::pair<int8_t, int8_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_UINT8_T:
        {
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            std::pair<uint8_t, uint8_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_UINT16_T:
        {
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            std::pair<uint16_t, uint16_t> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            const char missing_val{std::get<char>(missing_value)};
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            std::pair<char, char> mean_no_missing_vals{std::make_pair(0,0)};
            /* Calculate the sum and the number of no missing values */
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    mean_no_missing_vals.first += data_ptr[i];
                    ++(mean_no_missing_vals.second);
                }
            }
            /* If not all values are missing values, return the mean value */
            if (mean_no_missing_vals.second != 0)
            {
                return static_cast<double>(mean_no_missing_vals.first / mean_no_missing_vals.second);
            } else
            {
                return 0.0;
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }   
}

std::vector<double>
var_array_t::calculate_relative_deviations_w_missing_values(const size_t start_index, const size_t end_index, const cmc_universal_type_t& nominal_value, const cmc_universal_type_t& missing_value) const
{
    /* Is the data to compare of the same type */
    cmc_assert(data->type == static_cast<cmc_type>(nominal_value.index()));
    /* Do we have an increasing range */
    cmc_assert(end_index >= start_index);

    /* Declare an output vector containing the relative deviations of the nominal value to data from the array */
    std::vector<double> deviations;
    deviations.reserve(end_index - start_index + 1);

    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            const int32_t reference_val{std::get<int32_t>(nominal_value)};
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(data_ptr[i])));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                        }
                    }
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            const float reference_val{std::get<float>(nominal_value)};
            const float missing_val{std::get<float>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!cmc_approx_cmp(data_ptr[i], static_cast<float>(0.0)))
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(data_ptr[i])));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (cmc_approx_cmp(reference_val, static_cast<float>(0.0)))
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                        }
                    }
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            const double reference_val{std::get<double>(nominal_value)};
            const double missing_val{std::get<double>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!cmc_approx_cmp(data_ptr[i], 0.0))
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(std::abs(data_ptr[i] - reference_val) / std::abs(data_ptr[i]));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (cmc_approx_cmp(reference_val, 0.0))
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(std::abs(data_ptr[i] - reference_val) / ((std::abs(data_ptr[i]) + std::abs(reference_val)) * 0.5));
                        }
                    }
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            const int16_t reference_val{std::get<int16_t>(nominal_value)};
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(data_ptr[i])));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                        }
                    }
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            const int64_t reference_val{std::get<int64_t>(nominal_value)};
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(data_ptr[i])));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                        }
                    }
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            const uint64_t reference_val{std::get<uint64_t>(nominal_value)};
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        const uint64_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(data_ptr[i]));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            const uint64_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                        }
                    }
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            const uint32_t reference_val{std::get<uint32_t>(nominal_value)};
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        const uint32_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(data_ptr[i]));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            const uint32_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                        }
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            const int8_t reference_val{std::get<int8_t>(nominal_value)};
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(data_ptr[i])));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                        }
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            const uint8_t reference_val{std::get<uint8_t>(nominal_value)};
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        const uint8_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(data_ptr[i]));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            const uint8_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                        }
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            const uint16_t reference_val{std::get<uint16_t>(nominal_value)};
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        const uint16_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(data_ptr[i]));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            const uint16_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                        }
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            const char reference_val{std::get<char>(nominal_value)};
            const char missing_val{std::get<char>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] != 0)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(data_ptr[i])));
                    } else
                    {
                        /* Check if the numerator equals zero as well */
                        if (reference_val == 0)
                        {
                            /* Just save a zero in this case */
                            deviations.push_back(0.0);
                        } else
                        {
                            /* In case the reference value is zero, we use an indicator of relative difference */
                            deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                        }
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }

    return deviations;
}

std::vector<double>
var_array_t::calculate_absolute_deviations_w_missing_values(const size_t start_index, const size_t end_index, const cmc_universal_type_t& nominal_value, const cmc_universal_type_t& missing_value) const
{
    /* Is the data to compare of the same type */
    cmc_assert(data->type == static_cast<cmc_type>(nominal_value.index()));
    /* Do we have an increasing range */
    cmc_assert(end_index >= start_index);

    /* Declare an output vector containing the absolute deviations of the nominal value to data from the array */
    std::vector<double> deviations;
    deviations.reserve(end_index - start_index + 1);

    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            const int32_t reference_val{std::get<int32_t>(nominal_value)};
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            const float reference_val{std::get<float>(nominal_value)};
            const float missing_val{std::get<float>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            const double reference_val{std::get<double>(nominal_value)};
            const double missing_val{std::get<double>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(std::abs(data_ptr[i] - reference_val));
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            const int16_t reference_val{std::get<int16_t>(nominal_value)};
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            const int64_t reference_val{std::get<int64_t>(nominal_value)};
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            const uint64_t reference_val{std::get<uint64_t>(nominal_value)};
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
                } else 
                {
                    /* Missing values are not considered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            const uint32_t reference_val{std::get<uint32_t>(nominal_value)};
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            const int8_t reference_val{std::get<int8_t>(nominal_value)};
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            const uint8_t reference_val{std::get<uint8_t>(nominal_value)};
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            const uint16_t reference_val{std::get<uint16_t>(nominal_value)};
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            const char reference_val{std::get<char>(nominal_value)};
            const char missing_val{std::get<char>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    /* Calculate the absolute deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }

    return deviations;
}

std::vector<double>
var_array_t::calculate_absolute_deviations(const size_t start_index, const size_t end_index, const cmc_universal_type_t& nominal_value) const
{
    /* Is the data to compare of the same type */
    cmc_assert(data->type == static_cast<cmc_type>(nominal_value.index()));
    /* Do we have an increasing range */
    cmc_assert(end_index >= start_index);

    /* Declare an output vector containing the absolute deviations of the nominal value to data from the array */
    std::vector<double> deviations;
    deviations.reserve(end_index - start_index + 1);

    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            const int32_t reference_val{std::get<int32_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
            }
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            const float reference_val{std::get<float>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
            }
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            const double reference_val{std::get<double>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(std::abs(data_ptr[i] - reference_val));
            }
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            const int16_t reference_val{std::get<int16_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
            }
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            const int64_t reference_val{std::get<int64_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
            }
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            const uint64_t reference_val{std::get<uint64_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
            }
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            const uint32_t reference_val{std::get<uint32_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
            }
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            const int8_t reference_val{std::get<int8_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
            }
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            const uint8_t reference_val{std::get<uint8_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
            }
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            const uint16_t reference_val{std::get<uint16_t>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>((data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i])));
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            const char reference_val{std::get<char>(nominal_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                /* Calculate the absolute deviation */
                deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)));
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }

    return deviations;
}


#if 0
std::vector<double>
var_array_t::calculate_relative_deviations_w_missing_values(const size_t start_index, const size_t end_index, const cmc_universal_type_t& nominal_value, const cmc_universal_type_t& missing_value) const
{
    cmc_assert(data->type == static_cast<cmc_type>(nominal_value.index()));
    cmc_assert(end_index >= start_index);

    std::vector<double> deviations;
    deviations.reserve(end_index - start_index + 1);

    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            const int32_t reference_val{std::get<int32_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                    } else
                    {
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            const float reference_val{std::get<float>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const float missing_val{std::get<float>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                    } else
                    {
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            const double reference_val{std::get<double>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const double missing_val{std::get<double>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(std::abs(data_ptr[i] - reference_val) / std::abs(data_ptr[i]));
                    } else
                    {
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(std::abs(data_ptr[i] - reference_val) / ((std::abs(data_ptr[i]) + std::abs(reference_val)) * 0.5));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            const int16_t reference_val{std::get<int16_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                    } else
                    {
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            const int64_t reference_val{std::get<int64_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                    } else
                    {
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            const uint64_t reference_val{std::get<uint64_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        const uint64_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                    } else
                    {
                        const uint64_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            const uint32_t reference_val{std::get<uint32_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        const uint32_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                    } else
                    {
                        const uint32_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            const int8_t reference_val{std::get<int8_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                    } else
                    {
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            const uint8_t reference_val{std::get<uint8_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        const uint8_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                    } else
                    {
                        const uint8_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            const uint16_t reference_val{std::get<uint16_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        const uint16_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                    } else
                    {
                        const uint16_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            const char reference_val{std::get<char>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            const char missing_val{std::get<char>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!reference_val_equal_to_zero)
                    {
                        /* Calculate the relative deviation */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                    } else
                    {
                        /* In case the reference value is zero, we use an indicator of relative difference */
                        deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                    }
                } else 
                {
                    /* Missing values are not cosnidered, therefore there is no deviation */
                    deviations.push_back(0.0);
                }
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }

    return deviations;
}
#endif

std::vector<double>
var_array_t::calculate_relative_deviations(const size_t start_index, const size_t end_index, const cmc_universal_type_t& nominal_value) const
{
    cmc_assert(data->type == static_cast<cmc_type>(nominal_value.index()));
    cmc_assert(end_index >= start_index);

    std::vector<double> deviations(end_index - start_index + 1);

    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            const int32_t reference_val{std::get<int32_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {

                if (!reference_val_equal_to_zero)
                {
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                } else
                {
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                }
            }
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            const float reference_val{std::get<float>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                } else
                {
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                }

            }
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            const double reference_val{std::get<double>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    /* Calculate the relative deviation */
                    deviations.push_back(std::abs(data_ptr[i] - reference_val) / std::abs(reference_val));
                } else
                {
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(std::abs(data_ptr[i] - reference_val) / ((std::abs(data_ptr[i]) + std::abs(reference_val)) * 0.5));
                }
            }
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            const int16_t reference_val{std::get<int16_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                } else
                {
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                }
            }
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            const int64_t reference_val{std::get<int64_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                } else
                {
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                }
            }
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            const uint64_t reference_val{std::get<uint64_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    const uint64_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                } else
                {
                    const uint64_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                }
            }
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            const uint32_t reference_val{std::get<uint32_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    const uint32_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                } else
                {
                    const uint32_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                }
            }
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            const int8_t reference_val{std::get<int8_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                } else
                {
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                }
            }
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            const uint8_t reference_val{std::get<uint8_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    const uint8_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                } else
                {
                    const uint8_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                }
            }
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            const uint16_t reference_val{std::get<uint16_t>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    const uint16_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(dividend) / static_cast<double>(reference_val));
                } else
                {
                    const uint16_t dividend = (data_ptr[i] >= reference_val ? data_ptr[i] - reference_val : reference_val - data_ptr[i]);
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(dividend) / ((static_cast<double>(data_ptr[i] + reference_val) * 0.5)));
                }
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            const char reference_val{std::get<char>(nominal_value)};
            const bool reference_val_equal_to_zero{reference_val == 0};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!reference_val_equal_to_zero)
                {
                    /* Calculate the relative deviation */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / static_cast<double>(std::abs(reference_val)));
                } else
                {
                    /* In case the reference value is zero, we use an indicator of relative difference */
                    deviations.push_back(static_cast<double>(std::abs(data_ptr[i] - reference_val)) / (static_cast<double>((std::abs(data_ptr[i]) + std::abs(reference_val))) * 0.5));
                }
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }

    return deviations;
}

bool
var_array_t::check_range_fullfills_deviation_threshold_from_value(const size_t start_index, const size_t end_index, const double deviation, const cmc_universal_type_t& value, const cmc_universal_type_t& missing_value) const
{
    cmc_assert(end_index > start_index && end_index < data->num_elements && deviation >= 0.0);

   switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            const int32_t reference_val{std::get<int32_t>(value)};
            const int32_t missing_val{std::get<int32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!(std::abs(data_ptr[i] - reference_val) < std::abs(data_ptr[i]) * deviation))
                    {
                        /* The threshold is not fullfilled */
                        return false;
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            const float reference_val{std::get<float>(value)};
            const float missing_val{std::get<float>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!(std::abs(data_ptr[i] - reference_val) < std::abs(data_ptr[i]) * deviation))
                    {
                        /* The threshold is not fullfilled */
                        return false;
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            const double reference_val{std::get<double>(value)};
            const double missing_val{std::get<double>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!(std::abs(data_ptr[i] - reference_val) < std::abs(data_ptr[i]) * deviation))
                    {
                        /* The threshold is not fullfilled */
                        return false;
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            const int16_t reference_val{std::get<int16_t>(value)};
            const int16_t missing_val{std::get<int16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!(std::abs(data_ptr[i] - reference_val) < std::abs(data_ptr[i]) * deviation))
                    {
                        /* The threshold is not fullfilled */
                        return false;
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            const int64_t reference_val{std::get<int64_t>(value)};
            const int64_t missing_val{std::get<int64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!(std::abs(data_ptr[i] - reference_val) < std::abs(data_ptr[i]) * deviation))
                    {
                        /* The threshold is not fullfilled */
                        return false;
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            const uint64_t reference_val{std::get<uint64_t>(value)};
            const uint64_t missing_val{std::get<uint64_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] >= reference_val)
                    {
                        if (!(data_ptr[i] - reference_val < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    } else
                    {
                        if (!(reference_val - data_ptr[i] < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            const uint32_t reference_val{std::get<uint32_t>(value)};
            const uint32_t missing_val{std::get<uint32_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] >= reference_val)
                    {
                        if (!(data_ptr[i] - reference_val < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    } else
                    {
                        if (!(reference_val - data_ptr[i] < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            const int8_t reference_val{std::get<int8_t>(value)};
            const int8_t missing_val{std::get<int8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!(std::abs(data_ptr[i] - reference_val) < std::abs(data_ptr[i]) * deviation))
                    {
                        /* The threshold is not fullfilled */
                        return false;
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            const uint8_t reference_val{std::get<uint8_t>(value)};
            const uint8_t missing_val{std::get<uint8_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] >= reference_val)
                    {
                        if (!(data_ptr[i] - reference_val < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    } else
                    {
                        if (!(reference_val - data_ptr[i] < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    }
                }
            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            const uint16_t reference_val{std::get<uint16_t>(value)};
            const uint16_t missing_val{std::get<uint16_t>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (data_ptr[i] >= reference_val)
                    {
                        if (!(data_ptr[i] - reference_val < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    } else
                    {
                        if (!(reference_val - data_ptr[i] < data_ptr[i] * deviation))
                        {
                            /* The threshold is not fullfilled */
                            return false;
                        }
                    }
                }
            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot add data of type byte.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            const char reference_val{std::get<char>(value)};
            const char missing_val{std::get<char>(missing_value)};
            for (size_t i{start_index}; i <= end_index; ++i)
            {
                if (!cmc_approx_cmp(data_ptr[i], missing_val))
                {
                    if (!(std::abs(data_ptr[i] - reference_val) < std::abs(data_ptr[i]) * deviation))
                    {
                        /* The threshold is not fullfilled */
                        return false;
                    }
                }

            }
            /* If all checks have passed, all values (except missing values) comply to the given threshold */
            return true;
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
            return false;
    }
}

cmc_type
var_array_t::get_data_type() const
{
    return data->type;
}

cmc_universal_type_t
var_array_t::operator[](const size_t index) const
{
    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_BYTE:
        {
            std::byte* data_ptr{static_cast<std::byte*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            return cmc_universal_type_t{data_ptr[index]};
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
    return cmc_universal_type_t{0};
}

void
var_array_t::print_data() const
{  
    switch(data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot print data of type byte to standard ostream.");
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                std::cout << data_ptr[i] << ", ";
            }
            std::cout << std::endl;
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void
var_array_t::copy_from_to(const var_array_t& source_array, const size_t src_start_index, const size_t src_end_index, const size_t dest_start_index)
{
    cmc_assert(src_end_index >= src_start_index);
    cmc_assert(source_array.size() > src_end_index);
    cmc_assert(data->num_elements - dest_start_index >= src_end_index - src_start_index);
    cmc_assert(data->type == source_array.get_data_type());

    switch (data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* dest_data_ptr{static_cast<int32_t*>(data->initial_data_ptr) + dest_start_index};
            int32_t* src_data_ptr{static_cast<int32_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            int32_t* src_data_end_ptr{static_cast<int32_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_FLOAT:
        {
            float* dest_data_ptr{static_cast<float*>(data->initial_data_ptr) + dest_start_index};
            float* src_data_ptr{static_cast<float*>(source_array.get_initial_data_ptr()) + src_start_index};
            float* src_data_end_ptr{static_cast<float*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_DOUBLE:
        {
            double* dest_data_ptr{static_cast<double*>(data->initial_data_ptr) + dest_start_index};
            double* src_data_ptr{static_cast<double*>(source_array.get_initial_data_ptr()) + src_start_index};
            double* src_data_end_ptr{static_cast<double*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* dest_data_ptr{static_cast<int16_t*>(data->initial_data_ptr) + dest_start_index};
            int16_t* src_data_ptr{static_cast<int16_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            int16_t* src_data_end_ptr{static_cast<int16_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* dest_data_ptr{static_cast<int64_t*>(data->initial_data_ptr) + dest_start_index};
            int64_t* src_data_ptr{static_cast<int64_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            int64_t* src_data_end_ptr{static_cast<int64_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* dest_data_ptr{static_cast<uint64_t*>(data->initial_data_ptr) + dest_start_index};
            uint64_t* src_data_ptr{static_cast<uint64_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            uint64_t* src_data_end_ptr{static_cast<uint64_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* dest_data_ptr{static_cast<uint32_t*>(data->initial_data_ptr) + dest_start_index};
            uint32_t* src_data_ptr{static_cast<uint32_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            uint32_t* src_data_end_ptr{static_cast<uint32_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* dest_data_ptr{static_cast<int8_t*>(data->initial_data_ptr) + dest_start_index};
            int8_t* src_data_ptr{static_cast<int8_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            int8_t* src_data_end_ptr{static_cast<int8_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* dest_data_ptr{static_cast<uint8_t*>(data->initial_data_ptr) + dest_start_index};
            uint8_t* src_data_ptr{static_cast<uint8_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            uint8_t* src_data_end_ptr{static_cast<uint8_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* dest_data_ptr{static_cast<uint16_t*>(data->initial_data_ptr) + dest_start_index};
            uint16_t* src_data_ptr{static_cast<uint16_t*>(source_array.get_initial_data_ptr()) + src_start_index};
            uint16_t* src_data_end_ptr{static_cast<uint16_t*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_BYTE:
        {
            std::byte* dest_data_ptr{static_cast<std::byte*>(data->initial_data_ptr) + dest_start_index};
            std::byte* src_data_ptr{static_cast<std::byte*>(source_array.get_initial_data_ptr()) + src_start_index};
            std::byte* src_data_end_ptr{static_cast<std::byte*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        case CMC_CHAR:
        {
            char* dest_data_ptr{static_cast<char*>(data->initial_data_ptr) + dest_start_index};
            char* src_data_ptr{static_cast<char*>(source_array.get_initial_data_ptr()) + src_start_index};
            char* src_data_end_ptr{static_cast<char*>(source_array.get_initial_data_ptr()) + src_end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void
var_array_t::crop_to(const size_t start_index, const size_t end_index)
{
    cmc_assert(end_index >= start_index);
    cmc_assert(data->num_elements > end_index);

    switch (data->type)
    {
        case CMC_INT32_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            int32_t* dest_data_ptr{static_cast<int32_t*>(data_new->initial_data_ptr)};
            int32_t* src_data_ptr{static_cast<int32_t*>(data->initial_data_ptr) + start_index};
            int32_t* src_data_end_ptr{static_cast<int32_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_FLOAT:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            float* dest_data_ptr{static_cast<float*>(data_new->initial_data_ptr)};
            float* src_data_ptr{static_cast<float*>(data->initial_data_ptr) + start_index};
            float* src_data_end_ptr{static_cast<float*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_DOUBLE:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            double* dest_data_ptr{static_cast<double*>(data_new->initial_data_ptr)};
            double* src_data_ptr{static_cast<double*>(data->initial_data_ptr) + start_index};
            double* src_data_end_ptr{static_cast<double*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_INT16_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            int16_t* dest_data_ptr{static_cast<int16_t*>(data_new->initial_data_ptr)};
            int16_t* src_data_ptr{static_cast<int16_t*>(data->initial_data_ptr) + start_index};
            int16_t* src_data_end_ptr{static_cast<int16_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_INT64_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            int64_t* dest_data_ptr{static_cast<int64_t*>(data_new->initial_data_ptr)};
            int64_t* src_data_ptr{static_cast<int64_t*>(data->initial_data_ptr) + start_index};
            int64_t* src_data_end_ptr{static_cast<int64_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_UINT64_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            uint64_t* dest_data_ptr{static_cast<uint64_t*>(data_new->initial_data_ptr)};
            uint64_t* src_data_ptr{static_cast<uint64_t*>(data->initial_data_ptr) + start_index};
            uint64_t* src_data_end_ptr{static_cast<uint64_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_UINT32_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            uint32_t* dest_data_ptr{static_cast<uint32_t*>(data_new->initial_data_ptr)};
            uint32_t* src_data_ptr{static_cast<uint32_t*>(data->initial_data_ptr) + start_index};
            uint32_t* src_data_end_ptr{static_cast<uint32_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_INT8_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            int8_t* dest_data_ptr{static_cast<int8_t*>(data_new->initial_data_ptr)};
            int8_t* src_data_ptr{static_cast<int8_t*>(data->initial_data_ptr) + start_index};
            int8_t* src_data_end_ptr{static_cast<int8_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_UINT8_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            uint8_t* dest_data_ptr{static_cast<uint8_t*>(data_new->initial_data_ptr)};
            uint8_t* src_data_ptr{static_cast<uint8_t*>(data->initial_data_ptr) + start_index};
            uint8_t* src_data_end_ptr{static_cast<uint8_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_UINT16_T:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            uint16_t* dest_data_ptr{static_cast<uint16_t*>(data_new->initial_data_ptr)};
            uint16_t* src_data_ptr{static_cast<uint16_t*>(data->initial_data_ptr) + start_index};
            uint16_t* src_data_end_ptr{static_cast<uint16_t*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_BYTE:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            std::byte* dest_data_ptr{static_cast<std::byte*>(data_new->initial_data_ptr)};
            std::byte* src_data_ptr{static_cast<std::byte*>(data->initial_data_ptr) + start_index};
            std::byte* src_data_end_ptr{static_cast<std::byte*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        case CMC_CHAR:
        {
            var_array* data_new = new var_array{end_index - start_index + 1, data->type};
            char* dest_data_ptr{static_cast<char*>(data_new->initial_data_ptr)};
            char* src_data_ptr{static_cast<char*>(data->initial_data_ptr) + start_index};
            char* src_data_end_ptr{static_cast<char*>(data->initial_data_ptr) + end_index + 1};
            std::copy(src_data_ptr, src_data_end_ptr, dest_data_ptr);
            delete data;
            data = data_new;
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
}

void*
var_array_t::get_initial_data_ptr() const
{
    return data->initial_data_ptr;
}


/*********************************************/
/*****    VAR_DYNAMIC_ARRAY_T STUFF       ****/

struct var_dynamic_array
{
public:
    var_dynamic_array(const cmc_type _type)
    {
        /* If no size is specified, we just allocate for 8 values */
        //array = new var_array_t{8, _type};
        array = new var_array(8, _type);

    };
    var_dynamic_array(const size_t _num_elems, const cmc_type _type)
    : capacity{_num_elems}{
        array = new var_array(_num_elems, _type);
    };
    ~var_dynamic_array(){
        if (array != nullptr)
        {
            delete array;
            array = nullptr;
        }
    };

    size_t capacity{8}; //!< The capacity of the dynamic array, if no capacity is specified, by default the array is capable of holding 8 elements
    size_t size{0}; //!< The current size of the dynamic array
    var_array* array{nullptr};
};

void
var_dynamic_array_t::_init_array(const size_t num_elements, const cmc_type _type)
{
    data = new var_dynamic_array(num_elements, _type);
}

void
var_dynamic_array_t::_destroy_array()
{
    if (data != nullptr)
    {
        delete data;
        data = nullptr;
    }
}

void*
var_dynamic_array_t::get_initial_data_ptr() const
{
    return data->array->initial_data_ptr;
}

size_t
var_dynamic_array_t::size() const
{
    return data->size;
}

size_t
var_dynamic_array_t::capacity() const
{
    return data->capacity;
}

cmc_type
var_dynamic_array_t::get_type() const
{
    return data->array->type;
}

void 
var_dynamic_array_t::resize(const size_t num_elements)
{
    if (num_elements > data->capacity)
    {
        /* Create a new array with the given number of elements */
        var_array* resized_data = new var_array(num_elements, data->array->type);
        /* Copy the current array entries over */
        memcpy(resized_data->initial_data_ptr, data->array->initial_data_ptr, data->size * cmc_type_to_bytes[data->array->type]);
        /* The current size stays the same, but the capacity changes */
        data->capacity = num_elements;
        /* Delete the old data */
        data->array->~var_array();
        /* Save the new resized array */
        data->array = resized_data;
    }
}

void
var_dynamic_array_t::push_back(const cmc_universal_type_t& value)
{
    /* Check if the dynamic array is still large enough */
    if (data->size >= data->capacity)
    {
        /* Resize the array (Double the current memory) */
        this->resize(2 * data->capacity);
    }

    /* Push back the supplied value based on the type of the dynamic array */
    switch (data->array->type)
    {
        case CMC_INT32_T:
        {
            int32_t* data_ptr{static_cast<int32_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<int32_t>(value);
        }
        break;
        case CMC_FLOAT:
        {
            float* data_ptr{static_cast<float*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<float>(value);
        }
        break;
        case CMC_DOUBLE:
        {
            double* data_ptr{static_cast<double*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<double>(value);
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* data_ptr{static_cast<int16_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<int16_t>(value);
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* data_ptr{static_cast<int64_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<int64_t>(value);
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* data_ptr{static_cast<uint64_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<uint64_t>(value);
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* data_ptr{static_cast<uint32_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<uint32_t>(value);
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* data_ptr{static_cast<int8_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<int8_t>(value);
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* data_ptr{static_cast<uint8_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<uint8_t>(value);
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* data_ptr{static_cast<uint16_t*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<uint16_t>(value);
        }
        break;
        case CMC_BYTE:
        {
            std::byte* data_ptr{static_cast<std::byte*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<std::byte>(value);
        }
        break;
        case CMC_CHAR:
        {
            char* data_ptr{static_cast<char*>(data->array->initial_data_ptr)};
            data_ptr[data->size] = std::get<char>(value);
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }

    /* Increment the size */
    ++(data->size);
}

/*********************************************/
/*****      VAR_VECTOR_T STUFF            ****/

struct var_vector
{
private:

public:
    var_vector(){};
    ~var_vector(){
        for (size_t i{0}; i < vector.size(); ++i)
        {
            delete vector[i];
        }
    };
    std::vector<var_array_t*> vector;
};

/* Move constructor */
var_vector_t::var_vector_t(var_vector_t&& vector)
{
    /* Create a new vector */
    vec = new var_vector();
    //vec->vector.reserve(vector.vec->vector.size());
    /* Move the vector with the pointers to the heap allocated arrays */
    vec->vector = std::move(vector.vec->vector);
};

void
var_vector_t::_init_vector()
{
    vec = new var_vector();
}

void
var_vector_t::_destroy_vector()
{
    delete vec;
}

void
var_vector_t::clear()
{
    for (size_t i{0}; i < vec->vector.size(); ++i)
    {
        delete vec->vector[i];
    }
    vec->vector.clear();
}
size_t
var_vector_t::size() const
{
    return vec->vector.size();
}

size_t
var_vector_t::capacity() const
{
    return vec->vector.capacity();
}

auto
var_vector_t::begin() const
{
    return vec->vector.begin();
}

auto
var_vector_t::end() const
{
    return vec->vector.end();
}

var_array_t& var_vector_t::operator[](const size_t index)
{
    return *(vec->vector[index]);
}

void
var_vector_t::create_and_push_back(const size_t num_elements, const cmc_type data_type)
{
    vec->vector.push_back(new var_array_t{num_elements, data_type});
}

void
var_vector_t::push_back(var_array_t* array)
{
    vec->vector.push_back(array);
}

void
var_vector_t::push_back_placeholder()
{
    vec->vector.push_back(new var_array_t{1, CMC_INT32_T});
}

void
var_vector_t::pop_back()
{
    vec->vector.pop_back();
}

void
var_vector_t::reserve(const size_t num_arrays)
{
    vec->vector.reserve(num_arrays);
}

std::vector<double>
var_array_t::retrieve_double_vector() const
{
    std::vector<double> transformed_vector;
    transformed_vector.reserve(data->num_elements);

    switch(data->type)
    {
        case CMC_INT32_T:
        {
            int32_t* ptr{static_cast<int32_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_FLOAT:
        {
            float* ptr{static_cast<float*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_DOUBLE:
        {
            double* ptr{static_cast<double*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_INT16_T:
        {
            int16_t* ptr{static_cast<int16_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_INT64_T:
        {
            int64_t* ptr{static_cast<int64_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_UINT64_T:
        {
            uint64_t* ptr{static_cast<uint64_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_UINT32_T:
        {
            uint32_t* ptr{static_cast<uint32_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_INT8_T:
        {
            int8_t* ptr{static_cast<int8_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_UINT8_T:
        {
            uint8_t* ptr{static_cast<uint8_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_UINT16_T:
        {
            uint16_t* ptr{static_cast<uint16_t*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        case CMC_BYTE:
            cmc_err_msg("Cannot transform byte.");
        break;
        case CMC_CHAR:
        {
            char* ptr{static_cast<char*>(data->initial_data_ptr)};
            for(size_t i{0}; i < data->num_elements; ++i)
            {
                transformed_vector.push_back(static_cast<double>(ptr[i]));
            }
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }

    return transformed_vector;
}
