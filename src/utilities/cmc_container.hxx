#ifndef CMC_CONTAINER_HXX
#define CMC_CONTAINER_HXX

#include "cmc_container.h"


/* Evaulation Functions of cmc_universal_data */
template<cmc_type d_type>
auto
eval(const cmc_universal_type_t& datum)
    -> std::enable_if_t<d_type == CMC_INT32_T, int32_t>
{
    return std::get<d_type>(datum);
}

template<cmc_type d_type>
auto
eval(const cmc_universal_type_t& datum)
    -> std::enable_if_t<d_type == CMC_FLOAT, float>
{
    return std::get<d_type>(datum);
}

template<cmc_type d_type>
auto
eval(const cmc_universal_type_t& datum)
    -> std::enable_if_t<d_type == CMC_DOUBLE, double>
{
    return std::get<d_type>(datum);
}

template<typename T>
void eval(T& dest, void* data)
{
    dest = static_cast<T>(data);
}
/* Evaluation of data_src, the result is written in data_dest */
auto
cmc_eval(void* data_dest, void* data_src, const cmc_type d_type);

/* Evaluation of data passed as void pointer */
template<cmc_type d_type>
auto eval(void* data)
    -> std::enable_if_t<d_type == CMC_INT32_T, int32_t>
{
    return *(static_cast<int32_t*>(data));
}

template<cmc_type d_type>
auto eval(void* data)
    -> std::enable_if_t<d_type == CMC_FLOAT, float>
{
    return *(static_cast<float*>(data));
}

template<cmc_type d_type>
auto eval(void* data)
    -> std::enable_if_t<d_type == CMC_DOUBLE, double>
{
    return *(static_cast<double*>(data));
}

/* Wrapper functions for internal operating functions */
template<typename T>
auto
var_array_t::axpy(T scale_factor, var_array_t& summand)
    -> std::enable_if_t<std::is_arithmetic_v<T>, void>
{
    _intern_axpy(static_cast<void*>(scale_factor), summand);
}

template<typename T>
auto
var_array_t::axpy_scalar(T scale_factor, T constant_summand)
    -> std::enable_if_t<std::is_arithmetic_v<T>, void>
{
    _intern_axpy_scalar(static_cast<void*>(&scale_factor), static_cast<void*>(&constant_summand));
}

template<typename T>
auto
var_array_t::scale(T scale_factor)
    -> std::enable_if_t<std::is_arithmetic_v<T>, void>
{
    _intern_scale(static_cast<void*>(&scale_factor));
}

template<typename T>
void
var_array_t::assign(const size_t index, T value)
{
    //asseert t ungleich data type
    _intern_assign(index, static_cast<void*>(&value));
}

template<typename T>
auto
var_array_t::add_const(T const_summand)
    -> std::enable_if_t<std::is_arithmetic_v<T>, void>
{
    _intern_add_const(static_cast<void*>(&const_summand));
}


template<typename T>
cmc_type
type_to_cmc_type(T val)
{
    if constexpr (std::is_same_v<std::remove_cv_t<T>, std::byte>)
    { 
        return CMC_BYTE;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, int8_t>)
    {
        return CMC_INT8_T;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, char>)
    {
        return CMC_CHAR;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, int16_t>)
    {
        return CMC_INT16_T;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, int32_t>)
    {
        return CMC_INT32_T;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, float>)
    {
        return CMC_FLOAT;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, double>)
    {
        return CMC_DOUBLE;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, uint8_t>)
    {
        return CMC_UINT8_T;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, uint16_t>)
    {
        return CMC_UINT16_T;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, uint32_t>)
    {
        return CMC_UINT32_T;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, int64_t>)
    {
        return CMC_INT64_T;
    }
    else if constexpr (std::is_same_v<std::remove_cv_t<T>, uint64_t>)
    {
        return CMC_UINT64_T;
    }
    else
    {
        return CMC_UNDEFINED;
    }
}


#endif /* CMC_CONTAINER_HXX */