#ifndef CMC_CONTAINER_H
#define CMC_CONTAINER_H

#include "cmc.h"
#include "utilities/cmc_util.h"
#include "utilities/cmc_log_functions.h"

/* Array holding the data_size (in bytes) of each variable data type */ 
inline constexpr std::array<size_t, cmc_type::CMC_NUM_TYPES> cmc_type_to_bytes{sizeof(std::byte), sizeof(int8_t), sizeof(char), sizeof(int16_t), sizeof(int32_t), sizeof(float), sizeof(double), sizeof(uint8_t), sizeof(uint16_t), sizeof(uint32_t), sizeof(int64_t), sizeof(uint64_t)};

/* Define a standard type for the data, in case of conflicts, data will be casted to this type */
typedef double cmc_standard_type;

class var_array_t
{
public:
    var_array_t(const size_t num_elements, const cmc_type _type){
        _init_array(num_elements, _type);
    };
    ~var_array_t(){
        _destroy_array();
    };

    cmc_universal_type_t operator[](const size_t index) const;
    var_array_t& operator=(const var_array_t& array);
    var_array_t(const var_array_t& array);

    size_t size() const;
    cmc_type get_data_type() const;
    void print_data() const;
    void* get_initial_data_ptr() const;
    
    cmc_universal_type_t sum_over_range(const size_t start_index, const size_t end_index) const;
    cmc_universal_type_t mean_value(const size_t start_index, const size_t end_index) const;
    cmc_universal_type_t mean_value_same_type_wo_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const;
    double mean_value_wo_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const;
    
    bool check_range_fullfills_deviation_threshold_from_value(const size_t start_index, const size_t end_index, const double deviation, const cmc_universal_type_t& value, const cmc_universal_type_t& missing_value) const;
    
    void assign(const size_t index, const cmc_universal_type_t& value);
    /* This function performs a type check of the data inside 'value' */
    void assign_value(const size_t index, const cmc_universal_type_t& value);
    
    void scale(const cmc_universal_type_t& scale_factor);
    void scale_with_missing_vals(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& missing_value);
    
    void add(var_array_t& summand);
    void add_const(const cmc_universal_type_t& constant_summand);
    void add_const_with_missing_vals(const cmc_universal_type_t& offset, const cmc_universal_type_t& missing_value);

    void axpy_scalar(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& add_offset);    
    void axpy_scalar_with_missing_vals(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& add_offset, const cmc_universal_type_t& missing_value);
    
    // Opaque pointer to the internal data structure
    struct var_array* data{nullptr};
private:
    /* Private functions are only considered for internal usage and management of the array */
    void _init_array(const size_t, const cmc_type);
    void _destroy_array();
};

class var_vector_t
{
public:
    var_vector_t(){
        _init_vector();
    };
    var_vector_t(var_vector_t&& vector);
    ~var_vector_t(){
        _destroy_vector();
    };

    size_t size() const;
    size_t capacity() const;
    auto begin() const;
    auto end() const;
    var_vector_t& operator=(var_vector_t&&) = default;	
    var_array_t& operator[](const size_t index);

    void create_and_push_back(const size_t num_elements, const cmc_type data_type);
    /* Ownership of the resources of the pointer are taken (it is assumed, that the resources have been allocated via 'new'. Deallocation will be performed by the destructor of the var_vector) */
    void push_back(var_array_t* array);
    void push_back_placeholder();
    void pop_back();
    void reserve(const size_t num_arrays);
    void clear();
    // Opaque pointer to internal data structure */
    struct var_vector* vec{nullptr};
private:
    void _init_vector();
    void _destroy_vector();
};


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

template<typename T>
auto convert_to_universal_type(const cmc_type dest_type, const T value)
 -> std::enable_if_t<std::is_arithmetic_v<T>, cmc_universal_type_t>
{
    switch (dest_type)
    {
        case CMC_INT32_T:
        {
            return cmc_universal_type_t(static_cast<int32_t>(value));
        }
        break;
        case CMC_FLOAT:
        {
            return cmc_universal_type_t(static_cast<float>(value));
        }
        break;
        case CMC_DOUBLE:
        {
            return cmc_universal_type_t(static_cast<double>(value));
        }
        break;
        case CMC_INT16_T:
        {
            return cmc_universal_type_t(static_cast<int16_t>(value));
        }
        break;
        case CMC_INT64_T:
        {
            return cmc_universal_type_t(static_cast<int64_t>(value));
        }
        break;
        case CMC_UINT64_T:
        {
            return cmc_universal_type_t(static_cast<uint64_t>(value));
        }
        break;
        case CMC_UINT32_T:
        {
            return cmc_universal_type_t(static_cast<uint32_t>(value));
        }
        break;
        case CMC_INT8_T:
        {
            return cmc_universal_type_t(static_cast<int8_t>(value));
        }
        break;
        case CMC_UINT8_T:
        {
            return cmc_universal_type_t(static_cast<uint8_t>(value));
        }
        break;
        case CMC_UINT16_T:
        {
            return cmc_universal_type_t(static_cast<uint16_t>(value));
        }
        break;
        case CMC_BYTE:
        cmc_err_msg("Cannot convert value to type byte.");
        break;
        case CMC_CHAR:
        {
            return cmc_universal_type_t(static_cast<char>(value));
        }
        break;
        default:
            cmc_err_msg("An unknown cmc data type has been supplied.");
    }
} 
#endif /* CMC_CONTAINER_H */
