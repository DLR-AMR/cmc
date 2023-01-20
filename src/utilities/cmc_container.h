#ifndef CMC_CONTAINER_H
#define CMC_CONTAINER_H

#include "cmc.h"
#include "utilities/cmc_util.h"

/* Array holding the data_size (in bytes) of each variable data type */ 
inline constexpr std::array<size_t, cmc_type::CMC_NUM_TYPES> cmc_type_to_bytes{sizeof(std::byte), sizeof(int8_t), sizeof(char), sizeof(int16_t), sizeof(int32_t), sizeof(float), sizeof(double), sizeof(uint8_t), sizeof(uint16_t), sizeof(uint32_t), sizeof(int64_t), sizeof(uint64_t)};

/* Define a typedef of variant able to hold severeal data types */
typedef std::variant<std::byte, int8_t, char, int16_t, int32_t, float, double, uint8_t, uint16_t, uint32_t, int64_t, uint64_t> cmc_universal_type_t;

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
    cmc_universal_type_t sum_over_range(const size_t start_index, const size_t end_index) const;
    cmc_universal_type_t mean_value(const size_t start_index, const size_t end_index) const;
    double mean_value_wo_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const;
    cmc_universal_type_t mean_value_same_type_wo_missing_vals(const size_t start_index, const size_t end_index, const cmc_universal_type_t& missing_value) const;
    bool check_range_fullfills_deviation_threshold_from_value(const size_t start_index, const size_t end_index, const double deviation, const cmc_universal_type_t& value, const cmc_universal_type_t& missing_value) const;
    template<typename T> void assign(const size_t index, T value);
    void assign(const size_t index, const cmc_universal_type_t& value);
    /* This function performs a type check of the data inside 'value' */
    void assign_value(const size_t index, const cmc_universal_type_t& value);
    template<typename T> auto axpy(T scale_factor, var_array_t& summand) -> std::enable_if_t<std::is_arithmetic_v<T>, void>;
    template<typename T> auto axpy_scalar(T scale_factor, T constant_summand) -> std::enable_if_t<std::is_arithmetic_v<T>, void>;
    template<typename T> auto scale(T scale_factor) -> std::enable_if_t<std::is_arithmetic_v<T>, void>;
    void scale_with_missing_vals(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& missing_value);
    template<typename T> auto add_const(T constant_summand) -> std::enable_if_t<std::is_arithmetic_v<T>, void>;
    void add_const_with_missing_vals(const cmc_universal_type_t& offset, const cmc_universal_type_t& missing_value);
    void axpy_scalar_with_missing_vals(const cmc_universal_type_t& scale_factor, const cmc_universal_type_t& add_offset, const cmc_universal_type_t& missing_value);
    void add(var_array_t& summand);
    
    void print_data() const;
    void* get_initial_data_ptr() const;

    // Opaque pointer to the internal data structure
    struct var_array* data{nullptr};
private:
    /* Private functions are only considered for internal usage and management of the array */
    void _init_array(const size_t, const cmc_type);
    void _destroy_array();
    void _intern_axpy(void*, var_array_t&);
    void _intern_axpy_scalar(void*, void*);
    void _intern_scale(void*);
    void _intern_assign(const size_t, void*);
    void _intern_add_const(void*);
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


#endif /* CMC_CONTAINER_H */