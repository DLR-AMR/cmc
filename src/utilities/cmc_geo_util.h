#ifndef CMC_GEO_UTIL_H
#define CMC_GEO_UTIL_H

#include "cmc.h"
#include "cmc_util.h"
#include "cmc_constants_definitions.h"
#include "cmc_container.h"

#define CMC_VAR_NOT_CONSIDERED -1
typedef struct cmc_var* cmc_var_t;

struct cmc_var
{
private:

public:
    cmc_var(){};
    cmc_var(std::string _name)
    : name{_name}{};
    cmc_var(const size_t num_elements, const cmc_type type){
        data = new var_array_t{num_elements, type};
    };
    ~cmc_var(){
        if (data != nullptr)
        {
            delete data;
        }
    };
    /* Meta information about the variable */
    cmc_universal_type_t missing_value; 
    cmc_universal_type_t add_offset{static_cast<double>(0.0)};
    cmc_universal_type_t scale_factor{static_cast<double>(1.0)};
    bool missing_value_present{false};
    bool applied_offset_scaling{false};

    /* ID (for example considered for netCDF files) */
    int var_id{CMC_VAR_NOT_CONSIDERED};
    int var_type{CMC_VAR_NOT_CONSIDERED};

    /* Dimensionality of the data and their ids */
    int num_dimensions{0};
    std::vector<int> dimension_ids;
    std::vector<size_t> dim_lengths;

    /* A pointer to the data array of the variable */
    var_array_t* data{nullptr};
    var_array_t* data_new{nullptr};
    
    /* The ordering of the data (e.g. lat x lon x lev) */
    std::vector<int> axis_ordering;
    /* The corresponding ordering sheme (e.g. linear, z-curve) */
    CMC_DATA_ORDERING_SCHEME data_scheme{CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_SCHEME_UNDEFINED};
    DATA_LAYOUT data_layout{DATA_LAYOUT::CMC_LAYOUT_UNDEFINED};

    /* The name of the variable */
    std::string name{};

    /* In a parallel setting, we are saving the start and count data pointer (for a blocked parallel distribution of the variable's data */
    std::vector<uint32_t> start_ptr;
    std::vector<uint32_t> count_ptr;

    void switch_data();
    void get_data_layout_from_axis_ordering();
    void assign_an_arbitrary_missing_value();
    bool is_equal_to_missing_value(const size_t value) const;
};

std::string
get_coord_name(const CMC_COORD_IDS coord_id);

bool
cmc_value_equal_to_missing_value(const cmc_universal_type_t& value, const cmc_universal_type_t& missing_value);

#endif /* CMC_GEO_UTIL_H */