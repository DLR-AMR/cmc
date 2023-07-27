#ifndef CMC_GEO_UTIL_H
#define CMC_GEO_UTIL_H

#include "cmc.h"
#include "cmc_util.h"
#include "cmc_constants_definitions.h"
#include "cmc_container.h"

#define CMC_VAR_NOT_CONSIDERED -1
typedef struct cmc_var* cmc_var_t;
typedef struct cmc_global_coordinate_system* cmc_global_coordinate_system_t;
typedef struct cmc_ref_coordinates* cmc_ref_coordinates_t;

enum cmc_coordinate_type {CMC_COORDINATES_UNDEFINED = 0, CMC_CARTESIAN_COORDINATES, CMC_MORTON_INDEX, CMC_GEO_DOMAIN_DEFINED_BY_BOX};

struct cmc_global_coordinate_system
{
private:
    int reference_count{0};
public:
    int dimensionality{0};
    var_vector_t coords;

    void ref();
    void unref();
    void reset_reference_count();

    cmc_universal_type_t get_coordinate_value_at_dim(const CMC_COORD_IDS coord_id, const size_t reference_position);
    cmc_universal_type_t get_coordinate_value_at_longitude(const size_t reference_position);
    cmc_universal_type_t get_coordinate_value_at_latitude(const size_t reference_position);
    cmc_universal_type_t get_coordinate_value_at_elevation(const size_t reference_position);
};

struct cmc_coordinate
{
private:
    enum current_representation {_CMC_CURRENT_REPR = 0, _CMC_CART_COORD, _CMC_MORTON_IDX};
    current_representation repr{current_representation::_CMC_CURRENT_REPR};
public:
    enum cmc_coordinate_id {_CMC_LON = 0, _CMC_LAT = 1, _CMC_LEV = 2};

    cmc_coordinate(){};
    cmc_coordinate(uint64_t morton_idx)
    : morton_index{morton_index}{};
    cmc_coordinate(std::tuple<uint64_t, uint64_t, uint64_t>&& cart_coordinate)
    : cartesian_coordinate{cart_coordinate}{};
    cmc_coordinate(std::tuple<uint64_t, uint64_t, uint64_t> cart_coordinate)
    : cartesian_coordinate{cart_coordinate}{};

    uint64_t morton_index{0};
    std::tuple<uint64_t, uint64_t, uint64_t> cartesian_coordinate{std::make_tuple(0UL, 0UL, 0UL)};
    
    uint64_t cmc_get_ref_coord(const cmc_coordinate_id _coord_id) const;
    uint64_t cmc_get_cart_coord_by_id(const cmc_coordinate_id _id) const;

    uint64_t get_longitude_ref_coordinate() const;
    uint64_t get_latitude_ref_coordinate() const;
    uint64_t get_elevation_ref_coordinate() const;

    uint64_t get_longitude_coordinate(var_vector_t* global_coordinate_domain) const;
    uint64_t get_latitude_coordinate(var_vector_t* global_coordinate_domain) const;
    uint64_t get_elevation_coordinate(var_vector_t* global_coordinate_domain) const;

    uint64_t get_morton_index() const;
};

struct cmc_ref_coordinates
{
private:
    int reference_count{0};
public:
    cmc_coordinate_type coordinate_representation{cmc_coordinate_type::CMC_COORDINATES_UNDEFINED}; //!< This enum indicates which representation is currently defined in @var coordinates
    //std::vector<cmc_coordinate> coordinates; //!< This vector stores the reference coordinates (compliant to the ordering of the data)
    std::vector<std::tuple<uint64_t, uint64_t, uint64_t>> cartesian_coordinates;
    std::vector<uint64_t> morton_indices;
    int refinement_level_of_morton_indices{-1};
    //var_vector_t* coords{nullptr}; //!< This var vector saves the (process-local) coordinate_values
    std::vector<uint64_t> start_ptr; //!< If a blocked domain is present, this vector indicates the global start of each dimension
    std::vector<uint64_t> count_ptr; //!< If a blocked domain is present, this vector indicates the count in each dimension

    void ref();
    void unref();
};


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

    /********** General Variable Data **********/
    /* Meta information about the variable */
    std::string name{}; //!< The name of the variable
    cmc_universal_type_t missing_value; //!< A vcalue indicating that at this position, no 'real' value exists (this value is just a placeholder for a missing value)
    cmc_universal_type_t add_offset{static_cast<double>(0.0)}; //!< Indicating the offset of the data if it was shifted
    cmc_universal_type_t scale_factor{static_cast<double>(1.0)}; //!< Indicating the scaling of the data 
    bool missing_value_present{false}; //!< Flag indicating whether or not missing values exists within the data
    bool applied_offset_scaling{false}; //!< Flag indicating whether the offset and scaling has been applied to data or not

    /* A pointer to the data array of the variable */
    var_array_t* data{nullptr}; //!< The 'actual' data of the variable
    var_array_t* data_new{nullptr}; //!< A convenience variable which may be used if the data is transformed

    /* Dimensionality of the data */
    int num_dimensions{0}; //!< Dimensionality of the data
    std::vector<uint64_t> dim_lengths; //!< Local dim lengths (in a 't8_..'-setting)

    /* The ordering of the data (e.g. lat x lon x lev) */
    std::vector<int> axis_ordering; //!< The ordering of the data (e.g. 1. Longitude, 2. Latitude ... ) (if the data is given in a linear/cartesian representation)
    /* The corresponding ordering sheme (e.g. linear, z-curve) */
    CMC_DATA_ORDERING_SCHEME data_scheme{CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_SCHEME_UNDEFINED}; //!< An enum indicating if the data is for example given in cartesian coordinates or in Space Filling Curve (i.e. Morton) index
    DATA_LAYOUT data_layout{DATA_LAYOUT::CMC_LAYOUT_UNDEFINED}; //!< An enum for the explicit representation of the data (e.g. CMC_LAT_LON_LEV)

    /* A pointer to the global domain */
    cmc_global_coordinate_system* coordinates{nullptr};

    /* A pointer to a struct collecting the underlying coordinate*s representation of the data */
    cmc_ref_coordinates* ref_coordinates;

    /********** NetCDF Specific Data **********/
    /* ID (for example considered for netCDF files) */
    int var_id{CMC_VAR_NOT_CONSIDERED}; //!> NetCDF varibale id
    int var_type{CMC_VAR_NOT_CONSIDERED}; //!< NetCDF variable type
    std::vector<int> dimension_ids; //!< NetCDF dimension ids

    /* Describing the netCDF data hyperslab */
    std::vector<uint64_t> start_ptr; //!< Coordinate start values of the data
    std::vector<uint64_t> count_ptr; //!< Coordinate length values of the data

    void switch_data();
    void get_data_layout_from_axis_ordering();
    void assign_an_arbitrary_missing_value();
    bool is_equal_to_missing_value(const size_t value) const;
    
};

std::string
get_coord_name(const CMC_COORD_IDS coord_id);

bool
cmc_value_equal_to_missing_value(const cmc_universal_type_t& value, const cmc_universal_type_t& missing_value);

bool
cmc_value_equal_to_zero(const cmc_universal_type_t& value);

double
calculate_two_step_relative_max_deviation(const double previous_max_deviation, const cmc_universal_type_t& previous_mean, const cmc_universal_type_t& current_mean);


#endif /* CMC_GEO_UTIL_H */
