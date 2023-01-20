#ifndef CMC_T8CODE_DATA_HXX
#define CMC_T8CODE_DATA_HXX
/**
 * @file cmc_t8code_data.hxx
 * @brief This file collects all structs which are used by the 'cmc_t8_...'-functions in order to perform a lossy AMR compression
 */

#include "cmc_t8code_data.h"
#include "mpi/cmc_mpi.h"
#include "utilities/cmc_geo_util.h"
#include "utilities/cmc_container.h"

/**
 * @brief This enum describes the current status of a \struct cmc_t8_data.
 * 
 * The status is checked before current functions are applied in order to ensure the correct ordering of function applications
 * during the lossy compression.
 * 
 */
enum CMC_T8_STATUS {CMC_STATUS_UNDEFINED = 0, SETUP_MESH, SETUP_VARS, VARS_ORDERED, COMPRESSED, DECOMPRESSED};


/**
 * @brief The struct 'cmc_t8_geo_data' saves information about the underlying geo-spatial coordinate system on which the variables' data is defined.
 * 
 * The members of this struct (except @var coords) are used in case all variables are defined ob the same geo-spatial domain.
 * For example, if a 'One for All'-compression-mode is choosen @see @var compression_mode in @struct cmc_t8_data.
 * 
 */
struct cmc_t8_geo_data
{
private:

public:
    cmc_t8_geo_data(){};
    ~cmc_t8_geo_data(){
        if (coords != nullptr)
        {
            delete coords;
        }
    };

    int initial_refinement_lvl{-1}; //!< The (smallest) uniform refinement level of the forest which encloses the geo-spatial domain describes by @var coords.
    int dim{0}; //!< If the data (as well as the forest mesh) is of dimension two or three. 
    int element_anchor_max_lvl{0}; //!< The maximum possible refinement level in case of a 2D or 3D mesh (respective of the @var dim).
    var_vector_t* coords{nullptr}; //!< A pointer to a @struct var_array_t holding all geo-spatial coordinate data (e.g. longitude, latitude, elevation and time coordinates).

    size_t get_coord_length(const CMC_COORD_IDS cmc_coord_id) const;
};

/**
 * @brief The struct 'cmc_t8_assets' holds a @var forest and @var initial_refinement_lvl, on which some (geo-spatial) variable(s) are defined.
 * 
 * @see @struct cmc_t8_var and @struct cmc_t8_data for the use cases fot this @struct cmc_t8_assets.
 */
struct cmc_t8_assets
{
private:

public:
    cmc_t8_assets(){};
    cmc_t8_assets(int initial_ref_lvl)
    : initial_refinement_lvl{initial_ref_lvl}{};
    cmc_t8_assets(t8_forest_t _forest, int initial_ref_lvl)
    : forest{_forest}, initial_refinement_lvl{initial_ref_lvl}{};
    ~cmc_t8_assets(){
        if (forest != nullptr)
        {
            t8_forest_unref(&forest);
        }
    };

    t8_forest_t forest{nullptr}; //!< The forest on which the data of a variable @struct cmc_t8_var is defined
    int initial_refinement_lvl{0}; //!< The initial uniform refinement level of the (smallest) forest which is capable of enclosing the geo-spatial domain of a variable @struct cmc_t8_var
};

/**
 * @brief The struct 'cmc_t8_var' holds a variable (and its data) as well as a pointer to the @struct cmc_t8_assets defining the underlying forest mesh of the variable.
 */
struct cmc_t8_var
{
private:

public:
    cmc_t8_var(){};
    cmc_t8_var(std::string&& name)
    : var{new cmc_var(std::move(name))}{};
    /** \note This constructor tales the full ownership of the 'cmc_var_t' pointer */
    cmc_t8_var(cmc_var_t _var)
    : var{_var} {};
    ~cmc_t8_var(){
        if(assets != nullptr)
        {
            if(assets_allocated)
            {
                delete assets;
            }
        }
        if (var != nullptr)
        {
            delete var;
        }
    };

    cmc_var_t var{nullptr}; //!< A pointer to the variable and it's data
    cmc_t8_assets_t assets{nullptr}; //!< A pointer to the underlying forest
    bool assets_allocated{false}; //!< A flag whether this @struct cmc_t8_var has allocated the @var @struct cmc_t8_assets itself or not

    /** \fn cmc_type get_type() const
     * @return The @enum cmc_type describing the data type of this variable
    */
    cmc_type get_type() const;

    /** \fn size_t get_data_size() const
     * @return The size in bytes of a single datum of this variable
    */
    size_t get_data_size() const;

    /** \fn void* get_initial_data_ptr() const
     * @return A void pointer to the start of a dynamic array holding the actual data of the varibale
    */
    void* get_initial_data_ptr() const;

    /** \fn void* get_initial_data_new_ptr() const
     * @return A void pointer to the start of a second dynamic array which may be allocated and used for calculating adapted data.
    */
    void* get_initial_data_new_ptr() const;

    /** \fn DATA_LAYOUT get_data_layout() const
     * @return An enum resembling the data layout of the variable's data, describing which is the slowest and fastest varying dimension of the data (i.e. the orderinng of the data, for examplae '2D_Lat_Lon')
     */
    DATA_LAYOUT get_data_layout() const;
};

struct cmc_amr_compression_settings
{
public:
    CMC_T8_COMPRESSION_CRITERIUM compression_criterium{CMC_T8_COMPRESSION_CRITERIUM::CMC_CRITERIUM_UNDEFINED};
    double max_err{0.01};
};

/**
 * @brief The struct 'cmc_t8_data' is the general object holding all the variables and all data necessary in order to perform the AMR lossy compression based on t8code.
 */
struct cmc_t8_data
{
private:

public:
    /** Standard destructor to cmc_t8_data */
    cmc_t8_data(){};
    /** Standard constructor with a supplied MPI communicator @var comm.
     * This communicator is used to distribute the variable's data in a parallel case
    */
    cmc_t8_data(MPI_Comm _comm)
    : comm{_comm} {};
    /** Standard destructor deallocating all previously allocated structs and forests */
    ~cmc_t8_data(){
        if (geo_data != nullptr)
        {
            delete geo_data;
        }
        if(assets != nullptr && assets_allocated)
        {
            delete assets;
        }
        if (initial_forest != nullptr)
        {
            t8_forest_unref(&initial_forest);
        }
        for (size_t i{0}; i < vars.size(); ++i)
        {
            delete vars[i];
        }
    };

    cmc_t8_geo_data_t geo_data{nullptr}; //!< A pointer to a @struct cmc_t8_geo_data (@see @struct cmc_t8_geo_data)
    std::vector<cmc_t8_var_t> vars; //!< A vector collecting pointers to all variables which will be compressed (@see @struct cmc_t8_var)

    cmc_t8_assets_t assets{nullptr}; //!< A pointer to a @struct cmc_t8_assets. These general assets are used in case all varibales are defined on the same geo-spatial domain (@see @struct cmc_t8_assets)
    bool assets_allocated{false}; //!< A flag indicating whether these general assets have been allocated (or not if for example some variables are not defined on the same geo-spatial domain)

    t8_forest_t initial_forest{nullptr}; //!< The initial forest enclosing the geo-spatial domain on which ALL variables are defined (only if this is the case, e.g. 'One for All'-compression otherwise it stays a 'nullptr')
    var_vector_t initial_data;  //!< A @struct var_vector_t storing the initial (the uncompressed) data of each variable (this may be used in order to determine an aposteriori error which was introduced by the lossy compression)

    MPI_Comm comm; //!< The communicator to use in a parallel environment 

    CMC_T8_COMPRESSION_MODE compression_mode{CMC_T8_COMPRESSION_UNDEFINED}; //!< An @enum CMC_T8_COMPRESSION_MODE indicating which compression scheme will be used (@see @enum CMC_T8_COMPRESSION_MODE)

    bool variables_are_defined_on_the_same_domain{true}; //!< A flag indicating whether the domain of all variables is the same or if the variables are defined on different geo-spatial domains (e.g. var1 is defined on 'lat x lon'; var2 is defined on 'lon x lev')

    cmc_amr_compression_settings settings{};
};

/**
 * @brief The struct 'cmc_t8_adapt_data' is used by t8code's adapt function in order to supply information during the adaptation step
 */
struct cmc_t8_adapt_data
{
public:
    cmc_t8_adapt_data(){};
    cmc_t8_adapt_data(cmc_t8_data_t _t8_data)
    : t8_data{_t8_data}{};
    cmc_t8_adapt_data(cmc_t8_data_t _t8_data, int _var_id)
    : t8_data{_t8_data}, current_var_id{_var_id}{};
    ~cmc_t8_adapt_data(){};

    cmc_t8_data_t t8_data{nullptr}; //!< A pointer to @struct cmc_t8_data holding all information about the variables and their underlying geo-spatial domains
    int current_var_id{-1}; //!< In case a 'One for One' compression mode is chosen (different meshes are considered for the different variables; this @var current_var_ids helps distinguishing which variable is currently adapted)
    int adapt_step{0}; //!< A counter indicating the amount of adaptation steps that have been applied previously
    
    size_t coarsening_counter{0}; //!< A counter counting the times element families have been coarsened before (during the same adaptation step)
    var_vector_t* adapted_data{nullptr}; //!< A pointer to an universal array which may be used to save data that has been calculated during the adaptation, but may also be used during the interpolation (if for example the same calculations have to be performed)
    
    std::vector<std::unordered_map<t8_locidx_t, t8_locidx_t>> initial_ref_lvl_ids; //!< A hash table saving for each key (the key is id of the coarse element in the current mesh) the id of the first initial element id (in the initial mesh before any coarsening has been introduced) laying within the physical area of the coarse element as a value. This hash table helps to keep track of all initial data points during the coarsening in order to fullfill predefined error threshold criteria.
    int _counter; //!< This counter is used for accesing the hash table @var initial_ref_lvl_ids at the right positions
    int _counter_nxt_lvl; //!< This counter is used for accesing the hash table @var initial_ref_lvl_ids at the right position

    cmc_amr_compression_settings* settings{nullptr}; //!< A pointer to @struct cmc_amr_compression_settings which saves information about the compression criterium
};

/**
 * @brief The struct 'cmc_t8_interpolation_data' is used by t8code's replace function in order to supply information during the interpolation step between two different meshes
 */
struct cmc_t8_interpolation_data
{
public:
    cmc_t8_interpolation_data(){};
    cmc_t8_interpolation_data(cmc_t8_data_t _t8_data)
    : t8_data{_t8_data}{};
    ~cmc_t8_interpolation_data(){};

    cmc_t8_data_t t8_data{nullptr}; //!< A pointer to @struct cmc_t8_data holding all information about the variables and their underlying geo-spatial domains
    int current_var_id{-1}; //!< In case a 'One for One' compression mode is chosen (different meshes are considered for the different variables; this @var current_var_ids helps distinguishing which variable is currently interpolated)
    size_t coarsening_counter{0}; //!< A counter counting the times element families have been coarsened before (during the same adaptation step) (It is used for accessing the @var adapted_data)
    var_vector_t* adapted_data{nullptr}; //!< A pointer to the data calculated and saved during an adaptation step (@see @var adapted_data in @struct cmc_t8_adapt_data)
};

#endif /* CMC_T8CODE_DATA_HXX */