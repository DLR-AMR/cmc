#ifndef CMC_AMR_COMPRESSOR_H
#define CMC_AMR_COMPRESSOR_H
/** @file cmc_amr_compressor.h
 * All functions in order to perform the compression/decompression of the AMR lossy compressor lie within this file
 */

#include "cmc.h"
#include "utilities/cmc_constants_definitions.h"
#include "mpi/cmc_mpi.h"
//#include "netcdf/cmc_netcdf.h"
#include "messy/cmc_messy.h"
#include "t8code/cmc_t8_adapt_callbacks.h"
#include "t8code/cmc_t8_replace_callbacks.h"


/** Opaque pointer declaration of @struct cmc_amr_data */
typedef struct cmc_amr_data* cmc_amr_data_t;

/** Typedef for @enum CMC_LOSSY_COMPRESSION_MODE describing the compression mode (definition in @file cmc_constants_definitions.h) */
typedef enum CMC_LOSSY_COMPRESSION_MODE CMC_AMR_COMPRESSION_MODE;

/** Macro indicating that all variables should be written in a netCDF file (and not just a single variable; @see \fn void cmc_amr_write_netcdf_file(cmc_amr_data_t amr_data, const char* path, const int var_id) ) */
#define CMC_AMR_WRITE_ALL_VARS_TO_NETCDF -1

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif

/**
 * @brief This function creates a @struct cmc_amr_data based on the data collected previously in @var nc_data from a netCDF file
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param comm The MPI communicator to use in a parallel environment
 * @return A pointer the newly created @struct cmc_amr_data
 */
cmc_amr_data_t
cmc_create_amr_compression_data(cmc_nc_data_t nc_data, const MPI_Comm comm);

#if 0
//TODO: update messy fucntion
cmc_amr_data_t
cmc_create_amr_compression_data_messy(cmc_messy_data_t messy_data, const MPI_Comm comm);
#endif

/**
 * @brief This function converts a 3D variable into several 2D variables by splitting up one coordinate dimension.
 * For example, a 3D variable defined on "longitude x latitude x elevation" can be uncoupled in #eleveation "longitude x latitude" variables
 * 
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param var_id Either a specific id, or the macro 'CMC_APPLY_TRANSFORMATION_3D_TO_2D_TO_ALL_VARS' in order to convert all 3D variables to 2D variables
 * @param preferred_data_layot An @enum DATA_LAYOUT setting the data layout the new 2D variables will have (e.g. CMC_2D_LAT_LON @see @enum DATA_LAYOUT)
 */
void
cmc_amr_pre_setup_split_3D_variable(cmc_amr_data_t amr_data, const int var_id, const enum DATA_LAYOUT preferred_data_layout);

/**
 * @brief This function sets an error threshold compression criterium for the lossy AMR compressor. The data of all variables will be compressed compliant to the this threshold.
 * Precisly, the introduced relative data loss during the compression will not excced @var maximum_error_tolerance.
 * 
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param maximum_rel_error_tolerance Sets the maximum relative data loss for the compression in percent. The value should be a decimal value, e.g. 0.01 (-> equals a maximum data loss of 1%)
 */
void
cmc_amr_pre_setup_set_compression_criterium_relative_error_threshold(cmc_amr_data_t amr_data, const double maximum_rel_error_tolerance);

/**
 * @brief This function sets an error threshold compression criterium for the lossy AMR compressor. The data of all variables will be compressed compliant to the this threshold.
 * Precisly, the introduced data absolute loss during the compression will not excced @var maximum_error_tolerance.
 * 
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param maximum_abs_error_tolerance Sets the maximum absolute data loss for the compression
 */
void
cmc_amr_pre_setup_set_compression_criterium_absolute_error_threshold(cmc_amr_data_t amr_data, const double maximum_abs_error_tolerance);

/**
 * @brief Resets all the previous supplied compression settings which were set via the 'cmc_amr_pre_setup_...()'-functions
 *
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain which will be compressed and information about the coordinate system
 */
void
cmc_amr_pre_setup_reset_compression_settings(cmc_amr_data_t amr_data);

/**
 * @brief This function sets whether or not the initial data is saved or lost during the compression.
 * 
 * @param flag_keep_initial_data A Value unequal to zero indicates that the data should be kept during the compression
 *
 * @note By default, the initial data is will vanish after the first adaptation/compression step and therefore is not kept.
 * @note A benefit from keeping the data would be for example if a compression is directly followed by a decompression, we could afterwards verify the introduced data inaccuracy if the data was kept
 */
void
cmc_amr_pre_setup_set_flag_in_order_to_keep_the_initial_data(cmc_amr_data_t amr_data, const int flag_keep_initial_data);

/**
 * @brief This function sets an error threshold compression criterium for the lossy AMR compressor. The data of all variables will be compressed compliant to the this threshold.
 * Precisly, the introduced data loss during the compression will not excced @var maximum_error_tolerance.
 * 
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param coord_id The corrdinate which should be restriceted (e.g. longitude: CMC_COORD_IDS::CMC_LON)
 * @param start_value An actual coordinate value specifiying a start (in the given coordinate dimension @var coord_id) of an area which will be excluded during the compression 
 * @param end_value An actual coordinate value specifiying a end (in the given coordinate dimension @var coord_id) of an area which will be excluded during the compression 
 */
#ifdef __cplusplus
void
cmc_amr_pre_setup_set_compression_criterium_exclude_area(cmc_amr_data_t amr_data, const enum CMC_COORD_IDS coord_id, const cmc_universal_type_t& start_value, const cmc_universal_type_t& end_value);
#else
/* In order to keep a C-compliant interface, the start and end value needs to default to double values */
void
cmc_amr_pre_setup_set_compression_criterium_exclude_area(cmc_amr_data_t amr_data, const enum CMC_COORD_IDS coord_id, const double start_value, const double end_value);
#endif

/**
 * @brief This function sets the variables' data up fpr the lossy AMR compression based on the given @var compression_mode.
 * The @var compression_mode could be 'One for One' or 'One for All' in 2D or 3D (@see @enum CMC_AMR_COMPRESSION_MODE).
 * The main difference between the two approaches is, that in case of a 'One for One' compression, each variable is compressed independetly with the supplied compression criterium.
 * Ergo, each varibale is defined on it's own mesh (One mesh for one variable).
 * In case of a 'One for All' compression, all variables are defined on the same mesh and will be compressed simultaneously, e.g. the compression criterium must be fulfilled for each variable.
 * Therefore, the compresseion of the variables is dependent on one another (One mesh for all variables).
 * 
 * @note By default, an error threshold criterium is used (if no else has been set via a 'cmc_amr_pre_setup_...'-function).
 * @note Any 'cmc_amr_pre_setup_...' have to be called before \fn void cmc_amr_setup_compression(cmc_amr_data_t amr_data, CMC_AMR_COMPRESSION_MODE compression_mode)
 *
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param compression_mode Sets the compression mode for the lossy AMR compressor
 */
void
cmc_amr_setup_compression(cmc_amr_data_t amr_data, CMC_AMR_COMPRESSION_MODE compression_mode);


#ifdef __cplusplus
/**
 * @brief This function performs the lossy AMR compression based on the data which was aquired by a call to one of the 'cmc_create_amr_compression_data...'-functions and the settings which have been previously set via \fn void cmc_amr_setup_compression(cmc_amr_data_t amr_data, CMC_AMR_COMPRESSION_MODE compression_mode)
 * 
 * @note It is possible to pass nullptr's to the adapt- and interpolation-functions. This results in a default case using error threshold adaptation and chooses the arithmetic mean as interpolation function.
 * @note In a C++ environment this is the default case. However, to have a C-compliant interface, thsi function is defined below without default values.

 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param adapt_function A function pointer describing the adaptation/compression which will be used (@see @file cmc_t8_adapt_callbacks.h for all opportunities)
 * @param interpolation_function A function pointer describing which interpolation will be used in order to map the fine data on a coarser domain (@see @file cmc_t8_replace_callbacks.h for all opportunities)
 */
void
cmc_amr_compress(cmc_amr_data_t amr_data, const t8_forest_adapt_t adapt_function = nullptr, cmc_t8_forest_interpolation_t interpolation_function = nullptr);
#else
/**
 * @brief This function performs the lossy AMR compression based on the data which was aquired by a call to one of the 'cmc_create_amr_compression_data...'-functions and the settings which have been previously set via \fn void cmc_amr_setup_compression(cmc_amr_data_t amr_data, CMC_AMR_COMPRESSION_MODE compression_mode)
 * 
 * @note It is possible to pass nullptr's to the adapt- and interpolation-functions. This results in a default case using error threshold adaptation and chooses the arithmetic mean as interpolation function.
 *
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param adapt_function A function pointer describing the adaptation/compression which will be used (@see @file cmc_t8_adapt_callbacks.h for all opportunities)
 * @param interpolation_function A function pointer describing which interpolation will be used in order to map the fine data on a coarser domain (@see @file cmc_t8_replace_callbacks.h for all opportunities)
 */
void
cmc_amr_compress(cmc_amr_data_t amr_data, const t8_forest_adapt_t adapt_function, cmc_t8_forest_interpolation_t interpolation_function);
#endif

/**
 * @brief This function decompresses the (previously) lossy AMR compressed data
 * 
 * @note Currently, it is only possible to decompress data after the initial data has been compressed (because it is currently not possible to store the compressed data)
 *
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param adapt_function A function pointer describing the adaptation/compression which will be used (@see @file cmc_t8_adapt_callbacks.h for all opportunities)
 * @param interpolation_function A function pointer describing which interpolation will be used in order to map the fine data on a coarser domain (@see @file cmc_t8_replace_callbacks.h for all opportunities)
 */
void
cmc_amr_decompress(cmc_amr_data_t amr_data);

void
cmc_amr_write_netcdf_file(cmc_amr_data_t amr_data, const char* path, const int var_id);

/**
 * @brief This function writes a vtk-file with the variables (initial data, compressed data or decompressed data) as well as the underlying meshes to one or severeal vtk-files
 * 
 * @param amr_data A pointer to a @struct cmc_amr_data holding the variables defined on a geo-spatial domain and information about the coordinate system
 * @param file_prefix The prefix the vtk-files will be given
 */
void
cmc_amr_write_vtk_file(cmc_amr_data_t amr_data, const char* file_prefix);

/**
 * @brief This function destroys a @struct cmc_amr_data and deallocates it's members
 * 
 * @param amr_data A pointer to a @struct cmc_amr_data which will be destroyed/deallocated
 */
void
cmc_amr_destroy(cmc_amr_data_t amr_data);


#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif


#endif /* CMC_AMR_COMPRESSOR_H */
