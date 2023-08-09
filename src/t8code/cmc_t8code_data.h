#ifndef CMC_T8CODE_DATA_H
#define CMC_T8CODE_DATA_H
/**
 * @file cmc_t8code_data.h
 * @brief This file collects typedefs for pointers to all structs defined in @file cmc_t8code_data.hxx
 *        Furthermore, it includes all necessary t8code header files needed for the AMR lossy compression.
 */

#include "cmc.h"
#include "utilities/cmc_constants_definitions.h"

/* Include all necessary t8code dependent header files */
#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h>
#include <t8_element_c_interface.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_vtk.h>
#include <t8_forest/t8_forest_partition.h>
#include <p4est.h>
#include <p8est.h>
#endif

/** The maximum number of children in t8code's standard quad/hex scheme */
#define CMC_T8_MAX_CHILDREN_HEX_QUAD 8

/** Typedef of @enum CMC_LOSSY_COMPRESSION_MODE (defined in @file cmc_constants_definitions.h) of possible compression mode (i.e. how the adaptation and interpolation functions should be applied) */
typedef enum CMC_LOSSY_COMPRESSION_MODE CMC_T8_COMPRESSION_MODE;

/** This enum describes all possible compression/coarsening criteria for the AMR lossy compression */
enum CMC_T8_COMPRESSION_CRITERIUM {CMC_CRITERIUM_UNDEFINED = 0, CMC_REL_ERROR_THRESHOLD, CMC_ABS_ERROR_THRESHOLD, CMC_EXCLUDE_AREA, CMC_COMBINED_CRITERION};

/** Forward declarations of all structs (defined in @file cmc_t8code_data.hxx */
typedef struct cmc_t8_geo_data* cmc_t8_geo_data_t;
typedef struct cmc_t8_assets* cmc_t8_assets_t;
typedef struct cmc_t8_var* cmc_t8_var_t;
typedef struct cmc_t8_data* cmc_t8_data_t;
typedef struct cmc_t8_adapt_data* cmc_t8_adapt_data_t;
typedef struct cmc_t8_interpolation_data* cmc_t8_interpolation_data_t;
typedef struct cmc_amr_compression_settings* cmc_amr_compression_settings_t;
typedef struct cmc_t8_forest_interpolate* cmc_t8_forest_interpolate_t;
/** End of forward declarations */

/* A typedef for an interpolation function */
typedef cmc_t8_forest_interpolate_t (*cmc_t8_forest_interpolation_t)(void);

/**
 * @brief This function writes out all variables in @var t8_data on their underlying forests to a VTK file
 * 
 * @note This function is applicable before, during and after the lossy AMR compression
 * @param t8_data A pointer to a @struct cmc_t8_data holding the variables
 * @param file_prefix The supposed name of the VTK file
 */
void
cmc_t8_write_forest_all_vars(cmc_t8_data_t t8_data, const char* file_prefix);

#endif /* CMC_T8CODE_DATA_H */
