#ifndef CMC_CONSTANTS_DEFINITIONS_H
#define CMC_CONSTANTS_DEFINITIONS_H

#ifdef __cplusplus
#include <climits>
#else
#include <limits.h>
#endif
#define CMC_MACRO_EXPANSION(x) #x
#define CMC_MACRO_EXPANSION2(x) CMC_MACRO_EXPANSION(x)
#define CMC_FILE_LOCATION __FILE__ ": " CMC_MACRO_EXPANSION2(__LINE__)

/* Do not change these definitions */
#define CMC_NUM_COORD_IDS 4 /* latitude, longitude, leverage and time */
#define CMC_NUM_GEO_COORDS 3 /* latitude, longitude and leverage */
#define CMC_COORDINATE_NOT_CONSIDERED -1

#define CMC_APPLY_ZCURVE_TO_ALL_VARS -1
#define CMC_APPLY_OFFSET_AND_SCALING_TO_ALL_VARS -1
#define CMC_APPLY_TRANSFORMATION_3D_TO_2D_TO_ALL_VARS -1

enum CMC_COORD_IDS {CMC_LAT, CMC_LON, CMC_LEV, CMC_TIME}; // CMC_COORD_FIRST = CMC_LAT, CMC_COORD_LAST = CMC_TIME};
enum CMC_DATA_ORDERING_SCHEME {CMC_GEO_DATA_SCHEME_UNDEFINED, CMC_GEO_DATA_MESSY, CMC_GEO_DATA_LINEAR, CMC_GEO_DATA_Z_CURVE};
enum DATA_LAYOUT {CMC_LAYOUT_UNDEFINED = 0, CMC_2D_LAT_LON, CMC_2D_LON_LAT, CMC_2D_LAT_LEV, CMC_2D_LEV_LAT, CMC_2D_LON_LEV, CMC_2D_LEV_LON, _INTERN_ID_END_2D_START_3D,
                  CMC_3D_LAT_LON_LEV, CMC_3D_LAT_LEV_LON, CMC_3D_LEV_LAT_LON, CMC_3D_LEV_LON_LAT, CMC_3D_LON_LEV_LAT, CMC_3D_LON_LAT_LEV};
/* Enumeration of all fully supported data types for 'cmc_var_vector_t' */
enum cmc_type {CMC_UNDEFINED = -1, CMC_BYTE, CMC_INT8_T, CMC_CHAR, CMC_INT16_T, CMC_INT32_T, CMC_FLOAT, CMC_DOUBLE, CMC_UINT8_T, CMC_UINT16_T, CMC_UINT32_T, CMC_INT64_T, CMC_UINT64_T, CMC_NUM_TYPES};

/* Enumeration of possible compression mode (i.e. how the adaptation and interpolation functions should be applied) */
/* If CMC is built with Fortran, we enforce the enum to be a 4-byte integer, in order to be compliant with Fortran's 'C_INT' */
//enum CMC_LOSSY_COMPRESSION_MODE {CMC_T8_COMPRESSION_UNDEFINED = -1, ONE_FOR_ALL_2D, ONE_FOR_ONE_2D, GROUPED_2D, ONE_FOR_ALL_3D, ONE_FOR_ONE_3D, GROUPED_3D
enum CMC_LOSSY_COMPRESSION_MODE {CMC_T8_COMPRESSION_UNDEFINED = -1, ONE_FOR_ALL, ONE_FOR_ONE, GROUPED
#ifdef CMC_ENABLE_FORTRAN
, CMC_DUMMY_INT32_VAL = INT_MAX
#endif
};


#endif /* CMC_CONSTANTS_DEFINITIONS_H */