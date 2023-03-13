/** This @file example/cmc_lzss.cxx implements an example for using the lossless LZSS compression and decompression routines */

#include "cmc.h"
#include "lossless/cmc_lzss.h"
#include "utilities/cmc_log_functions.h"

int
main()
{
    /* Initialize cmc */
    cmc_initialize();

    #if 0
    /* Compress the contents of a file, e.g. called 'example.txt' */
    /* The second argument is the name the output file will have, additionally suffixed by '.cmc_lzss' */
    cmc_lzss_compression("lorem_ipsum.txt", "example_compression");

    /* If we want to decompress the contents of a file, we are specifying the name of this file and pass the desired output file name as a second argument (it will be additionally prefixed).
     * The decompressed data will reside in the output file */
    cmc_lzss_decompression("lorem_ipsum.txt.cmc_lzss", "example_text.txt");
    #else
    cmc_global_msg("Currently, there is no example supplied for the LZSS compression/decompression.");
    #endif

    /* Finalize cmc */
    cmc_finalize();

    return 0;
}