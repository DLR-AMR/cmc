#ifndef CMC_LZSS_H
#define CMC_LZSS_H

/** This @file cmc_lzss.h holds the function declarations of an implementation of a variation of the LZSS lossless compression algorithm. For further details on the implementation, please take a look in @file src/lossless/cmc_lzss.cxx. */ 

#include "cmc.h"
#include "utilities/cmc_util.h"

/**
 * @brief This function performs a variation of the LZSS lossless compression algorithm.
 *
 * @param [in] in_file_name The path to a file which contents should be compressed
 * @param [in] output_name A name for the output file which will be generated containing the encoded/compressed bytes of the input file. If no output name is supplied (default value = nullptr), than a name will be created
 *
 * @note The algorithm writes an output file containing the compressed sequence. The name of the file is similar to @var in_file_name with an additional suffix '.cmc_lzss' or to the supplied @var output_name with the extension '.cmc_lzss'
 */
void
cmc_lzss_compression(const char* in_file_name, const char* output_name = nullptr);

/**
 * @brief This functions performs the decompression of file which was compressed using the cmc LZSS compression algorithm 
 * 
 * @param [in] in_file_name The file name of the file containing the compressed data
 * @param [in] output_name If not a nullptr, the output file name for the decompressed data  
 */
void
cmc_lzss_decompression(const char* in_file_name, const char* output_name = nullptr);

/**
 * @brief This function performs the actual LZSS lossless compression with optional parameters for a previous sliding window and a look ahead buffer.
 * 
 * @param [in] iterator A pointer to the beginning of an input sequence which will be encoded using a variation of the LZSS algorithm 
 * @param [in] num_bytes The length of the sequence which will be encoded (amount of unsigned chars/bytes)
 * @param [in] previous_sliding_window If not a nullptr, a pointer to an array holding the previous sliding window (@note has to be of length SLIDING_WINDOW_SIZE)
 * @param [in] look_ahead_buffer If not a nullptr, a pointer to an array holding the look ahead buffer (which may be additionally encoded if a match occurs which exceeds the bounds of the sequence @var iterator or if the sequence has come to an end but the current encoding block is not yet filled with 8 encoded values)
 * @param [in] max_bytes_in_look_ahead_buffer The amounts of unsigned chars/bytes in the @var look ahead buffer (if one is supplied) 
 * 
 * @return std::tuple<unsigned char*, size_t, ptrdiff_t> The return value is a tuple of three values. The first is an unsigned char pointer to the encoded sequence. The second value is the amount of bytes (length of) the encoded sequence
 * and the third value is (negative/positive) value resembling an offset into the look_ahead_buffer (e.g. how much bytes of the look ahead buffer already have been encoded or how many bytes of the sequence @var iterator have not been encoded due to an incomplete encoded block (this is only possible if a look ahead buffer has been supplied))
 * 
 * @note The returned unsigned char pointer has been allocated with 'std::malloc'. The ownership of this memory is given to the calling function. Therefore, this memoryy segment needs to be deallocated by calling 'std::free()' if not needed anymore
 * @note In case no look ahead buffer is supplied, the third value of the return tuple is always zero.
 */
std::tuple<unsigned char*, size_t, ptrdiff_t>
cmc_lzss_compress(unsigned char* iterator, const size_t num_bytes, unsigned char* previous_sliding_window = nullptr, unsigned char* look_ahead_buffer = nullptr, const size_t size_look_ahead_buffer = 0);

/**
 * @brief This function performs the actual LZSS lossless compression with optional parameters for a previous sliding window and a look ahead buffer.
 * 
 * @param [in] iterator A pointer to the beginning of an input sequence which will be decoded
 * @param [in] num_bytes The length of the sequence which will be decoded (amount of unsigned chars/bytes)
 * @param [in] previous_sliding_window If not a nullptr, a pointer to an array holding the previous sliding window with decompressed values (@note has to be of length SLIDING_WINDOW_SIZE)
 *
 * @return std::tuple<unsigned char*, size_t, ptrdiff_t> The return value is a tuple of three values. The first is an unsigned char pointer to the decoded sequence. The second value is the amount of bytes (length of) the decoded sequence
 * and the third value is a positive value resembling an offset (e.g. how much bytes of the the current memory segment have not been decoded and therefore had to be discarded. This positive value needs to be subtracted from the offset describing the current position of decompressed bytes when this function is called with the succeding memory segment)
 * 
 * @note The returned unsigned char pointer has been allocated with 'std::malloc'. The ownership of this memory is given to the calling function. Therefore, this memory segment needs to be deallocated by calling 'std::free()' if not needed anymore
 */
std::tuple<unsigned char*, size_t, size_t>
cmc_lzss_decompress(unsigned char* iterator, const size_t num_bytes, unsigned char* previous_sliding_window = nullptr);


#endif /* CMC_LZSS_H */