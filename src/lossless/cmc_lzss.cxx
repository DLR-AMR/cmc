/**
 * This @file implements a variation of the LZ77 lossless compression algorithm proposed by Abraham Lempel und Jacob Ziv.
 * This variation is similar to the LZSS algorithm; published by James A. Storer und Thomas G. Szymanski in 1982.
 * 
 * The compression algorithm implemented in this file uses a byte aligned version. In order to ensure the byte-alignment, we pack eight values per block and store the blocks subsequently.
 * Within the compression, a byte from the input stream can be either just copied or a sequence of subsequent bytes in the input stream (> 3 bytes) may be encoded using three bytes. At the beginning of each block
 * is one byte, each bit of this byte indicates whether its corresponding value in the block is encoded (=1) or just copied (=0) (e.g. if the first bit is 'zero' -> the first value in the block is just a plain copy of the byte from the input stream; if the second bit is 'one', the following three bytes are an encoding of a compressed sequence, ...)
 * We have a sliding window which indicates the length of the search buffer, in which we want to find a match as long as possible to resemble the current byte in the stream position and as many subsequent bytes in the look ahead buffer as possible.
 *
 * Example:       ....ABCGHDFSHDGDJHFGHSJFHGSDJFDGHJDHG          G           DJHFZDKJTNSD...................
 *                |                                   |          ^           |                             |
 *                |-----------search-buffer-----------|   current position   |------look-ahead-buffer------|
 *         Length:         SLIDING_WINDOW_SIZE                   1               LOOK_AHEAD_BUFFER_SIZE
 *
 * The three-byte-encoding describes reoccurring patterns in the data. The encoding consists of an 'offset-length' pair.
 * We are trying to find a maximum match to a sequence starting at the current position (into the look ahead buffer) which has already occured within the search buffer.
 * It is possible that such a reoccuring sequence starts in the search buffer but continues in the look ahead buffer (in case of repeating patterns), this is handled by the LZSS algorithm as well.
 * If we have a found a maximum match, it is encoded as a three-byte 'offset-length' pair. The offset uses two bytes and indicates the relative offset to the current position at which the match starts.
 * And the third byte is used for the length of the matching sequence.
 * These choices limit the length of the sliding window and the maximum possible match length (described by the macros below).
 *
 * As mentioned before, we need to keep track whether we encode a sequence of bytes or we just copy the byte from the input stream. In order to keep the bytes aligned, we pack eight values per block.
 * Such a block has the following layout:
 * Example:      | 1-Bit-Flags | 8 x {Plain Copy of Byte; Encoded Sequnece of Bytes}
 *        Length:    1 Byte      8 x {      1 Byte      ;          3 Bytes         }   -> At a Maximum the length of block may be 1 + 8 * 3 = 25 Bytes
 *
 * Example 2:    | Flags |
 *               |0010111| Encoding | Encoding | Encoding | Copy | Encoding | Copy | Copy |
 *                   ^
 *                   This flag at the beginning of the block would be followed by the displayed layout
 */

 /* TODOs:
  * Update the maximum possible match length, since an encoding will only appear when the minimum matching length is reached by sequence, we won't have a match length < MINIMUM_MATCHING_LENGTH. Therefore we cann just offset the maximum match length by this value -> MAXIMUM_MATCH_LENGTH = 255 + MINIMUM_MATCH_LENGTH.
  * Currently the sliding window is linearily searched, this may be optimized
  * Curerntly the length is stored as uint16_t (this makes the implementation dependant on the endianness of the system)
  * Common choices for the sliding window sizes are <= 32KB. therefore there is one bit left in our offset within the 'offest-length'-pair, which may be used somehow. (maybe indicating matches which start at very short offsets from the current position < 128. In this case one would not need two bytes for the offset, because one would suffice (one bit to indicate this case and 7 bits for the actual offset). If this is worthy is another question
  * Implement a version of the algorithm which only uses two bytes for the 'offset-length'-pair (for example: 12 Bits for the offset, 4 Bits for the matching length) for comparisons
*/
#include <utility>
#include <tuple>
#include <bitset>
#include <fstream>
#include <string>
#include <string_view>

#include "cmc_lzss.h"
#include "utilities/cmc_log_functions.h"

/** The following mmacros describe some basic properties of the implemented LZSS compression.
 * Therefore, they may not be changed. However, the sliding window size may be changed within the range of a uint16_t, but no matter which size is choosen, the encoding will always use two bytes for the offset in an encodig.
 */
/* The length of a maximum possible match */
#define MAXIMUM_MATCH_LENGTH 255
/* The maximum size of the look ahead buffer */
#define LOOK_AHEAD_BUFFER_SIZE (4 * MAXIMUM_MATCH_LENGTH)
/* The size of the sliding window (therefore, the size of the search buffer) */
#define SLIDING_WINDOW_SIZE 32767
#define FACTOR_FOR_NOW_SWITCH 1000
/* If the file is read in different memeory segments (and not as a whole), this macro defines the amount of sliding windows which will be read from the file */
#define NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK 2
/* The match has to be longer than 3 in order to 'save' storage, because the offset-length pair is encoded using three bytes (but within the comparison is a greater-equal sign (in order to find the maximum match which is nearest to the current position), therefore it is set to four here) */
#define MINIMUM_MATCHING_LENGTH 4
/* Thereare eight values/bytes sequences encoded within 'one block', in order comply with byte alignments */
#define ENCODED_VALUES_PER_BLOCK 8
/* The size needed for encoding a byte sequence (two bytes for an offset and one byte for the length of a match) */
#define ENCODED_VALUE_SIZE (sizeof(uint16_t) + sizeof(uint8_t))
/* The size of one unsigned char (should be one byte */
#define PLAIN_VALUE_BYTE_SIZE (sizeof(unsigned char))
/* The maximum amount of bytes an encoding for eight values may use */
#define MAXIMUM_BYTES_EIGHT_BLOCK_SEQUENCE  (ENCODED_VALUE_SIZE * ENCODED_VALUES_PER_BLOCK * PLAIN_VALUE_BYTE_SIZE + 1)

#if 0
/* The following two functions are currently unused, but will be used in the future. */
/* Overloading for putting out a byte as a sequence of 0's and 1's */
std::ostream& operator<< (std::ostream& os, std::byte b) {
    return os << std::bitset<8>(std::to_integer<int>(b));
}
/**
 * @brief This functions determines during runtime whether the system uses little endian or not.
 * 
 * @return true If little endian byte ordering is used
 * @return false Otherwise, if little endian is not used
 */
static bool
cmc_is_little_endian()
{
    union {
        uint32_t example_int;
        char example_char[4];
    } endianness_test = {0x01020304};

    if (endianness_test.example_char[0] == 1)
    {
        return false;
    } else
    {
        return true;
    }
}
#endif

/**
 * @brief In case we are encoding and current memory segment has come to an end, but the 'eight-block-sequence' is not yet completely filled, 
 * we try to fill it with further encoding of the look ahead buffer. But when even the look ahead buffer is exceeded, we have no chance to fill the 'eight-block-sequence'.
 * In this case, the functions determines how many bytes we have to discard. The amount of bytes is passed back to the calling function of the compression, which then can adjust
 * the sliding_window, the look ahead buffer, etc. for the next iteration.
 * This procdure is used in order to gurantee that the compression result will always be the same, whther the file is read as a whole or in any other arbitraty memory chunks.
 * 
 * @param block [in] The current block the compression tris to fill
 * @param block_count [in] The amount of values which are already encoded in the @var block
 * @return ptrdiff_t The number of bytes we needed to discard from the compression of the current memory segment. (@note This return value will always be a positive integer (but it has too be subtracted later))
 */
static
ptrdiff_t cmc_determine_how_many_bytes_to_discard_while_encoding(unsigned char* block, const size_t block_count)
{
    /* Counter for the amount of bytes which we have to discard */
    ptrdiff_t num_bytes_to_discard{0};

    /* A variable to keep track of the offset in the current encoding block */
    size_t block_offset{1};

    /* Iterate through the whole block */
    for (uint8_t id{0}; id < block_count; ++id)
    {
        /* Check whether the bit set or not */
        if ((((*block) >> id) & static_cast<uint8_t>(1)) != 0)
        {
            /* The bit is set, therefore we have to decode the value.
             * Therefore, the next two bytes describe the offset and the third byte the length */
            /* Obtain the length of the match (we do not need the offset from the match) and update the number of bytes to discard by it */
            num_bytes_to_discard += static_cast<size_t>(*(block + block_offset + 2));

            /* Update the byte offset within the current block */
            block_offset += ENCODED_VALUE_SIZE;

        } else
        {
            /* Increment the counter of bytes to discard by one (because this was just a plain copy) */
            ++num_bytes_to_discard;

            /* Update the byte offset within the current block */
            block_offset += PLAIN_VALUE_BYTE_SIZE;
        }
    }

    /* Return the calculated number */
    return num_bytes_to_discard;
}


/**
 * @brief This function performs a variation of the LZSS lossless compression algorithm.
 *
 * @param [in] in_file_name The path to a file which contents should be compressed
 * @param [in] output_name A name for the output file which will be generated containing the encoded/compressed bytes of the input file. If no output name is supplied (default value = nullptr), than a name will be created
 *
 * @note The algorithm writes an output file containing the compressed sequence. The name of the file is similar to @var in_file_name with an additional suffix '.cmc_lzss' or to the supplied @var output_name with the extension '.cmc_lzss'
 */
void
cmc_lzss_compression(const char* in_file_name, const char* output_name)
{
    cmc_debug_msg("The compression of the file ", in_file_name, " using a variation of the LZSS algorithm starts...");

    /* Variable for stream positions */
    std::streampos fbegin, fend;

    /* Ceate a input file stream and open the given file for input in binary mode */
    std::ifstream file_handle(in_file_name, std::ios::in | std::ios::binary);

    /* Check if the file has been opened */
    if (!(file_handle.is_open()))
    {
        cmc_err_msg("Unable to open the file: ", in_file_name);
    }

    /* Get the beginning of the file */
    fbegin = file_handle.tellg();

    /* Move to the end of the file */
    file_handle.seekg(0, std::ios::end);
    
    /* Get the end of the file */
    fend = file_handle.tellg();

    cmc_debug_msg("The size of the supplied file is: ", (fend - fbegin), " bytes");

    /* Size of the file */
    const auto size = fend - fbegin;

    cmc_assert(size > 0);

    /* Move to the beginning of the file */
    file_handle.seekg(0, std::ios::beg);

    /* Create a vairbale storing the encoded data as well as the amount of bytes needed for the encoding of the stream */
    std::vector<std::tuple<unsigned char*, size_t, ptrdiff_t>> encoded_blocks;

    /* If the file is smaller than this threshold, it is read in as a whole */
    if (size <= 8 * NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE)
    {
        /** Read in the whole data at once **/
        cmc_debug_msg("The file will be read in at once and compressed afterwards.");

        /* Dynamic char array for storing the data of the file */
        char* in_stream = (char*) std::malloc(sizeof(char) * size);

        /* Read the data */
        file_handle.read(in_stream, size);

        /* Close the file */
        file_handle.close();

        /* Compress the data */
        encoded_blocks.push_back(cmc_lzss_compress(reinterpret_cast<unsigned char*>(in_stream), size));

        /* Free the allocated memory */
        std::free(in_stream);
    } else
    {
        cmc_debug_msg("The file will be read in one part after another...\nThe seperate memory segments will be compressed sequentially.\n(Note: This has no impact on the compression result.)");

        /* Dynamic char array for storing the data of the file */
        char* in_stream = (char*) std::malloc(sizeof(char) * NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE);

        /* In this case, we obtain at least two blocks which has to be handled */
        const size_t number_of_encoded_blocks = size / (NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE) + 1;

        /* Allocate space for each encoded memory segment */
        encoded_blocks.reserve(number_of_encoded_blocks);

        /* Dynamic char array for storing the data of the file */
        char* look_ahead_buffer = (char*) std::malloc(sizeof(char) * LOOK_AHEAD_BUFFER_SIZE);
        
        /* Dynamic char array for storing the data of the file */
        char* previous_sliding_window = (char*) std::malloc(sizeof(char) * SLIDING_WINDOW_SIZE);

        /* Variable holding the information how many bytes we are able to read ahead (at most LOOK_AHEAD_BUFFER_SIZE) */
        size_t bytes_to_look_ahead{0};
        /* Counter of the bytes we have encoded so far */
        size_t encoded_bytes{0};

        /* We encode the file stream until we have reached it's end */
        while (encoded_bytes < static_cast<size_t>(size))
        {
            /* In the first iteration, we have no encoded bytes so far and no previous sliding window (to use as a dictionary) */
            if (encoded_bytes == 0)
            {
                /* Move to the beginning of the file */
                file_handle.seekg(0, std::ios::beg);

                /* Read the first memory segment */
                if (size >= NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE)
                {
                    /* Read the current memory segment */
                    file_handle.read(in_stream, NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE);

                    /* Check if we can fully read in a look ahead buffer or a just a little less */
                    bytes_to_look_ahead = (size - NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE >= LOOK_AHEAD_BUFFER_SIZE) ? LOOK_AHEAD_BUFFER_SIZE : size - NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE;

                    /* Read in the look ahead buffer, if there is the opportunity to do so */
                    if (bytes_to_look_ahead > 0)
                    {
                        file_handle.read(look_ahead_buffer, bytes_to_look_ahead);
                    }

                    /* Compress/Encode the current memory segment */
                    encoded_blocks.push_back(cmc_lzss_compress(reinterpret_cast<unsigned char*>(in_stream),
                                                               NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE,
                                                               nullptr,
                                                               (bytes_to_look_ahead > 0 ? reinterpret_cast<unsigned char*>(look_ahead_buffer) : nullptr), bytes_to_look_ahead));
                    
                    /* Update the encoded bytes */
                    encoded_bytes += NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE + std::get<2>(encoded_blocks.back());
                } else
                {
                    /* This case should normally not incide, because if the 'size' is small, the whole file will be read at once and not segmented. But for completeness, this possibility is listed here */
                    file_handle.read(in_stream, size);
                    encoded_blocks.push_back(cmc_lzss_compress(reinterpret_cast<unsigned char*>(in_stream), size, nullptr, nullptr, 0));
                }
            }
            /* In case we would encounter the end of the file with the next memory segment, there is no look ahead buffer anymore for finding possible matches */
            else if (encoded_bytes + NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE >= static_cast<size_t>(size))
            {
                /* In this case the last memory segment will be read (there is no succeding memory segment afterwards, i.e. there is no look ahead buffer) */
                /* Move one sliding window behind the currently encoded bytes and read the 'previous sliding window', the 'current memory segment' and the 'look ahead buffer' from there onwards */
                file_handle.seekg(encoded_bytes - SLIDING_WINDOW_SIZE, std::ios::beg);

                /* Read the previous sliding window */
                file_handle.read(previous_sliding_window, SLIDING_WINDOW_SIZE);

                /* Read the current memory segment */
                file_handle.read(in_stream, size - encoded_bytes);

                /* Compress/encode the last memory segment of the file */
                encoded_blocks.push_back(cmc_lzss_compress(reinterpret_cast<unsigned char*>(in_stream),
                                                           size - encoded_bytes,
                                                           reinterpret_cast<unsigned char*>(previous_sliding_window),
                                                           nullptr, 0));

                /* Assertion that there is no overlapping offset returned, when the stream has come to an end */
                cmc_assert(std::get<2>(encoded_blocks.back()) == 0);

                /* Update the encoded bytes */
                encoded_bytes = size;
            }
            /* In any other case, where a previous sliding window as well as a look ahead buffer is present */
            else
            {
                /* In these cases we have a previous sliding window as well as a look ahead buffer */
                /* Update the position pointer of the file */
                /* Move one sliding window behind the currently encoded bytes and read the 'previous sliding window', the 'current memory segment' and the 'look ahead buffer' from there onwards */
                file_handle.seekg(encoded_bytes - SLIDING_WINDOW_SIZE, std::ios::beg);
                
                /* Read the previous sliding window */
                file_handle.read(previous_sliding_window, SLIDING_WINDOW_SIZE);
                
                /* Read the current memory segment */
                file_handle.read(in_stream, NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE);
                
                /* Check the amount of bytes left in the stream */
                bytes_to_look_ahead = (size - encoded_bytes - NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE >= LOOK_AHEAD_BUFFER_SIZE) ? LOOK_AHEAD_BUFFER_SIZE : size - encoded_bytes - NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE;
                
                /* Read the look ahead buffer. if there is one */
                if (bytes_to_look_ahead > 0)
                {
                    file_handle.read(look_ahead_buffer, bytes_to_look_ahead);
                }
                
                /* Compress/Encode the current memory segment */
                encoded_blocks.push_back(cmc_lzss_compress(reinterpret_cast<unsigned char*>(in_stream),
                                                               NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE,
                                                               reinterpret_cast<unsigned char*>(previous_sliding_window),
                                                               (bytes_to_look_ahead > 0 ? reinterpret_cast<unsigned char*>(look_ahead_buffer) : nullptr), bytes_to_look_ahead));

                /* Update the encoded bytes */
                encoded_bytes += NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE + std::get<2>(encoded_blocks.back());                                                 
            }
        }

        /* Close the file */
        file_handle.close();

        /* Free the look ahead buffer */
        std::free(look_ahead_buffer);

        /* Free the previous sliding window */
        std::free(previous_sliding_window);

        /* Free the allocated memory */
        std::free(in_stream);
    }

    /** Write out the encoded sequnece to an output file **/
    /* Create an output name */
    std::string output_file_name;
    if (output_name != nullptr)
    {
        /* Use the supplied output file name with the custom extension */
        output_file_name = std::string(output_name) + ".cmc_lzss";
    } else
    {
        /* Just use the regular file name plus an extension */
        output_file_name = std::string(in_file_name) + ".cmc_lzss";
    }

    /* Create an output file */
    std::ofstream output_file_handle(output_file_name, std::ios::out | std::ios::binary);

    /* Check if the file has been opened */
    if (!(output_file_handle.is_open()))
    {
        cmc_err_msg("Unable to open the output file: ", output_file_name);
    } 

    /* Counter for the size of the compressed data */
    size_t compressed_data_size{0};

    /* Define an output pointer of type char*, since the read/write operations cannot use unsigned char pointers */
    char* output_ptr{nullptr};

    /* Itearte over all encoded blocks */
    for (auto& [encoded_memory_segment, num_bytes_encoding, overlapping_offsets] : encoded_blocks)
    {
        /* Count the encoded bytes */
        compressed_data_size += num_bytes_encoding;

        /* Cast the pointer for the write operation */
        output_ptr = reinterpret_cast<char*>(encoded_memory_segment);

        /* Write out the memory block to the output file */
        output_file_handle.write(output_ptr, num_bytes_encoding);
    }

    /* Close the output file */
    output_file_handle.close();

    cmc_debug_msg("The output file '", output_file_name, "' containing the compressed data has been written to disk.");
    cmc_debug_msg("The compressed data has a size of ", compressed_data_size, " bytes.");

    /* Free the encoded sequence */
    for (auto& [encoded_memory_segment, num_bytes_encoding, overlapping_offsets] : encoded_blocks)
    {
        /* Free the memory of each encoded segment */
        std::free(encoded_memory_segment);
    }
    
    cmc_debug_msg("The compression using a variation of the LZSS algorithm was successful.");
}

/**
 * @brief This functions performs the decompression of file which was compressed using the cmc-LZSS compression algorithm 
 * 
 * @param [in] in_file_name The file name of the file containing the compressed data
 * @param [in] output_name If not a nullptr, the output file name for the decompressed data
 */
void
cmc_lzss_decompression(const char* in_file_name, const char* output_name)
{
    cmc_debug_msg("The decompression of the file ", in_file_name, " starts...");

    /* Variable for stream positions */
    std::streampos fbegin, fend;

    /* Ceate a input file stream and open the given file for input in binary mode */
    std::ifstream file_handle(in_file_name, std::ios::in | std::ios::binary);

    /* Check if the file has been opened */
    if (!(file_handle.is_open()))
    {
        cmc_err_msg("Unable to open the file: ", in_file_name);
    }

    /* Get the beginning of the file */
    fbegin = file_handle.tellg();

    /* Move to the end of the file */
    file_handle.seekg(0, std::ios::end);

    /* Get the end of the file */
    fend = file_handle.tellg();

    cmc_debug_msg("The size of the file is: ", (fend - fbegin), " bytes");

    /* Size of the file */
    auto size = fend - fbegin;

    /* Move to the beginning of the file */
    file_handle.seekg(0, std::ios::beg);

    /* Create a vairbale storing the encoded data as well as the amount of bytes needed for the encoding of stream */
    std::vector<std::tuple<unsigned char*, size_t, size_t>> decoded_blocks;

    /* If the file is smaller than this threshold, it is read in as a whole */
    if (size <= 16 * NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE)
    {
        /** Read in the whole data at once **/

        /* Dynamic char array for storing the data of the file */
        char* in_stream = (char*) std::malloc(sizeof(char) * size);

        /* Read the data */
        file_handle.read(in_stream, size);

        /* Close the file */
        file_handle.close();

        /* Compress the data */
        decoded_blocks.push_back(cmc_lzss_decompress(reinterpret_cast<unsigned char*>(in_stream), size, nullptr));

        /* Free the allocated memory */
        std::free(in_stream);
    } else
    {
        /* In contrast to the compression, we will not use a look ahead buffer during the decompression. In case an overlapping match (into the look ahead buffer) would occur, we just discard the match and perform it's decompression during the next iteration */

        /* Dynamic char array for storing the data of the file */
        char* in_stream = (char*) std::malloc(sizeof(char) * NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE);

        /* In this case, we obtain at least two blocks which has to be handled */
        const size_t number_of_decoded_blocks = size / (NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE) + 1;

        /* Allocate space for each encoded memory segment */
        decoded_blocks.reserve(number_of_decoded_blocks);

        /* Dynamic char array for storing the data of the file */
        char* previous_sliding_window = (char*) std::malloc(sizeof(char) * SLIDING_WINDOW_SIZE);

        /* Counter of the bytes we have decoded so far */
        size_t decoded_bytes{0};

        /* We encode the file stream until we have reached it's end */
        while (decoded_bytes < static_cast<size_t>(size))
        {
            if (decoded_bytes == 0)
            {
                /* In case the beginning of the file will be decompressed */
                /* Read a memory segment from the file */
                file_handle.read(in_stream, (size >= NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE ? NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE : size));

                /* Just compress the current memory segment without supplying a previous sliding window */
                decoded_blocks.push_back(cmc_lzss_decompress(reinterpret_cast<unsigned char*>(in_stream), (size >= NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE ? NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE : size), nullptr));

                /* Update the counter of the decoded bytes by the bytes we would have desirably decoded minus the bytes we had discarded */
                decoded_bytes += (size >= NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE ? NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE : size) - std::get<2>(decoded_blocks.back());
            } else
            {
                /* In this case it is not the first iteration. Therefore, we can pass a sliding window to the decompression function */
                /* Seek the position in the file which is exactly one sliding window before the current position */
                file_handle.seekg(decoded_bytes, std::ios::beg);

                /* Read in the previous sliding window */
                file_handle.read(previous_sliding_window, SLIDING_WINDOW_SIZE);

                /* Afterwards read the current memory segment from the file which will be decoded */
                file_handle.read(in_stream, (size - decoded_bytes >= NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE ? NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE : size - decoded_bytes));

                /* Compress the current memory segment without supplying a previous sliding window */
                decoded_blocks.push_back(cmc_lzss_decompress(reinterpret_cast<unsigned char*>(in_stream), (size - decoded_bytes >= NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE ? NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE : size - decoded_bytes), nullptr));
            
                /* Update the counter of the decoded bytes by the bytes we would have desirably decoded minus the bytes we had discarded */
                decoded_bytes += (size - decoded_bytes >= NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE ? NUM_SLIDING_WINDOWS_PER_MEMORY_BLOCK * SLIDING_WINDOW_SIZE : size - decoded_bytes) - std::get<2>(decoded_blocks.back());
            }
        }

        /* Free the allocated memory */
        std::free(previous_sliding_window);
        std::free(in_stream);
    }

    /* Write the decompressed values/bytes to an output file */
    std::string output_file_name;

    /* Check if an output name for the resulting file was supplied */
    if (output_name != nullptr)
    {
        /* Save the output name */
        output_file_name = std::string(output_name);
    } else
    {
        /* Save the input file name */
        std::string partly_file_name{std::string(in_file_name)};

        /* Define the extension of the compression */
        constexpr std::string_view compression_extension{".cmc_lzss"};

        /* Search for the extension within the input file name */
        auto extension_ptr = std::search(partly_file_name.begin(), partly_file_name.end(), compression_extension.begin(), compression_extension.end());
        
        /* Check if the input file has the extension we have appended during the compression */
        if (extension_ptr != partly_file_name.end())
        {
            /* If the extension is found in the string, we will erase it for the output file name */
            partly_file_name.erase(extension_ptr, partly_file_name.end());
        } else
        {
            /* The input file name does not have the extension appended during the compression, therefore, print a warning */
            cmc_warn_msg("The file was eventually not compressed using cmc. Therefore, the decompression may fail or lead to undefined behaviour if this is the case.");
        }

        /* Create the output file name */
        output_file_name = "decompr_cmc_lzss_" + partly_file_name;
    }

    /* Create an output file */
    std::ofstream output_file_handle(output_file_name, std::ios::out | std::ios::binary);

    /* Check if the file has been opened */
    if (!(output_file_handle.is_open()))
    {
        cmc_err_msg("Unable to open the output file: ", output_file_name);
    } 

    /* Define an output pointer of type char*, since the read/write operations cannot use unsigned char pointers */
    char* output_ptr{nullptr};

    /* Itearte over all decoded blocks */
    for (auto& [decoded_memory_segment, num_bytes_decoding, discarded_bytes] : decoded_blocks)
    {
        /* Cast the pointer for the write operation */
        output_ptr = reinterpret_cast<char*>(decoded_memory_segment);

        /* Write the memory block to the output file */
        output_file_handle.write(output_ptr, num_bytes_decoding);
    }

    /* Close the output file */
    output_file_handle.close();

    /* Free the decoded sequence */
    for (auto& [decoded_memory_segment, num_bytes_decoding, discarded_bytes] : decoded_blocks)
    {
        /* Free the memory of each decoded segment */
        std::free(decoded_memory_segment);
    }

    cmc_debug_msg("The decompression was successful.");
}


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
 * @note The returned unsigned char pointer has been allocated with 'std::malloc'. The ownership of this memory is given to the calling function. Therefore, this memory segment needs to be deallocated by calling 'std::free()' if not needed anymore
 * @note In case no look ahead buffer is supplied, the third value of the return tuple is always zero.
 */
std::tuple<unsigned char*, size_t, ptrdiff_t>
cmc_lzss_compress(unsigned char* iterator, const size_t num_bytes, unsigned char* previous_sliding_window, unsigned char* look_ahead_buffer, const size_t max_bytes_in_look_ahead_buffer)
{
    /* Allocate memory (enough that even the worst case would be handled without reallocating the memory) (This results in each byte is just plain copied plus one byte per eight values -> resulting in 112.5% memory of 'num_bytes' */
    unsigned char* encoded_sequence = (unsigned char*) std::malloc(sizeof(unsigned char) * std::ceil(1.125 * num_bytes) + MAXIMUM_BYTES_EIGHT_BLOCK_SEQUENCE);
    /* Check if the allocation was succesfull */
    if (encoded_sequence == nullptr)
    {
        cmc_err_msg("Allocation of the block sequence has failed.");
    }
    /* Allocate a block holding a byte indicating which of the following values are encoded or just plain copied */
    unsigned char* block_sequence = (unsigned char*) std::malloc(MAXIMUM_BYTES_EIGHT_BLOCK_SEQUENCE);
    /* Check if the allocation was succesfull */
    if (block_sequence == nullptr)
    {
        cmc_err_msg("Allocation of the block sequence has failed.");
    }

    /* Initialize the one bit flags at the beginning of each block */
    block_sequence[0] = 0;

    /* Counter of values already encoded in the current block; a block is full once it has eight values (consisting of on byte copies or three-byte-encodings) */
    size_t block_count{0};

    /* Counter of the bytes which have been already filled within the current block; Since the first byte resembles the eight one bit flags, the first byte of a new block is already occupied */
    size_t block_occupied_bytes{1};

    /* Counter for current length of the encoded sequence */
    size_t encoded_bytes_count{0};

    /* Counter for the overlapping offset into the look ahead buffer, respectivly counter for the bytes which could not be encoded during this iteration of the function */
    ptrdiff_t overlapping_offset_into_look_ahead_buffer{0};

    /* Variable resembling the current length of the maximum match (during the current sliding window) */
    int max_match_length{MINIMUM_MATCHING_LENGTH};

    /* The length of the match */
    uint8_t current_match_length{1};

    /* Save the initial pointer */
    unsigned char* initial_ptr_backup{iterator};

    /* Assign an initial pointer based on whether a previous sliding window was supplied or not */
    unsigned char* initial_ptr{(previous_sliding_window == nullptr) ? iterator : previous_sliding_window};

    /* A pointer tfor the search buffer */
    unsigned char* search_buffer_ptr{nullptr};

    /* A pointer for the look ahead buffer (the bytes after the current position in the sequence) */
    unsigned char* look_ahead_buffer_ptr{nullptr};

    /* A pair for storing the current offset and length of the maximum match in the sliding window */
    std::pair<uint16_t, uint8_t> encoded_value{std::make_pair(0,0)};

    /* A flag indicating whether a match (worthy of encoding) has been found within the sliding window */
    bool byte_is_encoded{false};

    /* Offset Variables describing the start and length of the 'dictionary'/search buffer */
    size_t offset{0};
    size_t offset_until_current_byte_position{0};

    /* Iterate first until the sliding window is fully within the sequence */
    for (size_t current_byte{0}; current_byte < num_bytes;)
    {
        /* If there was a previous sliding window which is no longer in the buffer */
        if (previous_sliding_window != nullptr && current_byte < SLIDING_WINDOW_SIZE)
        {
            /* Update the offset in the previous sliding window */
            offset = current_byte;
            /* Update the end of the sliding window (located within the current memory segment) */
            offset_until_current_byte_position = SLIDING_WINDOW_SIZE + current_byte;
        } 
        /* If we no longer need the previous sliding window, we switch to the initial pointer of the current sliding window (pointed to by the iterator initially when the function was called) */
        else if (previous_sliding_window != nullptr && (current_byte - current_match_length < SLIDING_WINDOW_SIZE && current_byte >= SLIDING_WINDOW_SIZE)) 
        {
            /* Reset the initial_ptr to the original initial ptr to the start of the current buffer (initially pointed to by iterator when the function was called) */
            initial_ptr = initial_ptr_backup;
            /* Reset the offset to current byte position */
            offset_until_current_byte_position = SLIDING_WINDOW_SIZE;
            /* Reset the offset to the start of the current memory segment */
            offset = 0;
        }
        /* Calculate an offset, when we first excced the sliding window size or if there was a previous sliding window which is no longer in the buffer */
        else if (current_byte > SLIDING_WINDOW_SIZE)
        {
            /* Update the the end of the search buffer */
            offset_until_current_byte_position = current_byte;

            /* Update the offset in the current sliding window */
            offset = current_byte - SLIDING_WINDOW_SIZE;
        } else
        {
            /* In case we have not exceeded the first sliding window of the current memory segment (and there is no previous sliding window supplied), we will iterate at most up to the current byte in the search buffer */
            offset_until_current_byte_position = current_byte;
        }
        
        /* Reset the maximum match length */
        max_match_length = MINIMUM_MATCHING_LENGTH;

        /* Reset the flag indicating whether the value/byte is encoded or just plain copied */
        byte_is_encoded = false;

        /* Skip the last comparisons, because they are not worthy to encode and will be plain copied anyways */
        if (look_ahead_buffer == nullptr && num_bytes - current_byte <= MINIMUM_MATCHING_LENGTH)
        {
            /* We can skip iterations of the last few bytes, because an encoding will not be worthy, so they will be plain copied anyways */
            offset = offset_until_current_byte_position;
        }

        /* Find the longest match within the beginning right from the dictionary pointer */
        for (size_t iter{offset}; iter < offset_until_current_byte_position; ++iter)
        {
            if (previous_sliding_window != nullptr && current_byte <= SLIDING_WINDOW_SIZE && iter +1 > SLIDING_WINDOW_SIZE) // In der vorletzten Ungleichung stand eben groeßer gleich, macht anscheinedn keinen Untersch
            {
                /* Switch to the current memory segment when our match leaves the previous sliding window */
                search_buffer_ptr = iterator - current_byte - SLIDING_WINDOW_SIZE + iter;
            } else
            {
                /* Update the search buffer */
                search_buffer_ptr = initial_ptr + iter;
            }

            /* Check if the current position is the start of a match */
            if (*search_buffer_ptr == *iterator)
            {
                /* Reset the current match length */
                current_match_length = 1;

                /* Save the next pointer positions (after the first match) */
                if (previous_sliding_window != nullptr && current_byte <= SLIDING_WINDOW_SIZE && iter == SLIDING_WINDOW_SIZE)
                {
                    /* Switch to the current memory segment when our match leaves the previous sliding window */
                    search_buffer_ptr = initial_ptr_backup;
                } else 
                {
                    /* Just increment the search_buffer pointer to the byte after the match */
                    ++search_buffer_ptr;
                }

                /* If a look ahead buffer is supplied, check if the match overlappes the current memory block into the look ahead buffer */
                if (look_ahead_buffer != nullptr && current_byte >= num_bytes) 
                {
                    /* If we are exceeding the current memory block, we have to update the pointer to the supplied look ahead buffer (if there is a look ahead buffer succeding this memory block) */
                    look_ahead_buffer_ptr = look_ahead_buffer;
                } else 
                {
                    /* Just increment the look ahead buffer if it is still in a contiguous memory segment */
                    look_ahead_buffer_ptr = iterator + 1;
                }
                
                /* Check if we are allowed to dereference the pointer in the next byte position or if we have alredy exceeded the byte stream */
                if (current_byte + current_match_length + 1 <= num_bytes + max_bytes_in_look_ahead_buffer)
                {
                    /* Iterate until the match ends or the maximum lenght of the match has been reached or if the sequence has come to an end */
                    while (*search_buffer_ptr == *(look_ahead_buffer_ptr))
                    {
                        /* Increment the match length */
                        ++current_match_length;

                        /* Increment the search buffer pointer */
                        if (previous_sliding_window != nullptr && current_byte <= SLIDING_WINDOW_SIZE && iter + current_match_length == SLIDING_WINDOW_SIZE)
                        {
                            /* Switch to the current memory segment when our match leaves the previous sliding window */
                            search_buffer_ptr = initial_ptr_backup;

                        } else 
                        {
                            /* Just increment the search buffer pointer when we are in a contiguous memory segment */
                            ++search_buffer_ptr;
                        }
                        
                        /* Increment the look ahead buffer pointer */
                        if (look_ahead_buffer != nullptr && current_byte + current_match_length -1 > num_bytes)
                        {
                            /* If the match exceeds the current contiguous memory segment, switch to the supplied look ahead buffer */
                            look_ahead_buffer_ptr = look_ahead_buffer + current_match_length - (SLIDING_WINDOW_SIZE - offset);

                            /* Increment the overlapping offset count, if the match still continues */
                            ++overlapping_offset_into_look_ahead_buffer;
                        } else 
                        {
                            /* If the look ahead buffer pointer is in a contiguous memory segment, just increment it */
                            ++look_ahead_buffer_ptr;
                        }

                        /* Check if we are allowed to dereference the pointer in the next iteration or if we have alredy exceeded the byte stream */
                        if ((current_byte + current_match_length + 1) > num_bytes + max_bytes_in_look_ahead_buffer ||
                             current_match_length >= MAXIMUM_MATCH_LENGTH)
                        {
                            break;
                        }
                    }

                    /* Check if the newly found match is longer than the previous maximum match */
                    if (max_match_length <= current_match_length)
                    {
                        /* Update the maximum match length */
                        max_match_length = current_match_length;
                        /* Update the offset and length pair */
                        encoded_value.first = offset_until_current_byte_position - iter;
                        encoded_value.second = current_match_length;
                        /* Set the flag indicating that we have found a match which is worthy to encode */
                        byte_is_encoded = true;
                    }

                    /* Update the counters by the match length that we have found*/
                    iter += current_match_length;
                } else
                {
                    /* The sequence has come to an end */
                    break;
                }
            }
        }

        /* Check if the block of 8 encoded values is filled */
        if (block_count >= 8)
        {
            /* Copy the block to the encoded sequence */
            memcpy((encoded_sequence + encoded_bytes_count), &block_sequence[0], block_occupied_bytes);

            /* Update the number of bytes written */
            encoded_bytes_count += block_occupied_bytes;

            /* Reset the 1-Bit flags byte */
            block_sequence[0] = 0;

            /* Reset the block count */
            block_count = 0;

            /* Reset to one, because the first byte is always ocupied by the 1-Bit flags inidicating an offset-length pair or a byte as is */
            block_occupied_bytes = 1;
        }

        /* If an offset-length pair is given or not */
        if (byte_is_encoded)
        {
            /* Set the one bit flag for the current value (wheter or not the exact byte value is used or a offset-length pair) at the 'end' of the byte */
            /* The bit is set at the index within the byte corresponding to the index of the encoded value within the block */
            block_sequence[0] |= static_cast<unsigned char>(1) << block_count;

            /* Copy the offset and the length of the match */
            //TODO:: Da encoded value ein pair ist, sollte man die zwei Befehle zu einem memcpy befehl zusammen fügen können (this memcpy depends on the endianness)
            memcpy(static_cast<void*>(&block_sequence[block_occupied_bytes]), &(encoded_value.first), sizeof(uint16_t));
            memcpy(static_cast<void*>(&block_sequence[block_occupied_bytes + 2]), &(encoded_value.second), sizeof(uint8_t));

            /* Update the occupied bytes counter */
            block_occupied_bytes += ENCODED_VALUE_SIZE;

            /* Update the pointers and offsets */
            iterator += encoded_value.second;
            current_byte += encoded_value.second;

            /* Reset the encoded value */
            encoded_value.first = 0;
            encoded_value.second = 0;

            /* Reset the flag indicating that a byte sequence has been encoded */
            byte_is_encoded = false;
        } else
        {
            /* Copy the offset and the length of the match */
            memcpy(static_cast<void*>(&block_sequence[block_occupied_bytes]), static_cast<void*>(iterator), PLAIN_VALUE_BYTE_SIZE);
            
            /* Update the occupied bytes */
            block_occupied_bytes += PLAIN_VALUE_BYTE_SIZE;

            /* Move the iterater further */
            ++iterator;

            /* Update the position of the current byte (for the next iteration) */
            ++current_byte;
        }

        /* Increment the block_count */
        ++block_count;
    }

    /* Check if all values have been written to the encoded sequene */
    if (look_ahead_buffer != nullptr && block_count < 8) 
    {
        /* Define the current offset in the look ahead buffer */
        size_t current_look_ahead_byte{static_cast<size_t>(overlapping_offset_into_look_ahead_buffer)};

        /* Define a new iterator for completing the current block */
        unsigned char* extended_iterator{look_ahead_buffer + current_look_ahead_byte};

        /* Variable keeping track of the offset value for the transistion between the the current memory segment and the look ahead buffer */
        size_t iter_offset_correction{0};

        /* Iterate until the current block is complete or the look ahead buffer has come to an end */
        while (block_count < ENCODED_VALUES_PER_BLOCK && current_look_ahead_byte < max_bytes_in_look_ahead_buffer)
        {
            /* We assume that the search buffer will start in the current memroy segment and not in the look ahead buffer (since it should hold SLIDING_WINDOW_SIZE >> overlapping_offset_into_look_hahead_buffer) */
            cmc_assert(current_look_ahead_byte < SLIDING_WINDOW_SIZE);

            /* Therefore, we assign the initial pointer to the current memory segment */
            initial_ptr = initial_ptr_backup;

            /* We assume that the search  buffer starts somewhere in the current memory segment and not in the look ahead buffer */
            offset = num_bytes - SLIDING_WINDOW_SIZE + current_look_ahead_byte;
            offset_until_current_byte_position = num_bytes + current_look_ahead_byte;

            /* Reset the maximum match length */
            max_match_length = MINIMUM_MATCHING_LENGTH;

            /* Reset the flag indicating whether the value/byte is encoded or just plain copied */
            byte_is_encoded = false;

            /* Reset the iter_offset_correction variable */
            iter_offset_correction = 0;

            /* Find the longest match within the beginning right from the dictionary pointer */
            for (size_t iter{offset}; iter < offset_until_current_byte_position; ++iter)
            {
                /* We have to switch to the look ahead buffer, when we excced the current memory segment */
                if (iter == num_bytes + 1)
                {
                    /* Assign the beginning of the look ahead buffer */
                    initial_ptr = look_ahead_buffer;

                    /* Save the correction */
                    iter_offset_correction = iter;
                }

                /* Check if the current position is the start of a match */
                if (*(initial_ptr + iter - iter_offset_correction) == *extended_iterator)
                {
                    /* Reset the current match length */
                    current_match_length = 1;

                    /* Save the next pointer positions (after the first match) */
                    if (iter == num_bytes) 
                    {
                        /* Switch to the current memory segment when our match leaves the previous sliding window */
                        search_buffer_ptr = look_ahead_buffer;
                    } else 
                    {
                        /* Just increment the search_buffer pointer to the byte after the match */
                        search_buffer_ptr = initial_ptr + iter + 1;
                    }

                    /* Update the look ahead buffer */
                    look_ahead_buffer_ptr = extended_iterator + 1;

                    /* Check if we are allowed to dereference the pointer in the next byte position or if we have alredy exceeded the byte stream */
                    if (current_look_ahead_byte + current_match_length + 1 <= max_bytes_in_look_ahead_buffer)
                    {
                        /* Iterate until the match ends or the maximum lenght of the match has been reached or if the sequence has come to an end */
                        while (*search_buffer_ptr == *look_ahead_buffer_ptr)
                        {
                            /* Increment the match length */
                            ++current_match_length;

                            /* Increment the search buffer pointer */
                            if (iter + current_match_length == num_bytes)
                            {
                                /* Switch from the current memory segment to the look ahead buffer */
                                search_buffer_ptr = look_ahead_buffer;
                            } else 
                            {
                                /* Just increment the search buffer pointer when we are in a contiguous memory segment */
                                ++search_buffer_ptr;
                            }
                            
                            /* Increment the look ahead buffer pointer */
                            ++look_ahead_buffer_ptr;

                            /* Check if we are allowed to dereference the pointer in the next iteration or if we have alredy exceeded the byte stream */
                            if ((current_look_ahead_byte + current_match_length + 1) > max_bytes_in_look_ahead_buffer ||
                                current_match_length >= MAXIMUM_MATCH_LENGTH)
                            {
                                break;
                            }
                        }

                        /* Check if the newly found match is longer than the previous maximum match */
                        if (max_match_length <= current_match_length)
                        {
                            /* Update the maximum match length */
                            max_match_length = current_match_length;

                            /* Update the offset and length pair */
                            encoded_value.first = offset_until_current_byte_position - iter;
                            encoded_value.second = current_match_length;

                            /* Set the flag indicating that we have found a match which is worthy to encode */
                            byte_is_encoded = true;
                        }

                        /* Update the counters by the match length that we have found*/
                        iter += current_match_length;

                    } else
                    {
                        cmc_debug_msg("We need to discard some encoded bytes because the look ahead buffer has been exceeded and the current encoding block is not finished.\n(Note: These discarded bytes will be (re-)encoded in the next compression iteration.");

                        /* Get the amount of bytes encoded in this block */
                        ptrdiff_t bytes_encoded_in_block{cmc_determine_how_many_bytes_to_discard_while_encoding(block_sequence, block_count)};
                        
                        /* Free the allocated memory */
                        std::free(block_sequence);
                        
                        /* Return without the last encoding block */
                        return std::make_tuple(encoded_sequence, encoded_bytes_count, static_cast<ptrdiff_t>(current_look_ahead_byte) - bytes_encoded_in_block);
                    }
                }
            }

            /* If an offset-length pair is given or not */
            if (byte_is_encoded)
            {
                /* Set the one bit flag for the current value (wheter or not the exact byte value is used or a offset-length pair) at the 'end' of the byte */
                /* The bit is set at the index within the byte corresponding to the index of the encoded value within the block */
                block_sequence[0] |= static_cast<unsigned char>(1) << block_count;

                /* Copy the offset and the length of the match */
                //TODO:: Da encoded value ein pair ist, sollte man die zwei Befehle zu einem memcpy befehl zusammen fügen können
                memcpy(static_cast<void*>(&block_sequence[block_occupied_bytes]), &(encoded_value.first), sizeof(uint16_t));
                memcpy(static_cast<void*>(&block_sequence[block_occupied_bytes + 2]), &(encoded_value.second), sizeof(uint8_t));

                /* Update the occupied bytes counter */
                block_occupied_bytes += ENCODED_VALUE_SIZE;

                /* Update the pointers and offsets */
                extended_iterator += encoded_value.second;
                current_look_ahead_byte += encoded_value.second;

                /* Reset the encoded value */
                encoded_value.first = 0;
                encoded_value.second = 0;

                /* Reset the encoding flag */
                byte_is_encoded = false;
            } else
            {
                /* Copy the offset and the length of the match */
                memcpy(static_cast<void*>(&block_sequence[block_occupied_bytes]), static_cast<void*>(extended_iterator), PLAIN_VALUE_BYTE_SIZE);

                /* Update the occupied bytes */
                block_occupied_bytes += PLAIN_VALUE_BYTE_SIZE;

                /* Move the iterater further */
                ++extended_iterator;

                /* Update the position of the current byte (for the next iteration) */
                ++current_look_ahead_byte;
            }

            /* Increment the block_count */
            ++block_count;
        }

        /* Check if the block is now full */
        if (block_count == 8) 
        {
            /* If the block is now complete, copy it to the encoded sequence */
            memcpy((encoded_sequence + encoded_bytes_count), &block_sequence[0], block_occupied_bytes);

            /* Update the number of bytes written */
            encoded_bytes_count += block_occupied_bytes;
        } else
        {
            /* Otherwise the look ahead buffer has come to an end and th block is not full.
             * This means, we have to discard the current block and pass a negative offset into the look ahead buffer to the calling function.
             * Such that, this block can be filles during the nect call to 'cmc_lzss_compress(...)
             */
            /* Get the amount of bytes encoded in this block */
            ptrdiff_t bytes_encoded_in_block{cmc_determine_how_many_bytes_to_discard_while_encoding(block_sequence, block_count)};
            
            /* Free the allocated memory */
            std::free(block_sequence);

            /* Return without the last block */
            return std::make_tuple(encoded_sequence, encoded_bytes_count, static_cast<ptrdiff_t>(current_look_ahead_byte) - bytes_encoded_in_block);

        }

        /* Set the overall overlapping bytes into the look ahead buffer */
        overlapping_offset_into_look_ahead_buffer = current_look_ahead_byte;
    } else
    {
        /* If not, write the remaining encoded values in the output sequence */
        memcpy((encoded_sequence + encoded_bytes_count), &block_sequence[0], block_occupied_bytes);
        
        /* Update the number of bytes written */
        encoded_bytes_count += block_occupied_bytes;
    }

    /* Free the allocated memory */
    std::free(block_sequence);

    /* Return the pointer to the memory holding the encoded values as well as an overlapping offset if an look ahead buffer has been supplied to this function */
    return std::make_tuple(encoded_sequence, encoded_bytes_count, overlapping_offset_into_look_ahead_buffer);
}

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
cmc_lzss_decompress(unsigned char* iterator, const size_t num_bytes, unsigned char* previous_sliding_window)
{
    /* Keep track of the allocated memory. This is the initial memory allocation for the decoding sequence */
    size_t currently_allocated_bytes = num_bytes * 8;

    /* Allocate memory for storing the decompressed data */
    unsigned char* decoded_sequence = (unsigned char*) std::malloc(sizeof(unsigned char) * currently_allocated_bytes);
    /* Check if the allocation was succesfull */
    if (decoded_sequence == nullptr)
    {
        cmc_err_msg("Memory allocation for decoding the sequnece has failed.");
    }

    cmc_debug_msg("Decompression of ", num_bytes, " bytes");

    /* A counter which keeps track of the amount of bytes which already have been decoded */
    size_t decoded_bytes_count{0};

    /* If we obtain an encoded value, we use this pair to save the 'offset-length'-pair */
    std::pair<uint16_t, uint8_t> decoded_properties{std::make_pair(0,0)};

    /* A counter for the bytes of a block which already have been decoded. Since the first byte of each block contains the one bit flags, we are initializing this offset counter to one */
    size_t byte_offset_within_block{1};
    
    /* A counter for bytes which will eventually be discarded */
    size_t discarded_bytes_count{0};

    /* Iterate over all bytes in the encoded sequence */
    for (size_t current_byte{0}; current_byte < num_bytes;)
    {
        /* Check if it is possible to exceed the allocated storage with the decoding the next block. If so, reallocate the memory */
        if (decoded_bytes_count + 8 * MAXIMUM_MATCH_LENGTH >= currently_allocated_bytes)
        {
            /* Reallocate the memory */

            /* Double the allocated memory */
            currently_allocated_bytes *= 2;

            /* Reallocate the memory */
            if (unsigned char* tmp_ptr = (unsigned char*) std::realloc(decoded_sequence, currently_allocated_bytes * sizeof(unsigned char)))
            {
                /* Check if the reallocation was succesful */
                cmc_debug_msg("Reallocated memory: ", currently_allocated_bytes, " bytes.");

                /* Assign the reallocated memory */
                decoded_sequence = tmp_ptr;
            } else
            {
                /* If the reallocation has failed */
                cmc_err_msg("Reallocation has been failed.");

                /* Free the allocated memory */
                std::free(decoded_sequence);
            }
        }

        /* Reset the bytes offset from the current iterator position */
        byte_offset_within_block = 1;

        /* Check the first byte of the block describing which of the following eight values is encoded */
        for (uint8_t id{0}; id < ENCODED_VALUES_PER_BLOCK; ++id)
        {
            /* Check if we are excceding the memory */
            if (current_byte + byte_offset_within_block >=num_bytes)
            {
                /* We have reached the end of the stream */
                break;
            }

            /* Check whether the bit set or not */
            if ((((*iterator) >> id) & static_cast<uint8_t>(1)) != 0)
            {
                /* The bit is set, therefore we have to decode the value.
                 * Therefore, the next two bytes describe the offset and the third byte the length
                 */
                /* Obtain the offset fot the match */
                decoded_properties.first = *(reinterpret_cast<uint16_t*>(iterator + byte_offset_within_block));
                /* Obtain the length of the match */
                decoded_properties.second = static_cast<uint8_t>(*(iterator + byte_offset_within_block + 2));

                /* Copy all values to the decoded_sequence */
                /* Check if this match overlaps the current position. In this case we have to be careful with copying */
                if (decoded_properties.first < static_cast<uint16_t>(decoded_properties.second)) 
                {
                    /* Overlapping boundaries occur */
                    /* Copy the non-overlapping bytes over */
                    memcpy(static_cast<void*>(&(decoded_sequence[decoded_bytes_count])), static_cast<void*>(&(decoded_sequence[decoded_bytes_count]) - decoded_properties.first), decoded_properties.first);
                    /* Afterwards, copy the remaining bytes one by one */
                    for (int i{0}; i < decoded_properties.second - decoded_properties.first; ++i)
                    {
                        memcpy(static_cast<void*>(&(decoded_sequence[decoded_bytes_count + decoded_properties.first + i])), static_cast<void*>(&(decoded_sequence[decoded_bytes_count + i])), 1);
                    }
                } else 
                {
                    /* No overlapping boundaries */
                    memcpy(static_cast<void*>(&(decoded_sequence[decoded_bytes_count])), static_cast<void*>(&(decoded_sequence[decoded_bytes_count]) - decoded_properties.first), decoded_properties.second);
                }

                /* Increment the counter of bytes written to the decoded_sequence */
                decoded_bytes_count += decoded_properties.second;

                /* Update the byte offset within the current block */
                byte_offset_within_block += ENCODED_VALUE_SIZE;

            } else
            {
                /* The bit is not set, therefore we can just copy the original value (i.e. this byte was not encoded) */
                memcpy(static_cast<void*>(&(decoded_sequence[decoded_bytes_count])), static_cast<void*>(iterator + byte_offset_within_block), 1);

                /* Increment the counter of bytes written to the decoded_sequence */
                ++decoded_bytes_count;

                /* Update the byte offset within the current block */
                byte_offset_within_block += PLAIN_VALUE_BYTE_SIZE;
            }
        }

        /* Increment the iterator */
        iterator += byte_offset_within_block;
        /* Increment the current count byte */
        current_byte += byte_offset_within_block;
    }

    /* Return the pair consisting of the pointer to the decoded bytes as well as the length of the uncompressed byte sequence */
    return std::make_tuple(decoded_sequence, decoded_bytes_count, discarded_bytes_count);
}
