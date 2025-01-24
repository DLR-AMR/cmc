#ifndef CMC_NC_IO_CONVENTIONS_HXX
#define CMC_NC_IO_CONVENTIONS_HXX

#include <string>

namespace cmc
{

enum CompressionScheme {UndefinedCompression = -1, NoCompression = 0, AdaptiveCoarsening = 0x01, PrefixExtraction = 0x10, AdaptiveCoarseningPlusPrefixExtraction = 0x11, DiffCompression = 0x100};
inline const std::string kCompressionSchemeAttrName = "compression_scheme";
typedef int CompressionSchemeType;


}

#endif /* !CMC_NC_IO_CONVENTIONS_HXX */
