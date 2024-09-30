#ifndef CMC_COMPRESSION_CLASS_HXX
#define CMC_COMPRESSION_CLASS_HXX

#include <string>

namespace cmc
{

namespace lossy
{

class Compressor
{
public:
    Compressor() = default;
    virtual ~Compressor() = default;

    virtual void Setup() = 0;
    virtual void Compress() = 0;
    virtual void GatherCompressedData() = 0;
};

class Decompressor
{
public:
    Decompressor() = default;
    virtual ~Decompressor() = default;

    virtual void Decompress(const std::string& file_name) = 0;
};

}

}

#endif /* !CMC_COMPRESSION_CLASS_HXX */
