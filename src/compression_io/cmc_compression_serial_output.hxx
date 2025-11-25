#ifndef CMC_COMPRESSION_SERIAL_OUTPUT_HXX
#define CMC_COMPRESSION_SERIAL_OUTPUT_HXX

#include "cmc.hxx"
#include "utilities/cmc_iface_patch_compression_variable.hxx"

#include <string>
#include <vector>

namespace cmc::compression_io::serial
{

struct SerialOutputVariable
{
    SerialOutputVariable(std::string&& name, std::vector<uint8_t>&& byte_stream)
    : name_(std::move(name)), byte_stream_(std::move(byte_stream)) {}
    
    std::string name_;
    std::vector<uint8_t> byte_stream_;
};

class Writer
{
public:
    Writer(const std::string& file_prefix)
    : file_prefix_{file_prefix}{};

    template<typename T> void SetVariable(cmc::IPatchCompressionVariable<T>* variable);

    void Write();

private:
    template<typename T> SerialOutputVariable SetDataVariable(cmc::IPatchCompressionVariable<T>* variable);

    std::string file_prefix_;
    std::vector<SerialOutputVariable> variables_;
};

}

#include "compression_io/cmc_patch_byte_compression_serial_output.txx"

#endif /* !CMC_COMPRESSION_SERIAL_OUTPUT_HXX */
