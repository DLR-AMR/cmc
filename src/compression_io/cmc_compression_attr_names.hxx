#ifndef CMC_COMPRESSION_ATTR_NAMES_HXX
#define CMC_COMPRESSION_ATTR_NAMES_HXX

#include "utilities/cmc_compression_schema.hxx"

#include <string>

namespace cmc::compression_io
{
    const std::string id_attr("id");
    const std::string data_type_attr("data_type");
    const std::string mesh_id_attr("mesh_id");
    const std::string compression_schema_attr("compression_schema");
    const std::string num_compression_iterations_attr("num_compression_iterations");

    inline std::string GenerateMeshVariableName(const int id)
    {
        return std::string("mesh_") + std::to_string(id);
    }
}


#endif /* !CMC_COMPRESSION_ATTR_NAMES_HXX */
