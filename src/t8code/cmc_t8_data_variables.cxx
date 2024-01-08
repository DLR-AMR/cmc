#include "t8code/cmc_t8_data_variables.hxx"

namespace cmc
{

inline std::string
Variable::GetName() const
{
    return name_;
}

bool
Variable::IsValidForCompression() const
{
    //TODO: implement
    return true;
}

inline AmrMesh
Variable::GetAmrMesh() const
{
    return mesh_;
}

}