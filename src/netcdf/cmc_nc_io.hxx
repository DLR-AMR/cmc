#ifndef CMC_NC_IO_HXX
#define CMC_NC_IO_HXX

#include "cmc.h"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_geo_domain.hxx"

#ifdef CMC_WITH_NETCDF
#include <netcdf.h>
#endif
#ifdef CMC_WIT_NETCDF_PAR
#include "netcdf_par.h"
#endif

#include <string>
#include <vector>
#include <variant>

namespace cmc
{

class NcAttribute;
template<typename T>
class NcSpecificVariable;
class NcVariable;


class NcAttribute
{
public:
    NcAttribute() = default;
    ~NcAttribute() = default;

    NcAttribute(const NcAttribute& other) = default;
    NcAttribute& operator=(const NcAttribute& other) = default;
    NcAttribute(NcAttribute&& other) = default;
    NcAttribute& operator=(NcAttribute&& other) = default;
private:
    std::string name_;
    CmcUniversalType value_;
};

template<typename T>
class NcSpecificVariable
{
public:
    NcSpecificVariable() = default;
    ~NcSpecificVariable() = default;

    NcSpecificVariable(const NcSpecificVariable& other) = default;
    NcSpecificVariable& operator=(const NcSpecificVariable& other) = default;
    NcSpecificVariable(NcSpecificVariable&& other) = default;
    NcSpecificVariable& operator=(NcSpecificVariable&& other) = default;
private:
    std::string name_;
    std::vector<T> data_;
    T missing_value_;
    DataLayout layout_;
    DataFormat format_;
    std::vector<Hyperslab> hyperslabs_;
    std::vector<LinearIndex> linear_indices_;
    GeoDomain global_domain_;
};

using NcGeneralVariable = std::variant<NcSpecificVariable<int8_t>, NcSpecificVariable<char>, NcSpecificVariable<int16_t>, NcSpecificVariable<int32_t>, NcSpecificVariable<float>, NcSpecificVariable<double>,
                                       NcSpecificVariable<uint8_t>, NcSpecificVariable<uint16_t>, NcSpecificVariable<uint32_t>, NcSpecificVariable<int64_t>, NcSpecificVariable<uint64_t>>;


class NcVariable
{
public:
    NcVariable() = default;
    ~NcVariable() = default;

    NcVariable(const NcVariable& other) = default;
    NcVariable& operator=(const NcVariable& other) = default;
    NcVariable(NcVariable&& other) = default;
    NcVariable& operator=(NcVariable&& other) = default;
    
private:
    NcGeneralVariable variable_;
    std::vector<NcAttribute> attributes_;
};

}

#endif /* !CMC_NC_IO_HXX */
