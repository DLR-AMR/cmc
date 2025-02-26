#ifndef CMC_COMPRESSION_SETTINGS_HXX
#define CMC_COMPRESSION_SETTINGS_HXX

/** @file cmc_compression_settings.hxx
 */
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_error_domain.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <vector>

namespace cmc {

enum CompressionMode {CompressionModeUndefined, OneForAll, OneForOne};


class CompressionSettings
{
public:
    void SetGeneralErrorCriterion(const CompressionCriterion criterion, const double permitted_error)
    {
        error_domains_.emplace_back(PermittedError{criterion, permitted_error}, GeneralErrorCriterion);
    }
    void SetErrorDomain(const ErrorDomain& error_domain){error_domains_.push_back(error_domain);};
    void SetErrorDomain(ErrorDomain&& error_domain){error_domains_.push_back(std::move(error_domain));};
private:
    std::vector<ErrorDomain> error_domains_;
};

}

#endif /* !CMC_COMPRESSION_SETTINGS_HXX */
