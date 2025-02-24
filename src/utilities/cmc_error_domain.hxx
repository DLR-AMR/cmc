#ifndef CMC_ERROR_DOMAIN_HXX
#define CMC_ERROR_DOMAIN_HXX

#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <vector>

namespace cmc
{

struct PermittedError
{
    PermittedError() = delete;
    PermittedError(const CompressionCriterion etype, const double permitted_error)
    : criterion{etype}, error{permitted_error}{};

    const CompressionCriterion criterion{CompressionCriterion::CriterionUndefined};
    const double error{0.0};
};

struct ErrorCompliance
{
    ErrorCompliance() = delete;
    ErrorCompliance(const bool is_error_threshold_fulfilled, const double max_error)
    : is_error_threshold_satisfied{is_error_threshold_fulfilled}, max_introduced_error{max_error}{};

    const bool is_error_threshold_satisfied;
    const double max_introduced_error;
};


}

#endif /* !CMC_ERROR_DOMAIN_HXX */
