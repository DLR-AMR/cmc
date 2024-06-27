#include "cmc.h"
#include "test/cmc_test.hxx"
#include "utilities/cmc_utilities.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();
    
    cmc::CompressionSettings settings;

    const double abs_max_err = 1.0;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    const cmc::DomainIndex lon_length = 21;
    const cmc::DomainIndex lat_length = 13;
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, lon_length),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, lat_length));

    const double abs_specific_err = 0.75;  
    settings.SetCertainErrorForDomain(cmc::CompressionCriterion::AbsoluteErrorThreshold, abs_specific_err, domain, cmc::kErrorCriterionHoldsForAllVariables);

    cmc::CertainErrorDomain error_domain(cmc::CompressionCriterion::AbsoluteErrorThreshold, 0.5, 
                                         cmc::GeoDomain(cmc::DimensionInterval(cmc::Dimension::Lon, 2, 17),
                                                        cmc::DimensionInterval(cmc::Dimension::Lat, 5, 12)));
    
    /* Per default the settings are always all set for 'all variables' if no specific variable_id has been specified */
    settings.SetCertainErrorForDomain(std::move(error_domain));

    cmc::ExpectTrue(settings.AreTheSettingsValid());

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}
