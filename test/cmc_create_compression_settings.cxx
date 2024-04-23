#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "test/cmc_test.hxx"
#include "cmc.h"

int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();
    
    cmc::CompressionSettings settings;

    const double abs_max_err = 1.0;
    const double rel_max_err = 0.05;

    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);
    settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    const int var_id = 0;
    const cmc::Dimension dimension_to_split = cmc::Dimension::Lev;

    cmc::SplitVariable split(var_id, dimension_to_split);
    settings.SplitVariableByDimension(split);

    split.variable_id = 1;
    settings.SplitVariableByDimension(std::move(split));

    cmc::ExpectEQ(settings.AreThereVariablesToSplit());

    cmc::ExpectEQ(settings.AreTheSettingsValid());

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}
