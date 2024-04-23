/**
 * @file cmc_t8_adapt_track_inaccuracy.hxx
 */

#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"

namespace cmc {











#if 0


InaccuracyTracker::InaccuracyTracker (const TrackingOption tracking_option, int i)
: tracking_option_{tracking_option}{
    switch (tracking_option)
    {
        case TrackingOption::TrackFullInaccuracy:
            deviations = std::make_unique<FullInaccuracyTracker>();
            break;
        case TrackingOption::TrackMinimalWorkingInaccuracy:
            eviations = std::make_unique<MinimalInaccuracyTracker>();
            break;
        default:
            cmc_err_msg("A not supported tracking option for the inaccuracy has been supplied.");
    }
};

#endif



#if 0

double
AbsoluteInaccuracyTracker::ComputeInaccuracy(t8_forest_t current_forest, const int which_tree, const int lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    std::cout << "Compute absolute inacc" << std::endl;

    
    return 0.0;
};

double
RelativeInaccuracyTracker::ComputeInaccuracy()
{
    std::cout << "Compute relative inacc" << std::endl;
    return 0.0;
};

bool
FullInaccuracyTracker::CheckInaccuracy()
{
    std::cout << "Full inacc calls" << std::endl;
    return (inaccuracy_computer_->ComputeInaccuracy() != 0.0 ? true : false);
};






bool
MinimalInaccuracyTracker::CheckInaccuracy()
{
    std::cout << "Minimal inacc calls" << std::endl;
    return (inaccuracy_computer_->ComputeInaccuracy() != 0.0 ? true : false);
};


InaccuracyTracker::InaccuracyTracker (const TrackingOption tracking_option, int i)
: tracking_option_{tracking_option}{
    switch (tracking_option)
    {
        case TrackingOption::TrackFullInaccuracy:
            switch(i)
            {
                case 0:
                    deviations = std::make_unique<FullInaccuracyTracker>(std::make_unique<RelativeInaccuracyTracker>());
                break;
                case 1:
                    deviations = std::make_unique<FullInaccuracyTracker>(std::make_unique<AbsoluteInaccuracyTracker>());
                break;
                default:
                    std::cout << "error" << std::endl;
            }
        break;
        case TrackingOption::TrackMinimalWorkingInaccuracy:
            switch(i)
            {
                case 0:
                    deviations = std::make_unique<MinimalInaccuracyTracker>(std::make_unique<RelativeInaccuracyTracker>());
                break;
                case 1:
                    deviations = std::make_unique<MinimalInaccuracyTracker>(std::make_unique<AbsoluteInaccuracyTracker>());
                break;
                default:
                    std::cout << "error" << std::endl;
            }
        break;
        default:
            std::cout << "error" << std::endl;
    }
};


/**
 * @brief This function calculates (an approximation) of the deviation up to the initial data and decides based on the supplied error threshold whether the element's family shall be coarsened or not
 * 
 * @param adapt_data A pointer to the current @struct cmc_t8_adapt_data holding the status and data of the adaptation
 * @param current_forest The current forest which is built from the adaptation
 * @param lelement_id the first local element id of it's family members
 * @param num_elements The number of elements within this family
 * @return true The element's family complies with the relative error threshold and shall be coarsened 
 * @return false The element's family does not comply with the relative error threshold and therefore is not allowed to be coarsened
 *
 * @note In order to keep track of the deviations it is necessary that if this function returns true, the family strictly has to be coarsened.
 *       In case this functions returns false, the element's family is not allowed to be coarsened. Otherwise this criterion will break!
 */
bool
InaccuracyTracker::CheckInaccuracy(std::vector<Variable>& variables_to_check, std::vector<CertainErrorDomain>& error_domains, ... compression_settings, ... interpolation_struct /* tba */, t8_forest_t current_forest, const int which_tree, const int lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE

    for (auto var_iter = variables_to_check.begin(); var_iter != variables_to_check.end(); ++var_iter)
    {
        /* Calculate the interpolation result if this coarsening would happen */
        const cmc_universal_type_t interpolation_result = interpolation_data.Interpolate(*var_iter, current_forest, which_tree, lelement_id, ts, num_elements, elements);

        /* Compute the inaccuracy resulting from the potential coarsening */
        const double resulting_inaccuracy = deviations_.ComputeInaccuracy(...);

        /* Iterate over all error domains and check whether the deviations complies with the permitted deviation */
        for (auto err_domain_iter = error_domains.begin(); err_domain_iter != error_domains.end(); ++err_domain_iter)
        {
            /* Check whether the element lies geo-spatial domain */
            const bool is_inside_error_domain = IsMeshElementWithinGeoDomain(var_iter->GetAmrMesh(), ..., ts, err_domain_iter->GetGeoDomain());

            /* If so, check whether the inaccuracy complies with the permitted error tolerance */
            if (is_inside_error_domain)
            {
                /* Check the tolerance */
                if (resulting_inaccuracy <= err_domain_iter->allowed_error)
                {
                    /* If the potential coarsening complies, we will store the this resulting deviation (needed for further coarsening steps) */
                    deviations_->StoreResultingInaccuracy(...);
                } else
                {
                    /* Remove the current coarsening informations */
                    deviations_->PopCurrentCoarsening(...); 
                    /* The variable(s) does not comply with the error domain; therefore the error tolerance would be violated and we cannot coarsen the elements */
                    return false;
                }
            }
        }
    }

    /* If the the coarsening complies with each error domain for all variables to be considered, we are able to perform the coarsening step */
    return true;

    #else
    return CMC_ERR;
    #endif
}

#endif

}
