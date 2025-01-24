#ifndef CMC_T8_ADAPT_CALLBACKS_HXX
#define CMC_T8_ADAPT_CALLBACKS_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_span.hxx"
#include "utilities/cmc_prefix.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h>
#include <t8_element_c_interface.h>
#include <t8_forest/t8_forest_iterate.h>
#include <p4est.h>
#include <p8est.h>
#endif

#include <vector>

namespace cmc
{

#ifdef CMC_WITH_T8CODE

/* Helper functions for return values during the t8code adaptation call */
constexpr t8_locidx_t kCoarsenElements = -1;
constexpr t8_locidx_t kRefineElement = 1;
constexpr t8_locidx_t kLeaveElementUnchanged = 0;

struct AdaptDataInitialMesh;

struct RefinementBits;

struct DecompressionRefinementBits;

t8_locidx_t
RefineToInitialMesh (t8_forest_t forest,
                     t8_forest_t forest_from,
                     t8_locidx_t which_tree,
                     t8_locidx_t lelement_id,
                     t8_eclass_scheme_c * ts,
                     const int is_family,
                     const int num_elements,
                     t8_element_t * elements[]);

t8_locidx_t
PerformAdaptiveCoarseningOneForOne (t8_forest_t forest,
                                    t8_forest_t forest_from,
                                    int which_tree,
                                    int lelement_id,
                                    t8_eclass_scheme_c * ts,
                                    const int is_family,
                                    const int num_elements,
                                    t8_element_t * elements[]);

t8_locidx_t
PerformAdaptiveCoarseningOneForAll (t8_forest_t forest,
                                    t8_forest_t forest_from,
                                    int which_tree,
                                    int lelement_id,
                                    t8_eclass_scheme_c * ts,
                                    const int is_family,
                                    const int num_elements,
                                    t8_element_t * elements[]);

t8_locidx_t
PerformAdaptiveCoarseningOneForOneRegardingInitialData (t8_forest_t forest,
                                    [[maybe_unused]] t8_forest_t forest_from,
                                    [[maybe_unused]] int which_tree,
                                    int lelement_id_,
                                    t8_eclass_scheme_c * ts,
                                    const int is_family,
                                    const int num_elements,
                                    t8_element_t * elements[]);

t8_locidx_t
FindRefinementBits (t8_forest_t forest,
                    t8_forest_t forest_from,
                    int which_tree,
                    int lelement_id,
                    t8_eclass_scheme_c * ts,
                    const int is_family,
                    const int num_elements,
                    t8_element_t * elements[]);

t8_locidx_t
ApplyRefinementBits (t8_forest_t forest,
                     t8_forest_t forest_from,
                     int which_tree,
                     int lelement_id,
                     t8_eclass_scheme_c * ts,
                     const int is_family,
                     const int num_elements,
                     t8_element_t * elements[]);

t8_locidx_t
ExtractCommonPrefixes (t8_forest_t forest,
                       [[maybe_unused]] t8_forest_t forest_from,
                       [[maybe_unused]] int which_tree,
                       int lelement_id,
                       [[maybe_unused]] t8_eclass_scheme_c * ts,
                       const int is_family,
                       const int num_elements,
                       [[maybe_unused]] t8_element_t * elements[]);

t8_locidx_t
DecompressPrefixEncoding (t8_forest_t forest,
                          t8_forest_t forest_from,
                          int which_tree,
                          int lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[]);


t8_locidx_t
DecompressSuffixEncoding (t8_forest_t forest,
                          t8_forest_t forest_from,
                          int which_tree,
                          int lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[]);

t8_locidx_t
DecompressPlainSuffixEncoding (t8_forest_t forest,
                          t8_forest_t forest_from,
                          int which_tree,
                          int lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[]);

t8_locidx_t
ExtractMeanAndLeaveDiffs (t8_forest_t forest,
                       [[maybe_unused]] t8_forest_t forest_from,
                       [[maybe_unused]] int which_tree,
                       int lelement_id,
                       [[maybe_unused]] t8_eclass_scheme_c * ts,
                       const int is_family,
                       const int num_elements,
                       [[maybe_unused]] t8_element_t * elements[]);

t8_locidx_t
ExtractMean (t8_forest_t forest,
                       [[maybe_unused]] t8_forest_t forest_from,
                       [[maybe_unused]] int which_tree,
                       int lelement_id,
                       [[maybe_unused]] t8_eclass_scheme_c * ts,
                       const int is_family,
                       const int num_elements,
                       [[maybe_unused]] t8_element_t * elements[]);

t8_locidx_t
BuildFittingPyramid (t8_forest_t forest,
                       [[maybe_unused]] t8_forest_t forest_from,
                       [[maybe_unused]] int which_tree,
                       int lelement_id,
                       [[maybe_unused]] t8_eclass_scheme_c * ts,
                       const int is_family,
                       const int num_elements,
                       [[maybe_unused]] t8_element_t * elements[]);

t8_locidx_t
DecompressDiffEncoding (t8_forest_t forest,
                          t8_forest_t forest_from,
                          int which_tree,
                          int lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[]);

struct AdaptDataInitialMesh
{
public:
    AdaptDataInitialMesh() = delete;
    AdaptDataInitialMesh(const GeoDomain& domain, const int initial_refinement_lvl, const DataLayout layout)
    : global_domain{domain}, initial_refinement_level{initial_refinement_lvl}, initial_layout{layout}{};

    const GeoDomain& global_domain;
    const int initial_refinement_level;
    const DataLayout initial_layout;
};

struct RefinementBits
{
public:
    RefinementBits() = delete;
    RefinementBits(const int size_hint)
    {
        refinement_indicator.reserve(size_hint / CHAR_BIT + 1);
        refinement_indicator.emplace_back(0);
    }

    int current_bit_position{0};
    std::vector<uint8_t> refinement_indicator;
};

struct DecompressionRefinementBits
{
public:
    DecompressionRefinementBits() = delete;
    DecompressionRefinementBits(const VectorView<uint8_t>& encoded_refinement_for_current_level)
    : encoded_refinements{encoded_refinement_for_current_level}{};
    DecompressionRefinementBits(VectorView<uint8_t>&& encoded_refinement_for_current_level)
    : encoded_refinements{std::move(encoded_refinement_for_current_level)}{};

    int byte_position{0};
    int bit_position{0};
    VectorView<uint8_t> encoded_refinements;  
};


//TODO: Remove below
#if 0
template<int N>
t8_locidx_t
DecodePrefixEGU (t8_forest_t forest,
              [[maybe_unused]] t8_forest_t forest_from,
              [[maybe_unused]] int which_tree,
              int lelement_id,
              [[maybe_unused]] t8_eclass_scheme_c * ts,
              const int is_family,
              const int num_elements,
              [[maybe_unused]] t8_element_t * elements[])
{
    /* Get the adapt data from the forest */
    DecompressionPrefixAdaptDataEGU<N>* adapt_data = static_cast<DecompressionPrefixAdaptDataEGU<N>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    const int initial_refinement_level = 3;//

    /* Possibly, adapt the position indicators */
    if (adapt_data->prefix_indication_bit_position == CHAR_BIT)
    {
        ++(adapt_data->prefix_indication_byte_position);
        adapt_data->prefix_indication_bit_position = 0;
    }

    if (t8_element_level(ts, elements[0]) >= initial_refinement_level)
    {
        /* Copy the current element, if there is no prefix and the elements are outside of the geo domain */
        adapt_data->KeepCompressionValue(lelement_id);
        
        /* Set the position to the succeding bit for the next element */
        ++(adapt_data->prefix_indication_bit_position);
        return kLeaveElementUnchanged;
    }

    const uint8_t current_byte = adapt_data->prefix_indication_bits[adapt_data->prefix_indication_byte_position];
    
    /* IF the bit is not set, but we have not reached the initial refinement level, we need to refine and copy num_children times */
    if (CheckIfBitIsSet(current_byte, adapt_data->prefix_indication_bit_position))
    {
        /* If it is set, there is a prefix */
        if (IsMeshElementWithinGlobalDomain(elements[0], ts, adapt_data->domain, initial_refinement_level, adapt_data->current_layout))
        {
            const auto [decoded_prefix, num_bits] = GetNextPrefix(adapt_data->prefix_lengths, adapt_data->prefix_length_byte_position, adapt_data->prefix_length_bit_position,
                                                                  adapt_data->prefix_encodings,  adapt_data->prefix_encoding_byte_position, adapt_data->prefix_encoding_bit_position);
            //uint32_t temp_pref{0};
            //std::memcpy(&temp_pref, decoded_prefix.data(), decoded_prefix.size());
            //cmc_debug_msg("For elem_id: ", lelement_id, ", LengthPref: ", num_bits, " Bitset: ", std::bitset<8*4>(temp_pref));
            adapt_data->ApplyPrefixForRange(decoded_prefix, num_bits, lelement_id, t8_element_num_children(ts, elements[0]));

            ++(adapt_data->prefix_indication_bit_position);
            return kRefineElement;
        } else
        {
             const auto [decoded_prefix, num_bits] = GetNextPrefix(adapt_data->prefix_lengths, adapt_data->prefix_length_byte_position, adapt_data->prefix_length_bit_position,
                                                                  adapt_data->prefix_encodings,  adapt_data->prefix_encoding_byte_position, adapt_data->prefix_encoding_bit_position);

            adapt_data->ApplyPrefixSingle(decoded_prefix, num_bits, lelement_id);
            ++(adapt_data->prefix_indication_bit_position);
            return kLeaveElementUnchanged;
        }
    } else
    {
        /* IF the the element is inside the domain, it is a real decompression value, otherwise it might be a missing value */
        if (IsMeshElementWithinGlobalDomain(elements[0], ts, adapt_data->domain, initial_refinement_level, adapt_data->current_layout))
        {
            //cmc_debug_msg("For elem_id: ", lelement_id, ", Bit is not set");
            adapt_data->KeepAndCopyCompressionValue(lelement_id, t8_element_num_children(ts, elements[0]));
            ++(adapt_data->prefix_indication_bit_position);
            return kRefineElement;
        } else
        {
            adapt_data->KeepCompressionValue(lelement_id);
        
            /* Set the position to the succeding bit for the next element */
            ++(adapt_data->prefix_indication_bit_position);
            return kLeaveElementUnchanged;
        }
    }
}


template<int N>
t8_locidx_t
DecodeSuffixEGU (t8_forest_t forest,
              [[maybe_unused]] t8_forest_t forest_from,
              [[maybe_unused]] int which_tree,
              int lelement_id,
              [[maybe_unused]] t8_eclass_scheme_c * ts,
              const int is_family,
              const int num_elements,
              [[maybe_unused]] t8_element_t * elements[])
{
    /* Get the adapt data from the forest */
    DecompressionPrefixAdaptDataEGU<N>* adapt_data = static_cast<DecompressionPrefixAdaptDataEGU<N>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    const int initial_refinement_level = 3;//

    /* Possibly, adapt the position indicators */
    if (adapt_data->prefix_indication_bit_position == CHAR_BIT)
    {
        ++(adapt_data->prefix_indication_byte_position);
        adapt_data->prefix_indication_bit_position = 0;
    }

    const uint8_t current_byte = adapt_data->prefix_indication_bits[adapt_data->prefix_indication_byte_position];
    
    /* IF the bit is not set, but we have not reached the initial refinement level, we need to refine and copy num_children times */
    if (CheckIfBitIsSet(current_byte, adapt_data->prefix_indication_bit_position))
    {
        const auto [decoded_prefix, num_bits] = GetNextPrefix(adapt_data->prefix_lengths, adapt_data->prefix_length_byte_position, adapt_data->prefix_length_bit_position,
                                                              adapt_data->prefix_encodings,  adapt_data->prefix_encoding_byte_position, adapt_data->prefix_encoding_bit_position);

        cmc_debug_msg("Decoded prefix hat num_bits: ", num_bits);
        adapt_data->ApplyPrefixSingle(decoded_prefix, num_bits, lelement_id);
        ++(adapt_data->prefix_indication_bit_position);
        return kLeaveElementUnchanged;
    } else
    {
        adapt_data->KeepCompressionValue(lelement_id);
    
        /* Set the position to the succeding bit for the next element */
        ++(adapt_data->prefix_indication_bit_position);
        return kLeaveElementUnchanged;
    }
}
#endif

#endif /* CMC_WITH_T8CODE */

}

#endif /* !CMC_T8_ADAPT_CALLBACKS_HXX */
