#ifndef CMC_T8_ADAPT_CALLBACKS_H
#define CMC_T8_ADAPT_CALLBACKS_H

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_span.hxx"
#include "utilities/cmc_prefix.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
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

template<int N>
class PrefixDecompressionAdaptData;

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
FindPrefixBitsEGU (t8_forest_t forest,
                    [[maybe_unused]] t8_forest_t forest_from,
                    [[maybe_unused]] int which_tree,
                    int lelement_id,
                    [[maybe_unused]] t8_eclass_scheme_c * ts,
                    const int is_family,
                    const int num_elements,
                    [[maybe_unused]] t8_element_t * elements[]);


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

template<int N>
class PrefixDecompressionAdaptData
{
public:
    PrefixDecompressionAdaptData() = delete;
    PrefixDecompressionAdaptData(const VectorView<uint8_t>& span_prefix_indication_bits, const VectorView<uint8_t>& span_prefix_encoding_bytes)
    : prefix_indication_bits{span_prefix_indication_bits}, prefix_encoding_bytes{span_prefix_encoding_bytes} {};

    bool IsDeCompressionProgressing() const;
    void InitializeCompressionIteration(const int size_hint);
    void FinalizeCompressionIteration();

    void ApplyPrefixForRange(const std::vector<uint8_t>& encoded_prefix, const int num_prefix_bits, const int start_index, const int num_elemnts);
    void KeepCompressionValue(const int index);

    const VectorView<uint8_t>& prefix_indication_bits;
    const VectorView<uint8_t>& prefix_encoding_bytes;
    std::vector<CompressionValue<N>> reconstructed_variable;
    std::vector<CompressionValue<N>> reconstructed_variable_new;
    int byte_position{0};
    int bit_position{0};
    int prefix_encoding_bytes_used{0};
    int count_adaptation_step_{0};
};


template<int N>
bool
PrefixDecompressionAdaptData<N>::IsDeCompressionProgressing() const
{
    return (byte_position < prefix_indication_bits.size() ? true : false);
}

template<int N>
void 
PrefixDecompressionAdaptData<N>::InitializeCompressionIteration(const int size_hint)
{
    reconstructed_variable_new.reserve(size_hint);
}

template<int N>
void 
PrefixDecompressionAdaptData<N>::FinalizeCompressionIteration()
{
    std::swap(reconstructed_variable, reconstructed_variable_new);
    reconstructed_variable_new.clear();
    ++count_adaptation_step_;
}

template<int N>
void
PrefixDecompressionAdaptData<N>::ApplyPrefixForRange(const std::vector<uint8_t>& decoded_prefix, const int num_prefix_bits, const int start_index, const int num_elements)
{
    if (count_adaptation_step_ > 0)
    {
        /* Apply the prefix to the compression values of the range */
        reconstructed_variable_new.push_back(CompressionValue(decoded_prefix, num_prefix_bits));
    } else
    {
        /* Create a new compression value from the encoded prefix */
        CompressionValue<N> prefixed_value = reconstructed_variable[start_index];
        prefixed_value.ApplyPrefix(decoded_prefix, num_prefix_bits);
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, prefixed_value);
    }
}



template<int N>
class DecompressionPrefixAdaptData
{
public:
    DecompressionPrefixAdaptData() = delete;
    DecompressionPrefixAdaptData(const VectorView<uint8_t>& span_prefix_indication_bits, const VectorView<uint8_t>& span_prefix_lengths, const VectorView<uint8_t>& span_prefix_encodings)
    : prefix_indication_bits{span_prefix_indication_bits}, prefix_lengths{span_prefix_lengths}, prefix_encodings{span_prefix_encodings} {};

    bool IsDeCompressionProgressing() const;
    void InitializeCompressionIteration(const int size_hint);
    void FinalizeCompressionIteration();

    void ApplyPrefixForRange(const std::vector<uint8_t>& encoded_prefix, const int num_prefix_bits, const int start_index, const int num_elemnts);
    void KeepCompressionValue(const int index);
    void KeepAndCopyCompressionValue(const int start_index, const int num_elemnts);

    const VectorView<uint8_t>& prefix_indication_bits;
    const VectorView<uint8_t>& prefix_lengths;
    const VectorView<uint8_t>& prefix_encodings;

    std::vector<CompressionValue<N>> reconstructed_variable;
    std::vector<CompressionValue<N>> reconstructed_variable_new;

    int prefix_indication_byte_position{0};
    int prefix_indication_bit_position{0};

    int prefix_encoding_byte_position{0};
    int prefix_encoding_bit_position{0};

    int prefix_length_byte_position{0};
    int prefix_length_bit_position{CHAR_BIT};

    int count_adaptation_step_{0};
    GeoDomain domain;
    DataLayout current_layout{DataLayout::LayoutUndefined};
};

template<int N>
bool
DecompressionPrefixAdaptData<N>::IsDeCompressionProgressing() const
{
    return (prefix_indication_byte_position < prefix_indication_bits.size() ? true : false);
}

template<int N>
void 
DecompressionPrefixAdaptData<N>::InitializeCompressionIteration(const int size_hint)
{
    reconstructed_variable_new.reserve(size_hint);
}

template<int N>
void 
DecompressionPrefixAdaptData<N>::FinalizeCompressionIteration()
{
    std::swap(reconstructed_variable, reconstructed_variable_new);
    reconstructed_variable_new.clear();
    prefix_indication_bit_position = 0;
    ++prefix_indication_byte_position;
    prefix_length_bit_position = CHAR_BIT;
    ++prefix_length_byte_position;
    prefix_encoding_bit_position = 0;
    ++prefix_encoding_byte_position;
    ++count_adaptation_step_;
}

template<int N>
void
DecompressionPrefixAdaptData<N>::ApplyPrefixForRange(const std::vector<uint8_t>& decoded_prefix, const int num_prefix_bits, const int start_index, const int num_elements)
{
   // cmc_debug_msg("in apply prefix for range");
    if (count_adaptation_step_ > 0)
    {
        /* Apply the prefix to the compression values of the range */
        CompressionValue<N> prefixed_value = reconstructed_variable[start_index];
        prefixed_value.ApplyPrefix(decoded_prefix, num_prefix_bits);
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, prefixed_value);
    } else
    {
        /* Create a new compression value from the encoded prefix */
        const CompressionValue<N> prefixed_value(decoded_prefix, num_prefix_bits);
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, prefixed_value);
    }
}

template<int N>
void
DecompressionPrefixAdaptData<N>::KeepCompressionValue(const int index)
{
    if (count_adaptation_step_ > 0)
    {
        reconstructed_variable_new.push_back(reconstructed_variable[index]);
    } else
    {
        reconstructed_variable_new.push_back(CompressionValue<N>());
    }
}

template<int N>
void
DecompressionPrefixAdaptData<N>::KeepAndCopyCompressionValue(const int start_index, const int num_elements)
{
   // cmc_debug_msg("in apply prefix for range");
    if (count_adaptation_step_ > 0)
    {
        /* Set the compression value several times wihtin the new vector */
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, reconstructed_variable[start_index]);
    } else
    {
        /* Create a new compression value from the encoded prefix */
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, CompressionValue<N>());
    }
}



template<int N>
class DecompressionPrefixAdaptDataEGU
{
public:
    DecompressionPrefixAdaptDataEGU() = delete;
    DecompressionPrefixAdaptDataEGU(const VectorView<uint8_t>& span_prefix_indication_bits, const VectorView<uint8_t>& span_prefix_lengths, const VectorView<uint8_t>& span_prefix_encodings)
    : prefix_indication_bits{span_prefix_indication_bits}, prefix_lengths{span_prefix_lengths}, prefix_encodings{span_prefix_encodings} {};

    DecompressionPrefixAdaptDataEGU(const VectorView<uint8_t>& span_prefix_indication_bits, const VectorView<uint8_t>& span_prefix_lengths, const VectorView<uint8_t>& span_prefix_encodings, const std::vector<CompressionValue<N>>& pref_reconstructed_variable)
    : prefix_indication_bits{span_prefix_indication_bits}, prefix_lengths{span_prefix_lengths}, prefix_encodings{span_prefix_encodings}, reconstructed_variable{pref_reconstructed_variable}, count_adaptation_step_{1} {};

    bool IsDeCompressionProgressing() const;
    void InitializeCompressionIteration(const int size_hint);
    void FinalizeCompressionIteration();

    void ApplyPrefixForRange(const std::vector<uint8_t>& encoded_prefix, const int num_prefix_bits, const int start_index, const int num_elemnts);
    void KeepCompressionValue(const int index);
    void KeepAndCopyCompressionValue(const int start_index, const int num_elemnts);
    void ApplyPrefixSingle(const std::vector<uint8_t>& decoded_prefix, const int num_prefix_bits, const int start_index);
    const VectorView<uint8_t>& prefix_indication_bits;
    const VectorView<uint8_t>& prefix_lengths;
    const VectorView<uint8_t>& prefix_encodings;

    std::vector<CompressionValue<N>> reconstructed_variable;
    std::vector<CompressionValue<N>> reconstructed_variable_new;

    int prefix_indication_byte_position{0};
    int prefix_indication_bit_position{0};

    int prefix_encoding_byte_position{0};
    int prefix_encoding_bit_position{0};

    int prefix_length_byte_position{0};
    int prefix_length_bit_position{CHAR_BIT};

    int count_adaptation_step_{0};
    GeoDomain domain;
    DataLayout current_layout{DataLayout::LayoutUndefined};
};

template<int N>
bool
DecompressionPrefixAdaptDataEGU<N>::IsDeCompressionProgressing() const
{
    return (prefix_indication_byte_position < prefix_indication_bits.size() ? true : false);
}

template<int N>
void 
DecompressionPrefixAdaptDataEGU<N>::InitializeCompressionIteration(const int size_hint)
{
    reconstructed_variable_new.reserve(size_hint);
}

template<int N>
void 
DecompressionPrefixAdaptDataEGU<N>::FinalizeCompressionIteration()
{
    reconstructed_variable = reconstructed_variable_new;
    reconstructed_variable_new.clear();

    prefix_indication_bit_position = 0;
    ++prefix_indication_byte_position;

    prefix_length_bit_position = CHAR_BIT;
    ++prefix_length_byte_position;

    prefix_encoding_bit_position = 0;
    ++prefix_encoding_byte_position;

    ++count_adaptation_step_;
}

template<int N>
void
DecompressionPrefixAdaptDataEGU<N>::ApplyPrefixForRange(const std::vector<uint8_t>& decoded_prefix, const int num_prefix_bits, const int start_index, const int num_elements)
{
   // cmc_debug_msg("in apply prefix for range");
    if (count_adaptation_step_ > 0)
    {
        /* Apply the prefix to the compression values of the range */
        CompressionValue<N> prefixed_value = reconstructed_variable[start_index];
        prefixed_value.ApplyPrefix(decoded_prefix, num_prefix_bits);
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, prefixed_value);
    } else
    {
        /* Create a new compression value from the encoded prefix */
        const CompressionValue<N> prefixed_value(decoded_prefix, num_prefix_bits);
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, prefixed_value);
    }
}

template<int N>
void
DecompressionPrefixAdaptDataEGU<N>::ApplyPrefixSingle(const std::vector<uint8_t>& decoded_prefix, const int num_prefix_bits, const int start_index)
{
   // cmc_debug_msg("in apply prefix for range");
    if (count_adaptation_step_ > 0)
    {
        /* Apply the prefix to the compression values of the range */
        CompressionValue<N> prefixed_value = reconstructed_variable[start_index];

        prefixed_value.ApplyPrefix(decoded_prefix, num_prefix_bits);

        reconstructed_variable_new.push_back(prefixed_value);
    } else
    {
        /* Create a new compression value from the encoded prefix */
        const CompressionValue<N> prefixed_value(decoded_prefix, num_prefix_bits);
        reconstructed_variable_new.push_back(prefixed_value);
    }
}

template<int N>
void
DecompressionPrefixAdaptDataEGU<N>::KeepCompressionValue(const int index)
{
    if (count_adaptation_step_ > 0)
    {
        reconstructed_variable_new.push_back(reconstructed_variable[index]);
    } else
    {
        reconstructed_variable_new.push_back(CompressionValue<N>());
    }
}

template<int N>
void
DecompressionPrefixAdaptDataEGU<N>::KeepAndCopyCompressionValue(const int start_index, const int num_elements)
{
   // cmc_debug_msg("in apply prefix for range");
    if (count_adaptation_step_ > 0)
    {
        /* Set the compression value several times wihtin the new vector */
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, reconstructed_variable[start_index]);
    } else
    {
        /* Create a new compression value from the encoded prefix */
        std::fill_n(std::back_inserter(reconstructed_variable_new), num_elements, CompressionValue<N>());
    }
}

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
        if (IsMeshElementWithinGeoDomain(elements[0], ts, adapt_data->domain, initial_refinement_level, adapt_data->current_layout))
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
        if (IsMeshElementWithinGeoDomain(elements[0], ts, adapt_data->domain, initial_refinement_level, adapt_data->current_layout))
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

#endif /* CMC_WITH_T8CODE */

}

#endif /* !CMC_T8_ADAPT_CALLBACKS_H */
