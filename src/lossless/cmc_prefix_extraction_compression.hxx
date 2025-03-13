#ifndef CMC_PREFIX_EXTRACTION_COMPRESSION_HXX
#define CMC_PREFIX_EXTRACTION_COMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossless/cmc_byte_compression_variable.hxx"
//#include "utilities/cmc_iface_lossless_compression_data.hxx"


namespace cmc::lossless::prefix
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class PrefixAdaptData : public ICompressionAdaptData<T>
{
public:
    PrefixAdaptData() = delete;
    PrefixAdaptData(AbstractByteCompressionVariable<T>* variable)
    : ICompressionAdaptData<T>(variable) {};

    void InitializeExtractionIteration() override;
    void FinalizeExtractionIteration() override;
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::vector<uint8_t> EncodeLevelData(const std::vector<SerializedCompressionValue<sizeof(T)>>& level_byte_values) const override;
    
protected:
    ExtractionData<T> PerformExtraction(const int which_tree, const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values) override;
    UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    void IndicateCommonPrefixFound();
    void IndicateNoCommonPrefix();
    void IndicateCoarsening();
    void IndicateElementStaysUnchanged();

    bit_map::BitMap refinement_indications_;
    bit_map::BitMap prefix_indications_;
    int count_adaptation_step_{0};
};

template<typename T>
inline void
PrefixAdaptData<T>::IndicateCommonPrefixFound()
{
    prefix_indications_.AppendSetBit();
}

template<typename T>
inline void
PrefixAdaptData<T>::IndicateNoCommonPrefix()
{
    prefix_indications_.AppendUnsetBit();
}

template<typename T>
inline void
PrefixAdaptData<T>::IndicateCoarsening()
{
    refinement_indications_.AppendSetBit();
}

template<typename T>
inline void
PrefixAdaptData<T>::IndicateElementStaysUnchanged()
{
    refinement_indications_.AppendUnsetBit();
}

template <typename T>
void
PrefixAdaptData<T>::InitializeExtractionIteration()
{
    refinement_indications_ = bit_map::BitMap();
    prefix_indications_ = bit_map::BitMap();
}

template <typename T>
void
PrefixAdaptData<T>::FinalizeExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
PrefixAdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
PrefixAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template<typename T>
std::pair<bool, CompressionValue<T>>
EvaluateCommonPrefix(const VectorView<CompressionValue<T>>& compression_values)
{
    cmc_assert(compression_values.size() >= 2);

    /* Check if all elements are holding an actual prefix */
    for (auto cv_iter = compression_values.begin(); cv_iter != compression_values.end(); ++cv_iter)
    {
        if (cv_iter->IsEmpty())
        {
            /* Since this prefix value is empty, there cannot be a common prefix */
            return std::make_pair(false, CompressionValue<T>());
        }
    }

    /* Determine a common prefix of the first two values */
    CompressionValue<T> prefix = GetCommonPrefix<sizeof(T)>(compression_values[0], compression_values[1]);

    /* Check if there is common prefix between the first two values */
    if (prefix.GetCountOfSignificantBits() == 0)
    {
        /* There is no common prefix */
        return std::make_pair(false, CompressionValue<T>());
    }

    /* Check if there is a common prefix with the other values within the view */
    for (size_t index = 2; index < compression_values.size(); ++index)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix, compression_values[index]);

        if (prefix.GetCountOfSignificantBits() == 0)
        {
            /* There is no common prefix */
            return std::make_pair(false, CompressionValue<T>());
        }
    }

    /* If the function arrives here, we do have found a common prefix which can be extracted from the 'previous prefixes' */
    return std::make_pair(true, prefix);
}

template <typename T>
ExtractionData<T>
PrefixAdaptData<T>::PerformExtraction([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, [[maybe_unused]] const int num_elements, const VectorView<CompressionValue<T>> values)
{
    /* Since we perform an extraction, a family of elements is coarsened */
    IndicateCoarsening();

    /* Evaluate whether there is a common prefix to extract */
    auto [is_prefix_found, prefix] = EvaluateCommonPrefix<T>(values);

    if (is_prefix_found)
    {
        /* Indicate that we have found a common prefix */
        IndicateCommonPrefixFound();

        std::vector<CompressionValue<T>> fine_values;
        fine_values.reserve(values.size());

        /* We need to trim the previous prefixes by the extracted common prefix */
        for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
        {
            fine_values.push_back(*val_iter);
            fine_values.back().SetFrontBit(sizeof(T) * bit_map::kCharBit - prefix.GetTailBit());
        }

        /* Return the common prefix and the trimmed remaining errors */
        return ExtractionData<T>(prefix, std::move(fine_values));
    } else
    {
        /* Indicate that we have not found a common prefix */
        IndicateNoCommonPrefix();

        std::vector<CompressionValue<T>> fine_values;
        fine_values.reserve(values.size());

        /* We copy the fine values over and leave them unchanged */
        for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
        {
            fine_values.push_back(*val_iter);
        }

        /* We return an empty prefix and a the unchanged previous prefixes */
        return ExtractionData<T>(CompressionValue<T>(), std::move(fine_values));
    }
}

template <typename T>
UnchangedData<T>
PrefixAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* In case the element stays unchanged, we indicate that no coarsening is possible */
    IndicateElementStaysUnchanged();

    /* Therefore, it is not possible to extract a common prefix as well */
    IndicateNoCommonPrefix();

    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     *  as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return UnchangedData<T>(value, CompressionValue<T>());
}

template <typename T>
std::vector<uint8_t>
PrefixAdaptData<T>::EncodeLevelData(const std::vector<SerializedCompressionValue<sizeof(T)>>& level_byte_values) const
{
    return std::vector<uint8_t>();
}

template <typename T>
inline ICompressionAdaptData<T>*
CreatePrefixExtractionAdaptationClass(AbstractByteCompressionVariable<T>* abstract_var)
{
    return new PrefixAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyPrefixExtractionAdaptationClass(ICompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class CompressionVariable : public AbstractByteCompressionVariable<T>
{
public:
    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<T>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixExtractionAdaptationClass<T>;
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<SerializedCompressionValue<sizeof(T)>>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixExtractionAdaptationClass<T>;
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, std::vector<SerializedCompressionValue<sizeof(T)>>&& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(std::move(variable_data));
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixExtractionAdaptationClass<T>;
    };
};



}



#endif /* !CMC_PREFIX_EXTRACTION_COMPRESSION_HXX */
