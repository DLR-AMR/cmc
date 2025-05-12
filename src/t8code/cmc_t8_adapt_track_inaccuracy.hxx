#ifndef CMC_T8_ADAPT_TRACK_INACCURACY_HXX
#define CMC_T8_ADAPT_TRACK_INACCURACY_HXX
/**
 * @file cmc_t8_adapt_track_inaccuracy.hxx
 */

#include "t8code/cmc_t8_adapt_track_inaccuracy_forward.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_partition.h>
#endif

#include <cmath>
#include <memory>
#include <iostream>
#include <unordered_map>

#include "utilities/cmc_log_functions.hxx"
namespace cmc
{

template<class T>
class InaccuracyComputerSkipMissingValues
{
public:
    InaccuracyComputerSkipMissingValues() = delete;
    InaccuracyComputerSkipMissingValues(const ComputeInaccuracySkipMissingValues<T>& inaccuracy_computer, const ComputeSingleInaccuracySkipMissingValues<T>& single_inaccuracy_computer)
    : inaccuracy_computer_{inaccuracy_computer}, single_inaccuracy_computer_{single_inaccuracy_computer} {};
    
    ~InaccuracyComputerSkipMissingValues(){};
    
    InaccuracyComputerSkipMissingValues(const InaccuracyComputerSkipMissingValues& other) = default;
    InaccuracyComputerSkipMissingValues& operator=(const InaccuracyComputerSkipMissingValues& other) = default;
    InaccuracyComputerSkipMissingValues(InaccuracyComputerSkipMissingValues&& other) = default;
    InaccuracyComputerSkipMissingValues& operator=(InaccuracyComputerSkipMissingValues&& other) = default;

    std::vector<double> operator()(const VectorView<T>& values, const T& interpolated_value, const std::vector<double>& previous_deviation, const T& missing_value) const
    {
        return inaccuracy_computer_(values, interpolated_value, previous_deviation, missing_value);
    };

    double operator()(const T initial_value, const T interpolated_value, const T& missing_value) const
    {
        return single_inaccuracy_computer_(initial_value, interpolated_value, missing_value);
    };
private:
    ComputeInaccuracySkipMissingValues<T> inaccuracy_computer_;
    ComputeSingleInaccuracySkipMissingValues<T> single_inaccuracy_computer_;
};

template<typename T>
auto
ComputeRelativeDeviationSkipMissingValues(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        if (std::fpclassify(previous_absolute_max_deviation[index]) == FP_ZERO)
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                /* Calculate the relative deviation */
                deviations.push_back(static_cast<double>(std::abs(*iter - nominal_value)) / static_cast<double>(std::abs(*iter)));
            } else
            {
                deviations.push_back(0.0);
            }
        } else
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                const double zaehler = previous_absolute_max_deviation[index] + std::abs(static_cast<double>(*iter) - static_cast<double>(nominal_value));
                const double nenner = std::min(std::initializer_list<double>({static_cast<double>(*iter), static_cast<double>(*iter) - previous_absolute_max_deviation[index], static_cast<double>(*iter) + previous_absolute_max_deviation[index]}));
                deviations.push_back(zaehler / nenner);
            } else
            {
                deviations.push_back(0.0);
            }
        }
    }

    return deviations;
}


template<typename T>
auto
ComputeRelativeDeviationSkipMissingValues(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        if (std::fpclassify(previous_absolute_max_deviation[index]) == FP_ZERO)
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                /* Calculate the relative deviation */
                deviations.emplace_back(static_cast<double>((*iter > nominal_value ? *iter - nominal_value : nominal_value - *iter)) / static_cast<double>(*iter));
            } else
            {
                deviations.emplace_back(0.0);
            }
        } else
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                const double zaehler = (*iter > nominal_value ? static_cast<double>(*iter) - static_cast<double>(nominal_value) : static_cast<double>(nominal_value) - static_cast<double>(*iter)) + previous_absolute_max_deviation[index];
                const double nenner = std::min(std::initializer_list<double>({static_cast<double>(*iter), static_cast<double>(*iter) - previous_absolute_max_deviation[index], static_cast<double>(*iter) + previous_absolute_max_deviation[index]}));
                deviations.push_back(zaehler / nenner);
            } else
            {
                deviations.emplace_back(0.0);
            }
        }
    }

    return deviations;
}


template<typename T>
auto
ComputeRelativeDeviationSkipMissingValues(const T value, const T nominal_value, const double previous_absolute_max_deviation, const T missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (std::fpclassify(previous_absolute_max_deviation) == FP_ZERO)
    {
        if (!ApproxCompare(missing_value, value))
        {
            /* Calculate the relative deviation */
            return (static_cast<double>(std::abs(value - nominal_value)) / static_cast<double>(std::abs(value)));
        } else
        {
            return 0.0;
        }
    } else
    {
        if (!ApproxCompare(missing_value, value))
        {
            /* Estimate the relative deviation */
            const double zaehler = previous_absolute_max_deviation + std::abs(static_cast<double>(value) - static_cast<double>(nominal_value));
            const double nenner = std::min(std::initializer_list<double>({static_cast<double>(value), static_cast<double>(value) - previous_absolute_max_deviation, static_cast<double>(value) + previous_absolute_max_deviation}));
            return (zaehler / nenner);
        } else
        {
            return 0.0;
        }
    }

    return 0.0;
}


template<typename T>
auto
ComputeRelativeDeviationSkipMissingValues(const T value, const T nominal_value, const double previous_absolute_max_deviation, const T missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
        if (std::fpclassify(previous_absolute_max_deviation) == FP_ZERO)
        {
            if (!ApproxCompare(missing_value, value))
            {
                /* Calculate the relative deviation */
                return (static_cast<double>((value > nominal_value ? value - nominal_value : nominal_value - value)) / static_cast<double>(value));
            } else
            {
                return 0.0;
            }
        } else
        {
            if (!ApproxCompare(missing_value, value))
            {
                const double zaehler = (value > nominal_value ? static_cast<double>(value) - static_cast<double>(nominal_value) : static_cast<double>(nominal_value) - static_cast<double>(value)) + previous_absolute_max_deviation;
                const double nenner = std::min(std::initializer_list<double>({static_cast<double>(value), static_cast<double>(value) - previous_absolute_max_deviation, static_cast<double>(value) + previous_absolute_max_deviation}));
                return (zaehler / nenner);
            } else
            {
                return 0.0;
            }
        }

    return 0.0;
}

template<typename T>
auto
ComputeSingleRelativeDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the relative deviation */
        return static_cast<double>(std::abs(initial_value - nominal_value)) / static_cast<double>(std::abs(initial_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleRelativeDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the relative deviation */
        return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value)) / static_cast<double>(initial_value);
    } else
    {
        return 0.0;
    }
}


template<typename T>
auto
ComputeAbsoluteDeviationSkipMissingValues(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        if (!ApproxCompare(missing_value, *iter))
        {
            /* Calculate the absolute deviation */
            deviations.emplace_back(static_cast<double>(std::abs(*iter - nominal_value)) + previous_absolute_max_deviation[index]);
        } else
        {
            deviations.emplace_back(0.0);
        }
    }

    return deviations;
}

template<typename T>
auto
ComputeAbsoluteDeviationSkipMissingValues(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        if (!ApproxCompare(missing_value, *iter))
        {
            /* Calculate the absolute deviation */
            deviations.emplace_back(static_cast<double>((*iter > nominal_value ? *iter - nominal_value : nominal_value - *iter)) + previous_absolute_max_deviation[index]);
        } else
        {
            deviations.emplace_back(0.0);
        }
    }

    return deviations;
}

template<typename T>
auto
ComputeAbsoluteDeviationSkipMissingValues(const T value, const T nominal_value, const double previous_absolute_max_deviation, const T missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (!ApproxCompare(missing_value, value))
    {
        /* Calculate the absolute deviation */
        return (static_cast<double>(std::abs(value - nominal_value)) + previous_absolute_max_deviation);
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeAbsoluteDeviationSkipMissingValues(const T value, const T nominal_value, const double previous_absolute_max_deviation, const T missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    if (!ApproxCompare(missing_value, value))
    {
        /* Calculate the absolute deviation */
        return (static_cast<double>((value > nominal_value ? value - nominal_value : nominal_value - value)) + previous_absolute_max_deviation);
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleAbsoluteDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the absolute deviation */
        return static_cast<double>(std::abs(initial_value - nominal_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleAbsoluteDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the absolute deviation */
        return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value));
    } else
    {
        return 0.0;
    }
}

template<class T>
class InaccuracyComputer
{
public:
    InaccuracyComputer() = delete;
    InaccuracyComputer(const ComputeInaccuracy<T>& inaccuracy_computer, const ComputeSingleInaccuracy<T>& single_inaccuracy_computer)
    : inaccuracy_computer_{inaccuracy_computer}, single_inaccuracy_computer_{single_inaccuracy_computer} {};
    
    ~InaccuracyComputer(){};
    
    InaccuracyComputer(const InaccuracyComputer& other) = default;
    InaccuracyComputer& operator=(const InaccuracyComputer& other) = default;
    InaccuracyComputer(InaccuracyComputer&& other) = default;
    InaccuracyComputer& operator=(InaccuracyComputer&& other) = default;

    std::vector<double> operator()(const VectorView<T>& values, const T& interpolated_value, const std::vector<double>& previous_deviation) const
    {
        return inaccuracy_computer_(values, interpolated_value, previous_deviation);
    };

    double operator()(const T initial_value, const T interpolated_value) const
    {
        return single_inaccuracy_computer_(initial_value, interpolated_value);
    };
private:
    ComputeInaccuracy<T> inaccuracy_computer_;
    ComputeSingleInaccuracy<T> single_inaccuracy_computer_;
};

template<typename T>
auto
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        if (std::fpclassify(previous_absolute_max_deviation[index]) == FP_ZERO)
        {

            /* Calculate the relative deviation */
            deviations.push_back(static_cast<double>(std::abs(*iter - nominal_value)) / static_cast<double>(std::abs(*iter)));
        } else
        {
            const double zaehler = previous_absolute_max_deviation[index] + std::abs(static_cast<double>(*iter) - static_cast<double>(nominal_value));
            const double nenner = std::min(std::initializer_list<double>({static_cast<double>(*iter), static_cast<double>(*iter) - previous_absolute_max_deviation[index], static_cast<double>(*iter) + previous_absolute_max_deviation[index]}));
            deviations.push_back(zaehler / nenner);
        }
    }

    return deviations;
}


template<typename T>
auto
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        if (std::fpclassify(previous_absolute_max_deviation[index]) == FP_ZERO)
        {
            /* Calculate the relative deviation */
            deviations.emplace_back(static_cast<double>((*iter > nominal_value ? *iter - nominal_value : nominal_value - *iter)) / static_cast<double>(*iter));
        } else
        {
            const double zaehler = (*iter > nominal_value ? static_cast<double>(*iter) - static_cast<double>(nominal_value) : static_cast<double>(nominal_value) - static_cast<double>(*iter)) + previous_absolute_max_deviation[index];
            const double nenner = std::min(std::initializer_list<double>({static_cast<double>(*iter), static_cast<double>(*iter) - previous_absolute_max_deviation[index], static_cast<double>(*iter) + previous_absolute_max_deviation[index]}));
            deviations.push_back(zaehler / nenner);
        }
    }

    return deviations;
}


template<typename T>
auto
ComputeRelativeDeviation(const T value, const T nominal_value, const double previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (std::fpclassify(previous_absolute_max_deviation) == FP_ZERO)
    {
        /* Calculate the relative deviation */
        return (static_cast<double>(std::abs(value - nominal_value)) / static_cast<double>(std::abs(value)));
    } else
    {
        /* Estimate the relative deviation */
        const double zaehler = previous_absolute_max_deviation + std::abs(static_cast<double>(value) - static_cast<double>(nominal_value));
        const double nenner = std::min(std::initializer_list<double>({static_cast<double>(value), static_cast<double>(value) - previous_absolute_max_deviation, static_cast<double>(value) + previous_absolute_max_deviation}));
        return (zaehler / nenner);
    }

    return 0.0;
}


template<typename T>
auto
ComputeRelativeDeviation(const T value, const T nominal_value, const double previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
        if (std::fpclassify(previous_absolute_max_deviation) == FP_ZERO)
        {
            /* Calculate the relative deviation */
            return (static_cast<double>((value > nominal_value ? value - nominal_value : nominal_value - value)) / static_cast<double>(value));
        } else
        {
            const double zaehler = (value > nominal_value ? static_cast<double>(value) - static_cast<double>(nominal_value) : static_cast<double>(nominal_value) - static_cast<double>(value)) + previous_absolute_max_deviation;
            const double nenner = std::min(std::initializer_list<double>({static_cast<double>(value), static_cast<double>(value) - previous_absolute_max_deviation, static_cast<double>(value) + previous_absolute_max_deviation}));
            return (zaehler / nenner);
        }

    return 0.0;
}

template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the relative deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value)) / static_cast<double>(std::abs(initial_value));
}

template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the relative deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value)) / static_cast<double>(initial_value);
}


template<typename T>
auto
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        /* Calculate the absolute deviation */
        deviations.emplace_back(static_cast<double>(std::abs(*iter - nominal_value)) + previous_absolute_max_deviation[index]);
    }

    return deviations;
}

template<typename T>
auto
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>
{
    cmc_assert(values.size() == previous_absolute_max_deviation.size());

    std::vector<double> deviations;
    deviations.reserve(values.size());

    int index = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++index)
    {
        /* Calculate the absolute deviation */
        deviations.emplace_back(static_cast<double>((*iter > nominal_value ? *iter - nominal_value : nominal_value - *iter)) + previous_absolute_max_deviation[index]);
    }

    return deviations;
}

template<typename T>
auto
ComputeAbsoluteDeviation(const T value, const T nominal_value, const double previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the absolute deviation */
    return (static_cast<double>(std::abs(value - nominal_value)) + previous_absolute_max_deviation);
}

template<typename T>
auto
ComputeAbsoluteDeviation(const T value, const T nominal_value, const double previous_absolute_max_deviation)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the absolute deviation */
    return (static_cast<double>((value > nominal_value ? value - nominal_value : nominal_value - value)) + previous_absolute_max_deviation);
}

template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value));
}

template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value));
}

class InaccuracyContainer
{
public:
    virtual double GetInaccuracy(const int index) = 0;
    virtual std::vector<double> GetInaccuracyForRange(const int start_index, const int num_elements) = 0;
    virtual void StoreInaccuracy(const int index, const double inaccuracy) = 0;
    virtual void TransferPreviousDeviations(const int start_index_previous_values, const int num_previous_values) = 0;
    virtual void TransferPreviousDeviation(const int start_index_previous_values) = 0;
    virtual void SwitchDeviations() = 0;
    virtual void AllocateDeviationStorage(const int num_elements) = 0;
    virtual void RepartitionDeviations(t8_forest_t initial_forest, t8_forest_t partitioned_forest) = 0;

    virtual InaccuracyContainer* clone() const = 0;

    virtual ~InaccuracyContainer() = default;
};

/* Track the inaccurcy of all elements throughout the whole adaptation.
 * This allows (for example) for additional lossy comrpession steps after the AMR lossy compression */
class FullInaccuracyTracker : public InaccuracyContainer
{
public:
    FullInaccuracyTracker() = default;
    FullInaccuracyTracker(size_t size_hint)
    {
        if (size_hint != kInvalidSizeHintForInaccuracyContainer)
        {
            deviations_.reserve(size_hint);
        }
    };

    FullInaccuracyTracker(const FullInaccuracyTracker& other) = default;
    FullInaccuracyTracker& operator=(const FullInaccuracyTracker& other) = default;
    FullInaccuracyTracker(FullInaccuracyTracker&& other) = default;
    FullInaccuracyTracker& operator=(FullInaccuracyTracker&& other) = default;
    
    virtual FullInaccuracyTracker* clone() const {
        return new FullInaccuracyTracker(*this);
    };

    ~FullInaccuracyTracker() = default;

    double GetInaccuracy(const int id) override
    {
        if (not previous_deviations_.empty())
        {
            return previous_deviations_[id];
        } else
        {
            return static_cast<double>(0.0);
        }
    };

    std::vector<double> GetInaccuracyForRange(const int start_index, const int num_elements) override
    {
        std::vector<double> deviations;
        deviations.reserve(num_elements);

        if (not previous_deviations_.empty())
        {
            for (int iter = start_index; iter < start_index + num_elements; ++iter)
            {
                deviations.push_back(previous_deviations_[iter]);
            }
        } else
        {
            std::fill_n(std::back_inserter(deviations), num_elements, 0.0);
        }

        return deviations;
    }

    void StoreInaccuracy([[maybe_unused]] const int index, const double inaccuracy) override
    {
        deviations_.push_back(inaccuracy);
    };

    void TransferPreviousDeviations(const int start_index_previous_values, const int num_previous_values) override
    {
        std::copy(previous_deviations_.begin() + start_index_previous_values,
                  previous_deviations_.begin() + start_index_previous_values + num_previous_values,
                  std::back_inserter(deviations_));
    }
    void TransferPreviousDeviation(const int start_index_previous_values) override
    {
        deviations_.push_back(previous_deviations_[start_index_previous_values]);
    }

    void SwitchDeviations() override
    {
        std::swap(previous_deviations_, deviations_);
        deviations_.clear();
    }

    void AllocateDeviationStorage(const int num_elements) override
    {
        if (previous_deviations_.empty())
        {
            previous_deviations_ = std::vector<double>(num_elements, 0.0);
        }

        deviations_.reserve(num_elements);
    }

    void RepartitionDeviations(t8_forest_t initial_forest, t8_forest_t partitioned_forest) override
    {
        /* Create an sc_array_t wrapper of the current tracked deviations */
        sc_array_t* in_data = sc_array_new_data (static_cast<void*>(previous_deviations_.data()), sizeof(double), previous_deviations_.size());

        /* Allocate memory for the partitioned data */
        const t8_locidx_t new_num_elems = t8_forest_get_local_num_elements(partitioned_forest);
        deviations_ = std::vector<double>(new_num_elems);

        /* Create a wrapper for the freshly allocated partitioned deviations */
        sc_array_t* out_data = sc_array_new_data (static_cast<void*>(deviations_.data()), sizeof(double), deviations_.size());

        /* Partition the deviations compliant to the new partitioned forest */
        t8_forest_partition_data(initial_forest, partitioned_forest, in_data, out_data);

        /* Destroy the array wrappers */
        sc_array_destroy(in_data);
        sc_array_destroy(out_data);

        /* Set the variable's data to the newly partitioned data */
        SwitchDeviations();

        cmc_debug_msg("Repartition Deviations is finsihed");
    }

private:
    std::vector<double> previous_deviations_;
    std::vector<double> deviations_;
};

/* Track the inaccurcay of the last adaptation step only. This saves storage, but after the compression,
 * we are not able to check the deviation of each element. 
 */
//TODO: implement
class MinimalInaccuracyTracker : public InaccuracyContainer
{
public:
    MinimalInaccuracyTracker() = default;
    MinimalInaccuracyTracker(const size_t size_hint)
    {
        if (size_hint != kInvalidSizeHintForInaccuracyContainer)
        {
            deviations_.reserve(size_hint);
        }
    };
    
    MinimalInaccuracyTracker(const MinimalInaccuracyTracker& other) = default;
    MinimalInaccuracyTracker& operator=(const MinimalInaccuracyTracker& other) = default;
    MinimalInaccuracyTracker(MinimalInaccuracyTracker&& other) = default;
    MinimalInaccuracyTracker& operator=(MinimalInaccuracyTracker&& other) = default;
    
    virtual MinimalInaccuracyTracker* clone() const {
        return new MinimalInaccuracyTracker(*this);
    };

    ~MinimalInaccuracyTracker() = default;
    
    double GetInaccuracy([[maybe_unused]] const int id) override {return 0.0;};
    std::vector<double> GetInaccuracyForRange([[maybe_unused]] const int start_index, [[maybe_unused]] const int num_elements) override
    {
        return std::vector<double>();
    }
    void StoreInaccuracy([[maybe_unused]] const int index, [[maybe_unused]] const double inaccuracy) override {};
    void SwitchDeviations() override {};
    void TransferPreviousDeviations([[maybe_unused]] const int start_index_previous_values, [[maybe_unused]] const int num_previous_values) override
    {};
    void TransferPreviousDeviation([[maybe_unused]] const int start_index_previous_values) override {};
    void AllocateDeviationStorage([[maybe_unused]] const int num_elements) override {};
    void RepartitionDeviations([[maybe_unused]] t8_forest_t initial_forest, [[maybe_unused]] t8_forest_t partitioned_forest) override {};

private:
    std::unordered_map<int, double> deviations_;
};


}

#endif /* !CMC_T8_ADAPT_TRACK_INACCURACY_HXX */
