#ifndef CMC_T8_ADAPT_TRACK_INACCURACY_HXX
#define CMC_T8_ADAPT_TRACK_INACCURACY_HXX
/**
 * @file cmc_t8_adapt_track_inaccuracy.hxx
 */

#include "t8code/cmc_t8_adapt_track_inaccuracy_forward.hxx"

#include <cmath>
#include <memory>
#include <iostream>
#include <unordered_map>

namespace cmc
{

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

    std::vector<double> operator()(const VectorView<T>& values, const T& interpolated_value, const std::vector<double>& previous_deviation, const T& missing_value) const
    {
        return inaccuracy_computer_(values, interpolated_value, previous_deviation, missing_value);
    };

    double operator()(const T initial_value, const T interpolated_value, const T& missing_value) const
    {
        return single_inaccuracy_computer_(initial_value, interpolated_value, missing_value);
    };
private:
    ComputeInaccuracy<T> inaccuracy_computer_;
    ComputeSingleInaccuracy<T> single_inaccuracy_computer_;
};

template<typename T>
auto
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
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
                deviations.push_back(std::max(std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation[index] - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation[index]),
                                              std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation[index] - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation[index])));
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
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
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
                deviations.emplace_back(std::max(std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation[index] - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation[index]),
                                                 std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation[index] - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation[index])));
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
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value, const T& missing_value)
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
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value, const T& missing_value)
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
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
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
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value)
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
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value, const T& missing_value)
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
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value, const T& missing_value)
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
        if (previous_deviations_.size() > 0)
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

        for (int iter = start_index; iter < start_index + num_elements; ++iter)
        {
            deviations.push_back(previous_deviations_[iter]);
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
private:
    std::unordered_map<int, double> deviations_;
};






















#if 0


template<class T>
class InaccuracyComputer
{
public:
    InaccuracyComputer() = delete;
    InaccuracyComputer(const ComputeInaccuracy<T>& inaccuracy_computer)
    : inaccuracy_computer_{inaccuracy_computer}{};
    
    ~InaccuracyComputer(){};
    
    InaccuracyComputer(const InaccuracyComputer& other) = default;
    InaccuracyComputer& operator=(const InaccuracyComputer& other) = default;
    InaccuracyComputer(InaccuracyComputer&& other) = default;
    InaccuracyComputer& operator=(InaccuracyComputer&& other) = default;

    std::vector<double> operator()(const VectorView<T>& values, const T& interpolated_value, const double previous_deviation, const T& missing_value) const
    {
        return inaccuracy_computer_(values, interpolated_value, previous_deviation, missing_value);
    };

private:
    ComputeInaccuracy<T> inaccuracy_computer_;
};

template<typename T>
auto
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>
{
    std::vector<double> deviations;
    deviations.reserve(values.size());

    if (previous_absolute_max_deviation == 0.0)
    {
        for (auto iter = values.begin(); iter != values.end(); ++iter)
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                /* Calculate the relative deviation */
                deviations.emplace_back(static_cast<double>(std::abs(*iter - nominal_value)) / static_cast<double>(std::abs(*iter)));
            } else
            {
                deviations.emplace_back(0.0);
            }
        }
    } else
    {
        for (auto iter = values.begin(); iter != values.end(); ++iter)
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                deviations.emplace_back(std::max(std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation),
                                                 std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation)));
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
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>
{
    std::vector<double> deviations;
    deviations.reserve(values.size());

    if (previous_absolute_max_deviation == 0.0)
    {
        for (auto iter = values.begin(); iter != values.end(); ++iter)
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                /* Calculate the relative deviation */
                deviations.emplace_back(static_cast<double>((*iter > nominal_value ? *iter - nominal_value : nominal_value - *iter)) / static_cast<double>(*iter));
            } else
            {
                deviations.emplace_back(0.0);
            }
        }
    } else
    {
        for (auto iter = values.begin(); iter != values.end(); ++iter)
        {
            if (!ApproxCompare(missing_value, *iter))
            {
                deviations.emplace_back(std::max(std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation),
                                                 std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation)));
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
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>
{
    std::vector<double> deviations;
    deviations.reserve(values.size());

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        if (!ApproxCompare(missing_value, *iter))
        {
            /* Calculate the relative deviation */
            deviations.emplace_back(static_cast<double>(std::abs(*iter - nominal_value)) + previous_absolute_max_deviation);
        } else
        {
            deviations.emplace_back(0.0);
        }
    }

    return deviations;
}

template<typename T>
auto
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>
{
    std::vector<double> deviations;
    deviations.reserve(values.size());

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        if (!ApproxCompare(missing_value, *iter))
        {
            /* Calculate the relative deviation */
            deviations.emplace_back(static_cast<double>((*iter > nominal_value ? *iter - nominal_value : nominal_value - *iter)) + previous_absolute_max_deviation);
        } else
        {
            deviations.emplace_back(0.0);
        }
    }

    return deviations;
}

class InaccuracyContainer
{
public:
    virtual double GetInaccuracy(const int id) = 0;
    virtual std::vector<double> GetInaccuracyForRange(const int start_index, const int num_elements)
    virtual void StoreInaccuracy(const double inaccuracy, const int id) = 0;

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

    //void SetUpTrackerForElementsInForest(const size_t num_elements);

    double GetInaccuracy(const int id) override
    {
        if (previous_deviations_.size() > 0)
        {
            return previous_deviations_[id];
        } else
        {
            return static_cast<double>(0.0);
        }
    };

    std::vector<double> GetInaccuracyForRange(const int start_index, const int num_elements)
    {
        std::vector<double> deviations;
        deviations.reserve(num_elements);

        for (int iter = start_index; iter < start_index + num_elements; ++iter)
        {
            deviations.push_back(previous_deviations_[iter]);
        }

        return deviations;
    }

    void StoreInaccuracy(const double inaccuracy, const int index) override
    {
        cmc_assert(index >= deviations_.size());
        deviations_[index] = inaccuracy;
    };

private:
    std::vector<double> previous_deviations_;
    std::vector<double> deviations_;
};

/* Track the inaccurcay of the last adaptation step only. This saves storage, but after the compression,
 * we are not able to check the deviation of each element. 
 */
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
    
    double GetInaccuracy(const int id) override {return 0.0;};
    std::vector<double> GetInaccuracyForRange(const int start_index, const int num_elements)
    {
        return std::vector<double>();
    }
    void StoreInaccuracy(const double inaccuracy, const int id) override {};

private:
    std::unordered_map<int, double> deviations_;
};


#endif





















#if 0

template<typename T>
using ComputeInaccuracy = std::vector<double>(*)(const Variable<T>& variable, const std::span<T> values, const T& interpolated_value, const double previous_absolute_max_deviation = 0.0);

template<class T>
InaccuracyComputer
{
public:
    InaccuracyComputer() = delete;
    InaccuracyComputer(const ComputeInaccuracy<T>& inaccuracy_computer)
    : inaccuracy_computer_{inaccuracy_computer}{};
    
    ~InaccuracyComputer(){};
    
    InaccuracyComputer(const InaccuracyComputer& other) = default;
    InaccuracyComputer& operator=(const InaccuracyComputer& other) = default;
    InaccuracyComputer(InaccuracyComputer&& other) = default;
    InaccuracyComputer& operator=(InaccuracyComputer&& other) = default;

    std::vector<double> operator()(const Variable<T>& variable, const std::span<T> values, const T& interpolated_value, const double previous_deviation = 0.0) const
    {

        return inaccuracy_computer_(variable, values, interpolated_value, previous_deviation);
    };

private:
    ComputeInaccuracy inaccuracy_computer_;
};

//TODO: Add comparison to zero for denominator
template<typename T>
std::vector<double>
ComputeRelativeDeviation(const Variable<T>& variable, const std::span<T> values, const T& nominal_value, const double previous_absolute_max_deviation = 0.0)
{
    std::vector<double> deviations;
    deviations.reserve(values.size);

    const T missing_value = variable.GetMissingValue();

    if (previous_absolute_max_deviation == 0)
    {
        for (auto iter = values.begin(); iter != values.end(); ++iter)
        {
            if (!EqualToMissingValue(missing_value, *iter))
            {
                /* Calculate the relative deviation */
                deviations.emplace_back(static_cast<double>(std::abs(*iter - nominal_value)) / static_cast<double>(std::abs(*iter)));
            } else
            {
                deviations.emplace_back(0.0);
            }
        }
    } else
    {
        for (auto iter = values.begin(); iter != values.end(); ++iter)
        {
            if (!EqualToMissingValue(missing_value, *iter))
            {
                deviations.emplace_back(std::max(std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) + previous_absolute_max_deviation),
                                                 std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation - static_cast<double>(nominal_value)) / std::abs(static_cast<double>(*iter) - previous_absolute_max_deviation)));
            } else
            {
                deviations.emplace_back(0.0);
            }
        }
    }

    return deviations;
}


template<typename T>
std::vector<double>
ComputeAbsoluteDeviation(const Variable<T>& variable, const std::span<T> values, const T& nominal_value, const double previous_absolute_max_deviation = 0.0)
{
    std::vector<double> deviations;
    deviations.reserve(values.size);

    const T missing_value = variable.GetMissingValue();

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        if (!EqualToMissingValue(missing_value, *iter))
        {
            /* Calculate the relative deviation */
            deviations.emplace_back(static_cast<double>(std::abs(*iter - nominal_value)) + previous_absolute_max_deviation);
        } else
        {
            deviations.emplace_back(0.0);
        }
    }

    return deviations;
}

class InaccuracyContainer
{
public:
    virtual double GetInaccuracy(const size_t id) = 0;
    virtual void StoreInaccuracy(const double inaccuracy, const size_t id) = 0;

    virtual ~InaccuracyContainer(){
        std::cout << "InaccuracyContainer Destructor called\n";
    };
};

/* Track the inaccurcy of all elements throughout the whole adaptation.
 * This allows for additional lossy comrpession steps after the AMR lossy compression */
class FullInaccuracyTracker : public InaccuracyContainer
{
public:
    FullInaccuracyTracker() = default;
    FullInaccuracyTracker(const size_t size_hint)
    {
        deviations_.reserve(size_hint);
    };

    FullInaccuracyTracker(const FullInaccuracyTracker& other) = default;
    FullInaccuracyTracker& operator=(const FullInaccuracyTracker& other) = default;
    FullInaccuracyTracker(FullInaccuracyTracker&& other) = default;
    FullInaccuracyTracker& operator=(FullInaccuracyTracker&& other) = default;
    
    ~FullInaccuracyTracker(){
        std::cout << "FullInaccuracyTracker Destructor called\n";
    };

    double GetInaccuracy(const size_t id) override;
    void StoreInaccuracy(const double inaccuracy, const size_t id) override;

private:
    std::vector<double> deviations_;
};

/* Track the inaccurcay of the last adaptation step only. This saves storage, but after the compression,
 * we are not able to check the deviation of each element. 
 */
class MinimalInaccuracyTracker : public InaccuracyContainer
{
public:
    MinimalInaccuracyTracker() = default;
    MinimalInaccuracyTracker(const size_t size_hint)
    {
        deviations_.reserve(size_hint / 4);
    };
    
    MinimalInaccuracyTracker(const MinimalInaccuracyTracker& other) = default;
    MinimalInaccuracyTracker& operator=(const MinimalInaccuracyTracker& other) = default;
    MinimalInaccuracyTracker(MinimalInaccuracyTracker&& other) = default;
    MinimalInaccuracyTracker& operator=(MinimalInaccuracyTracker&& other) = default;
    
    ~MinimalInaccuracyTracker(){
        std::cout << "MinimalInaccuracyTracker Destructor called\n";
    };
    
    double GetInaccuracy(const size_t id) override;
    void StoreInaccuracy(const double inaccuracy, const size_t id) override;

private:
    std::unordered_map<size_t, double> deviations_;
};

enum TrackingOption {TrackFullInaccuracy = 0, TrackMinimalWorkingInaccuracy};

class InaccuracyTracker
{
public:
    InaccuracyTracker() = delete;
    InaccuracyTracker(const TrackingOption tracking_option, int i);
    
    InaccuracyTracker(const InaccuracyTracker& other) = default;
    InaccuracyTracker& operator=(const InaccuracyTracker& other) = default;
    InaccuracyTracker(InaccuracyTracker&& other) = default;
    InaccuracyTracker& operator=(InaccuracyTracker&& other) = default;
    
    ~InaccuracyTracker() = default;
    
    void CheckInaccuracy();

private:
    const TrackingOption tracking_option_;
    std::unique_ptr<InaccuracyContainer> deviations_;
};



#endif

















#if 0


//JUst use the tacker container for storing the deviations and use something else to calculate the deviations
enum TrackingOption {TrackFullInaccuracy = 0, TrackMinimalWorkingInaccuracy};

class InaccuracyContainer
{
public:
    virtual double GetInaccuracy(const size_t id) = 0;
    virtual void StoreInaccuracy(const double inaccuracy, const size_t id) = 0;
    virtual bool CheckInaccuracy(const double inaccuracy_to_compare, const size_t id) = 0;

    virtual ~InaccuracyContainer(){
        std::cout << "InaccuracyContainer Destructor called\n";
    };
};

enum TrackingOption {TrackFullInaccuracy = 0, TrackMinimalWorkingInaccuracy};

class InaccuracyContainer
{
public:
    virtual double GetInaccuracy(const size_t id) = 0;
    virtual void StoreInaccuracy(const double inaccuracy, const size_t id) = 0;
    virtual bool CheckInaccuracy(const double inaccuracy_to_compare, const size_t id) = 0;
    
    virtual ~InaccuracyContainer(){
        std::cout << "InaccuracyContainer Destructor called\n";
    };
};

class InaccuracyTrackerType
{
public:
    virtual std::unique_ptr<InaccuracyTrackerType> clone() const = 0;
    virtual double ComputeInaccuracy() = 0;
    virtual ~InaccuracyTrackerType(){
        std::cout << "InaccuracyTrackerType Destructor called\n";
    };
};

class AbsoluteInaccuracyTracker : public InaccuracyTrackerType
{
public:

    AbsoluteInaccuracyTracker() = default;
    
    AbsoluteInaccuracyTracker(const AbsoluteInaccuracyTracker& other) = default;
    AbsoluteInaccuracyTracker& operator=(const AbsoluteInaccuracyTracker& other) = default;
    AbsoluteInaccuracyTracker(AbsoluteInaccuracyTracker&& other) = default;
    AbsoluteInaccuracyTracker& operator=(AbsoluteInaccuracyTracker&& other) = default;
    
    ~AbsoluteInaccuracyTracker(){
        std::cout << "AbsoluteInaccuracyTracker Destructor called\n";
    };
    
    double ComputeInaccuracy() override;
    
    virtual std::unique_ptr<InaccuracyTrackerType> clone() const override
    {
        std::cout << "clone is called" << std::endl;
        return std::make_unique<AbsoluteInaccuracyTracker>( *this ); 
    };
};

class RelativeInaccuracyTracker : public InaccuracyTrackerType
{
public:

    RelativeInaccuracyTracker() = default;
    
    RelativeInaccuracyTracker(const RelativeInaccuracyTracker& other) = default;
    RelativeInaccuracyTracker& operator=(const RelativeInaccuracyTracker& other) = default;
    RelativeInaccuracyTracker(RelativeInaccuracyTracker&& other) = default;
    RelativeInaccuracyTracker& operator=(RelativeInaccuracyTracker&& other) = default;
    
    ~RelativeInaccuracyTracker(){
        std::cout << "RelativeInaccuracyTracker Destructor called\n";
    };
    
    double ComputeInaccuracy() override;
    
    virtual std::unique_ptr<InaccuracyTrackerType> clone() const override
    {
        std::cout << "clone is called" << std::endl;
        return std::make_unique<RelativeInaccuracyTracker>( *this ); 
    };
    
};

/* Track the inaccurcy of all elements throughout the whole adaptation.
 * This allows for additional lossy comrpession steps after the AMR lossy compression */
class FullInaccuracyTracker : public InaccuracyContainer
{
public:
    FullInaccuracyTracker() = delete;
    FullInaccuracyTracker(InaccuracyTrackerType&& inaccuracy_computer)
    : inaccuracy_computer_{inaccuracy_computer.clone()} {};
    FullInaccuracyTracker(std::unique_ptr<InaccuracyTrackerType>&& inaccuracy_computer)
    : inaccuracy_computer_{std::move(inaccuracy_computer)} {};

    bool CheckInaccuracy() override;
    
    FullInaccuracyTracker(const FullInaccuracyTracker& other) = default;
    FullInaccuracyTracker& operator=(const FullInaccuracyTracker& other) = default;
    FullInaccuracyTracker(FullInaccuracyTracker&& other) = default;
    FullInaccuracyTracker& operator=(FullInaccuracyTracker&& other) = default;
    
    ~FullInaccuracyTracker(){
        std::cout << "FullInaccuracyTracker Destructor called\n";
    };

private:
    std::unique_ptr<InaccuracyTrackerType> inaccuracy_computer_;
    std::vector<double> deviations_;
};

/* Track the inaccurcay of the last adaptation step only. This saves storage, but after the compression,
 * we are not able to check the deviation of each element. 
 */
class MinimalInaccuracyTracker : public InaccuracyContainer
{
public:
    MinimalInaccuracyTracker(InaccuracyTrackerType&& inaccuracy_computer)
    : inaccuracy_computer_{inaccuracy_computer.clone()} {}; 
    MinimalInaccuracyTracker(std::unique_ptr<InaccuracyTrackerType>&& inaccuracy_computer)
    : inaccuracy_computer_{std::move(inaccuracy_computer)} {};
    
    bool CheckInaccuracy() override;
    
    MinimalInaccuracyTracker(const MinimalInaccuracyTracker& other) = default;
    MinimalInaccuracyTracker& operator=(const MinimalInaccuracyTracker& other) = default;
    MinimalInaccuracyTracker(MinimalInaccuracyTracker&& other) = default;
    MinimalInaccuracyTracker& operator=(MinimalInaccuracyTracker&& other) = default;
    
    ~MinimalInaccuracyTracker(){
        std::cout << "MinimalInaccuracyTracker Destructor called\n";
    };
    
private:
    std::unique_ptr<InaccuracyTrackerType> inaccuracy_computer_;
    std::unordered_map<int, double> deviations;
};

class InaccuracyTracker
{
public:
    InaccuracyTracker() = delete;
    InaccuracyTracker(const TrackingOption tracking_option, int i);
    
    InaccuracyTracker(const InaccuracyTracker& other) = default;
    InaccuracyTracker& operator=(const InaccuracyTracker& other) = default;
    InaccuracyTracker(InaccuracyTracker&& other) = default;
    InaccuracyTracker& operator=(InaccuracyTracker&& other) = default;
    
    ~InaccuracyTracker() = default;
    
    void CheckInaccuracy();

private:
    const TrackingOption tracking_option_;
    std::unique_ptr<InaccuracyContainer> deviations_;
};

#endif


}

#endif /* !CMC_T8_ADAPT_TRACK_INACCURACY_HXX */

#if 0

int main()
{
    
    InaccuracyTracker tracker(TrackFullInaccuracy, 0);
    
    tracker.CheckInaccuracy();
    
    return 0;
}


#endif