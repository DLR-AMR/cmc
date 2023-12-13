#ifndef CMC_T8_ADAPT_TRACK_INACCURACY_HXX
#define CMC_T8_ADAPT_TRACK_INACCURACY_HXX
/**
 * @file cmc_t8_adapt_track_inaccuracy.hxx
 */

#include <vector>
#include <memory>
#include <iostream>
#include <unordered_map>

namespace cmc
{

enum TrackingOption {TrackFullInaccuracy = 0, TrackMinimalWorkingInaccuracy};

class InaccuracyContainer
{
public:
    virtual bool CheckInaccuracy() = 0;
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