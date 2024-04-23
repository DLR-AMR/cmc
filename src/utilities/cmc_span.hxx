#ifndef CMC_SPAN_HXX
#define CMC_SPAN_HXX

#include "cmc_utilities.hxx"

namespace cmc
{

template<class T>
class VectorView
{
public:
    VectorView() = delete;
    VectorView(const T* data, std::size_t size)
    : data_{data}, size_{size}{};

    constexpr VectorView(const VectorView& other) = default;
    constexpr VectorView& operator=(const VectorView& other) = default;
    
    template<typename U> auto operator[](U index) const -> std::enable_if_t<std::is_integral_v<U>, T> {return data_[index];};

    constexpr const T* data() const
    {
        return data_;
    };

    constexpr size_t size() const
    {
        return size_;
    };

    constexpr const T* begin() const
    {
        return data();
    };

    constexpr const T* end() const
    {
        return (data() + size());
    };

    [[nodiscard]] constexpr bool empty() const
    {
        return (size_ == 0);
    };
    
private:
    const T* data_;
    std::size_t size_;
};

}

#endif /* !CMC_SPAN_HXX */
