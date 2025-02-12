#ifndef CMC_VECTOR_VIEW_HXX
#define CMC_VECTOR_VIEW_HXX

#include "cmc_utilities.hxx"

#include <vector>

namespace cmc
{

template<class T>
class VectorView
{
public:
    VectorView() = delete;
    VectorView(const T* data, std::size_t size)
    : data_{data}, size_{size}{};
    VectorView(const std::vector<T>& vector)
    : data_{vector.data()}, size_{vector.size()} {};

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

#endif /* !CMC_VECTOR_VIEW_HXX */
