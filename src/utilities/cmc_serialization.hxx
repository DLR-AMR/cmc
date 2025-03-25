#ifndef CMC_SERIALIZATION_HXX
#define CMC_SERIALIZATION_HXX

#include "utilities/cmc_byte_value.hxx"

#include <array>
#include <type_traits>
#include <cstring>
#include <algorithm>
#include <cstddef>

namespace cmc
{

/**
 * @brief Pushes a value byte-wise in Big-Endian ordering back to a vector
 * 
 * @tparam T The data type of the \a value
 * @param byte_stream The stream to which the value is byte-wise appended
 * @param value The value to be pushed back
 */
template <typename T>
void
PushBackValueToByteStream(std::vector<uint8_t>& byte_stream, const T& value)
{
    /* Serialize the value to a collection of bytes (in big endian order) */
    const std::array<uint8_t, sizeof(T)> serialized_value = SerializeValue(value, Endian::Big);

    /* Append the bytes to the byte stream */
    std::copy_n(serialized_value.begin(), sizeof(T), std::back_insert_iterator(byte_stream));
}

/* A value is extracted from the byte stream starting at position \a pos to which the iterator points to.
 * The serialized value is stored in big endian order.  */
template <typename T, typename Iter>
auto GetValueFromByteStream(Iter pos)
    -> std::enable_if_t<std::is_fundamental_v<T>, T>
{
    /* Get the number of bytes for this data type */
    const size_t type_length = sizeof(T);

    /* Deserialize the given */
    const std::array<uint8_t, sizeof(T)> deserialized_value = DeserializeValue<T>(pos, Endian::Big);

    T value;

    /* Copy the bytes in the correct endianness to the 'value' */
    std::memcpy(static_cast<void*>(&value), static_cast<const void*>(deserialized_value.data()), type_length);
    
    return value;
}

}



#endif /* !CMC_SERIALIZATION_HXX */