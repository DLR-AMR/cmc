#ifndef CMC_ENDIAN_HXX
#define CMC_ENDIAN_HXX

namespace cmc
{

enum Endian 
{
    /* These macros are GCC macros. Other compilers might not work. */
    Little = __ORDER_LITTLE_ENDIAN__,
    Big = __ORDER_BIG_ENDIAN__,
    Native = __BYTE_ORDER__
};

constexpr bool IsLittleEndian = (Endian::Native == Endian::Little);
constexpr bool IsBigEndian = (Endian::Native == Endian::Big);
constexpr bool IsUnsupportedEndianness = (IsLittleEndian || IsBigEndian);

}

#endif /* !CMC_ENDIAN_HXX */
