#ifndef CMC_BYTE_ITERATION_HXX
#define CMC_BYTE_ITERATION_HXX

#include "utilities/cmc_endian.hxx"
#include "utilities/cmc_log_functions.hxx"

namespace cmc
{

constexpr int kByteIterationError = -1;

template<typename T>
static inline
int GetMSBytePosition(const T& value)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return sizeof(value) - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

template<typename T>
static inline
int GetLSBytePosition(const T& value)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return sizeof(value) - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

#if 0
template<int N>
inline
constexpr int GetMSBByteStart([[maybe_unused]] const Serialized<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return N - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

template<int N>
inline
constexpr int GetMSBByteEnd([[maybe_unused]] const Serialized<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return N - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

static inline
void MSBByteIncrement(int& iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        --iterator;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        ++iterator;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
    }
}

template<int N>
inline
bool MSBContinueIteration(const int iterator, [[maybe_unused]] const Serialized<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return iterator >= GetMSBByteEnd(prefix);
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return iterator <= GetMSBByteEnd(prefix);
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return false;
    }
}

template<int N>
inline
constexpr int GetLSBByteStart([[maybe_unused]] const Serialized<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return N - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

template<int N>
inline
constexpr int GetLSBByteEnd([[maybe_unused]] const Serialized<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return N - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

static inline
void LSBByteIncrement(int& iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        ++iterator;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        --iterator;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
    }
}

template<int N>
inline
bool LSBContinueIteration(const int iterator, [[maybe_unused]] const Serialized<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return iterator <= GetLSBByteEnd(prefix);
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return iterator >= GetLSBByteEnd(prefix);
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return false;
    }
}
#endif
template<int N>
inline
constexpr int GetMSBByteStart()
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return N - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

template<int N>
inline
constexpr int GetMSBByteEnd()
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return N - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

static inline
void MSBByteIncrement(int& iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        --iterator;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        ++iterator;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
    }
}

template<int N>
inline
bool MSBContinueIteration(const int iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return iterator >= GetMSBByteEnd<N>();
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return iterator <= GetMSBByteEnd<N>();
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return false;
    }
}


template<int N>
inline
constexpr int GetLSBByteStart()
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return N - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

template<int N>
inline
constexpr int GetLSBByteEnd()
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return N - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kByteIterationError;
    }
}

static inline
void LSBByteIncrement(int& iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        ++iterator;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        --iterator;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
    }
}

template<int N>
inline
bool LSBContinueIteration(const int iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return iterator <= GetLSBByteEnd<N>();
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return iterator >= GetLSBByteEnd<N>();
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return false;
    }
}

};

#endif /* !CMC_BYTE_ITERATION_HXX */
