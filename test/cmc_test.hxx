#ifndef CMC_TEST_HXX
#define CMC_TEST_HXX

namespace cmc
{

constexpr inline int CMC_TEST_SUCCESS = 0;
constexpr inline int CMC_TEST_SKIP = 77;
constexpr inline int CMC_HARD_FAILURE = 99;
constexpr inline int CMC_TEST_FAILURE = 1;

template<typename T>
inline 
void ExpectEQ(T val1, T val2)
{
    if (val1 != val2)
    {
        std::exit(CMC_TEST_FAILURE);
    }
}

inline
void ExpectTrue(const bool expr)
{
    (expr ? void(0) : std::exit(CMC_TEST_FAILURE));
}

}

#endif /* !CMC_TEST_HXX */
