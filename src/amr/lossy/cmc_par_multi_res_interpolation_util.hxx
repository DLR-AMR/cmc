#ifndef CMC_AMR_LOSSY_PAR_MULTI_RES_INTERPOLATION_UTIL_HXX
#define CMC_AMR_LOSSY_PAR_MULTI_RES_INTERPOLATION_UTIL_HXX

#include "cmc.hxx"

#include <cstdint>
#include <string>
#include <type_traits>
#include <limits>
#include <utility>
#include <array>
#include <vector>
#include <climits>
#include <algorithm>
#include <numeric>
#include <concepts>
#include <cmath>

namespace cmc::par::lossy::rbf::util
{

template<typename T>
concept FloatType = (std::is_floating_point_v<T> && std::is_arithmetic_v<T>);

constexpr int kNumHexControlPoints = 7;
constexpr int kNumHexPredictionPoints = 7;

constexpr int kNumQuadControlPoints = 5;
constexpr int kNumQuadPredictionPoints = 3;

constexpr int kNumExtendedQuadControlPoints = 9;
constexpr int kNumExtendedQuadPredictionPoints = 3;

constexpr float eps_general = 0.1;
constexpr float eps_gauss = 0.1;
constexpr float eps_mq = 0.1;

template <FloatType T>
constexpr inline
T Eps()
{
    return static_cast<T>(eps_general);
}

template <FloatType T>
constexpr inline
T EpsGauss()
{
    return static_cast<T>(eps_gauss);
}

template <FloatType T>
constexpr inline
T EpsMQ()
{
    return static_cast<T>(eps_mq);
}

template<FloatType T>
constexpr inline 
T
cmc_rbf_inverse_quadratic(const T dist)
{
    return 1.0 / (1.0 + (Eps<T>()  * dist) * (Eps<T>()  * dist));
}

template<FloatType T>
constexpr inline 
T
cmc_rbf_inverse_multiquardic(const T dist)
{
    return 1.0 / std::sqrt(1.0 + (Eps<T>()  * dist) * (Eps<T>()  * dist));
}

template<FloatType T>
constexpr inline 
T
cmc_rbf_multiquardic(const T dist)
{
    return std::sqrt(1.0 + (EpsMQ<T>()  * dist) * (EpsMQ<T>()  * dist));
}

template<FloatType T>
constexpr inline 
T
cmc_rbf_gaussian_eval(const T dist)
{
    return std::exp(-1.0 * (EpsGauss<T>() * dist) * (EpsGauss<T>() * dist));
}


template<FloatType T>
constexpr inline 
T
cmc_rbf_inverse_quadratic(const T dist, const T eps)
{
    return 1.0 / (1.0 + (eps  * dist) * (eps  * dist));
}

template<FloatType T>
constexpr inline 
T
cmc_rbf_inverse_multiquardic(const T dist, const T eps)
{
    return 1.0 / std::sqrt(1.0 + (eps  * dist) * (eps  * dist));
}

template<FloatType T>
constexpr inline 
T
cmc_rbf_multiquardic(const T dist, const T eps)
{
    return std::sqrt(1.0 + (eps  * dist) * (eps  * dist));
}

template<FloatType T>
constexpr inline 
T
cmc_rbf_gaussian_eval(const T dist, const T eps)
{
    return std::exp(-1.0 * (eps * dist) * (eps * dist));
}

template<FloatType T>
constexpr inline 
T
cmc_rbf(const T dist)
{
    constexpr float eps = 0.1;

    //Gauss
    //return std::exp(-1.0 * (eps * dist) * (eps * dist));
    
    //MultiQuardic
    return std::sqrt(1.0 + (eps  * dist) * (eps  * dist));
}

/** Hex Element Operators **/
template <FloatType T>
constexpr
std::array<std::array<T, kNumHexControlPoints>, kNumHexControlPoints>
GetHexMatrixInverse()
{
    const long double x = cmc_rbf<long double>(2.0);
    const long double y = cmc_rbf<long double>(4.0);
    const long double z = cmc_rbf<long double>(std::sqrt(8.0));

    const T a = static_cast<T>((y + 4.0*z + 1.0)/(-6.0*x*x + y + 4.0*z + 1.0));
    const T b = static_cast<T>(x/(6.0*x*x - y - 4.0*z - 1.0));
    const T c = static_cast<T>((-x*x*y + 6.0*x*x*z - 5.0*x*x + y - 4.0*z*z + 2.0*z + 1.0)/(6.0*x*x*y*y - 12.0*x*x*y*z + 12.0*x*x*z - 6.0*x*x - y*y*y - 2.0*y*y*z - y*y + 8.0*y*z*z + y - 8.0*z*z + 2.0*z + 1.0));
    const T d = static_cast<T>((5.0*x*x*y - 6.0*x*x*z + x*x - y*y - 2.0*y*z - y + 4.0*z*z)/(6.0*x*x*y*y - 12.0*x*x*y*z + 12.0*x*x*z - 6.0*x*x - y*y*y - 2.0*y*y*z - y*y + 8.0*y*z*z + y - 8.0*z*z + 2.0*z + 1.0));
    const T e = static_cast<T>((x*x - z)/(-6.0*x*x*y + 12.0*x*x*z - 6.0*x*x + y*y + 2.0*y*z + 2.0*y - 8.0*z*z + 2.0*z + 1.0));

    std::array<std::array<T, kNumHexControlPoints>, kNumHexControlPoints>
    M_inv{{
    {a,b,b,b,b,b,b},
    {b,c,d,e,e,e,e},
    {b,d,c,e,e,e,e},
    {b,e,e,c,d,e,e},
    {b,e,e,d,c,e,e},
    {b,e,e,e,e,c,d},
    {b,e,e,e,e,d,c},
    }};

    return M_inv;
}


template <FloatType T>
constexpr
std::array<std::array<T, kNumHexControlPoints>, kNumHexControlPoints>
GetHexMatrix()
{
    const T x = cmc_rbf<T>(2.0);
    const T y = cmc_rbf<T>(4.0);
    const T z = cmc_rbf<T>(std::sqrt(8));

    std::array<std::array<T, kNumHexControlPoints>, kNumHexControlPoints>
    M{{
    {1.0,x,x,x,x,x,x},
    {x,1.0,y,z,z,z,z},
    {x,y,1.0,z,z,z,z},
    {x,z,z,1.0,y,z,z},
    {x,z,z,y,1.0,z,z},
    {x,z,z,z,z,1.0,y},
    {x,z,z,z,z,y,1.0},
    }};

    return M;
}

template <FloatType T>
constexpr 
std::array<std::array<T, kNumHexControlPoints>, kNumHexPredictionPoints>
GetHexPredictionCoordsEvaluation()
{
    std::array<std::array<T, kNumHexControlPoints>, kNumHexPredictionPoints>
    eval_coords{{
    {cmc_rbf(1.0),            cmc_rbf(3.0),             cmc_rbf(1.0),            cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0)), cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0))},
    {cmc_rbf(1.0),            cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0)), cmc_rbf(3.0),             cmc_rbf(1.0),            cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0))},
    {cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(6.0)),  cmc_rbf(std::sqrt(6.0))},
    {cmc_rbf(1.0),            cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0)), cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0)), cmc_rbf(3.0),             cmc_rbf(1.0)},
    {cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(6.0)),  cmc_rbf(std::sqrt(6.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0))},
    {cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(6.0)),  cmc_rbf(std::sqrt(6.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0))},
    {cmc_rbf(std::sqrt(3.0)), cmc_rbf(std::sqrt(11.0)), cmc_rbf(std::sqrt(3.0)), cmc_rbf(std::sqrt(11.0)), cmc_rbf(std::sqrt(3.0)), cmc_rbf(std::sqrt(11.0)), cmc_rbf(std::sqrt(3.0))}
    }};

    return eval_coords;
}

template <FloatType T>
std::array<T, kNumHexControlPoints>
EvaluateHexWeights(const std::array<std::array<T, kNumHexControlPoints>, kNumHexControlPoints>& M_inv, const std::vector<T>& control_values)
{
    cmc_assert(control_values.size() == static_cast<size_t>(kNumHexControlPoints));

    std::array<T, kNumHexControlPoints> weights;

    /* Compute matrix vector multiplication in order to obtainthe weights */
    int idx = 0;
    for (auto row_iter = M_inv.begin(); row_iter != M_inv.end(); ++row_iter, ++idx)
    {
        weights[idx] = std::inner_product(row_iter->begin(), row_iter->end(), control_values.begin(), static_cast<T>(0.0));
    }

    return weights;
}

template <FloatType T>
std::vector<T>
EvaluateHexPrediction(const std::array<T, kNumHexControlPoints>& weights, const std::array<std::array<T, kNumHexControlPoints>, kNumHexPredictionPoints>& eval_prediction_coords)
{
    std::vector<T> predictions(kNumHexPredictionPoints);

    for (int idx{0}; idx < kNumHexPredictionPoints; ++idx)
    {
        const T prediction = std::inner_product(weights.begin(), weights.end(), eval_prediction_coords[idx].begin(), static_cast<T>(0.0));

        predictions[idx] = prediction;        
    }

    return predictions;
}

template <FloatType T>
std::array<T, kNumHexPredictionPoints>
PerformHexPrediction(const std::array<T, kNumHexControlPoints>& control_values)
{
    /* Get the inverse sytem matrix */
    constexpr std::array<std::array<T, kNumHexControlPoints>, kNumHexControlPoints> M_inv = GetHexMatrixInverse<T>();

    /* Compute the weights for the quad prediction */
    const std::array<T, kNumHexControlPoints> weights{
        M_inv[0][0] * control_values[0] + M_inv[0][1] * control_values[1] + M_inv[0][2] * control_values[2] + M_inv[0][3] * control_values[3] + M_inv[0][4] * control_values[4] + M_inv[0][5] * control_values[5] + M_inv[0][6] * control_values[6],
        M_inv[1][0] * control_values[0] + M_inv[1][1] * control_values[1] + M_inv[1][2] * control_values[2] + M_inv[1][3] * control_values[3] + M_inv[1][4] * control_values[4] + M_inv[1][5] * control_values[5] + M_inv[1][6] * control_values[6],
        M_inv[2][0] * control_values[0] + M_inv[2][1] * control_values[1] + M_inv[2][2] * control_values[2] + M_inv[2][3] * control_values[3] + M_inv[2][4] * control_values[4] + M_inv[2][5] * control_values[5] + M_inv[2][6] * control_values[6],
        M_inv[3][0] * control_values[0] + M_inv[3][1] * control_values[1] + M_inv[3][2] * control_values[2] + M_inv[3][3] * control_values[3] + M_inv[3][4] * control_values[4] + M_inv[3][5] * control_values[5] + M_inv[3][6] * control_values[6],
        M_inv[4][0] * control_values[0] + M_inv[4][1] * control_values[1] + M_inv[4][2] * control_values[2] + M_inv[4][3] * control_values[3] + M_inv[4][4] * control_values[4] + M_inv[4][5] * control_values[5] + M_inv[4][6] * control_values[6],
        M_inv[5][0] * control_values[0] + M_inv[5][1] * control_values[1] + M_inv[5][2] * control_values[2] + M_inv[5][3] * control_values[3] + M_inv[5][4] * control_values[4] + M_inv[5][5] * control_values[5] + M_inv[5][6] * control_values[6],
        M_inv[6][0] * control_values[0] + M_inv[6][1] * control_values[1] + M_inv[6][2] * control_values[2] + M_inv[6][3] * control_values[3] + M_inv[6][4] * control_values[4] + M_inv[6][5] * control_values[5] + M_inv[6][6] * control_values[6]
    };

    /* Get the evaluation coordinates */
    constexpr std::array<std::array<T, kNumHexControlPoints>, kNumHexPredictionPoints> eval_coords = GetHexPredictionCoordsEvaluation<T>();

    /* Perform the prediction */
    const std::array<T, kNumHexPredictionPoints> prediction{
        weights[0] * eval_coords[0][0] + weights[1] * eval_coords[0][1] + weights[2] * eval_coords[0][2] + weights[3] * eval_coords[0][3] + weights[4] * eval_coords[0][4] + weights[5] * eval_coords[0][5] + weights[6] * eval_coords[0][6],
        weights[0] * eval_coords[1][0] + weights[1] * eval_coords[1][1] + weights[2] * eval_coords[1][2] + weights[3] * eval_coords[1][3] + weights[4] * eval_coords[1][4] + weights[5] * eval_coords[1][5] + weights[6] * eval_coords[1][6],
        weights[0] * eval_coords[2][0] + weights[1] * eval_coords[2][1] + weights[2] * eval_coords[2][2] + weights[3] * eval_coords[2][3] + weights[4] * eval_coords[2][4] + weights[5] * eval_coords[2][5] + weights[6] * eval_coords[2][6],
        weights[0] * eval_coords[3][0] + weights[1] * eval_coords[3][1] + weights[2] * eval_coords[3][2] + weights[3] * eval_coords[3][3] + weights[4] * eval_coords[3][4] + weights[5] * eval_coords[3][5] + weights[6] * eval_coords[3][6],
        weights[0] * eval_coords[4][0] + weights[1] * eval_coords[4][1] + weights[2] * eval_coords[4][2] + weights[3] * eval_coords[4][3] + weights[4] * eval_coords[4][4] + weights[5] * eval_coords[4][5] + weights[6] * eval_coords[4][6],
        weights[0] * eval_coords[5][0] + weights[1] * eval_coords[5][1] + weights[2] * eval_coords[5][2] + weights[3] * eval_coords[5][3] + weights[4] * eval_coords[5][4] + weights[5] * eval_coords[5][5] + weights[6] * eval_coords[5][6],
        weights[0] * eval_coords[6][0] + weights[1] * eval_coords[6][1] + weights[2] * eval_coords[6][2] + weights[3] * eval_coords[6][3] + weights[4] * eval_coords[6][4] + weights[5] * eval_coords[6][5] + weights[6] * eval_coords[6][6],
    };

    return prediction;
}

/** END OF Hex Element Operators **/

/** Quad Element Operators **/
template <FloatType T>
constexpr inline 
std::array<std::array<T, kNumQuadControlPoints>, kNumQuadControlPoints>
GetQuadMatrixInverse()
{
    constexpr long double x = cmc_rbf<long double>(2.0);
    constexpr long double y = cmc_rbf<long double>(4.0);
    constexpr long double z = cmc_rbf<long double>(std::sqrt(8.0));

    constexpr T a = static_cast<T>((-y - 2.0*z - 1.0) / (4.0*x*x - y - 2.0*z - 1.0));
    constexpr T b = static_cast<T>(x / (4.0*x*x - y - 2.0*z - 1));
    constexpr T c = static_cast<T>((-x*x*y + 4.0*x*x*z - 3.0*x*x + y - 2.0*z*z + 1.0) / (4.0*x*x*y*y - 8.0*x*x*y*z + 8.0*x*x*z - 4.0*x*x - y*y*y - y*y + 4.0*y*z*z + y - 4.0*z*z + 1.0));
    constexpr T d = static_cast<T>((3.0*x*x*y - 4.0*x*x*z + x*x - y*y - y + 2.0*z*z)  / (4.0*x*x*y*y - 8.0*x*x*y*z + 8.0*x*x*z - 4.0*x*x - y*y*y - y*y + 4.0*y*z*z + y - 4.0*z*z + 1.0));
    constexpr T e = static_cast<T>((-x*x + z) / (4.0*x*x*y - 8.0*x*x*z + 4.0*x*x - y*y - 2.0*y + 4.0*z*z - 1.0));

    constexpr std::array<std::array<T, kNumQuadControlPoints>, kNumQuadControlPoints>
        M_inv{{
        {a,b,b,b,b},
        {b,c,d,e,e},
        {b,d,c,e,e},
        {b,e,e,c,d},
        {b,e,e,d,c},
        }};

    return M_inv;
}

template <FloatType T>
constexpr inline
std::array<std::array<T, kNumQuadControlPoints>, kNumQuadControlPoints>
GetQuadMatrix()
{
    constexpr T x = cmc_rbf<T>(2.0);
    constexpr T y = cmc_rbf<T>(4.0);
    constexpr T z = cmc_rbf<T>(std::sqrt(8));

    constexpr std::array<std::array<T, kNumQuadControlPoints>, kNumQuadControlPoints>
        M{{
        {1.0, x, x, x, x}, 
        {x, 1.0, y, z, z},
        {x, y, 1.0, z, z},
        {x, z, z, 1.0, y},
        {x, z, z, y, 1.0}
        }};

    return M;
}

template <FloatType T>
constexpr inline
std::array<std::array<T, kNumQuadControlPoints>, kNumQuadPredictionPoints>
GetQuadPredictionCoordsEvaluation()
{
    constexpr std::array<std::array<T, kNumQuadControlPoints>, kNumQuadPredictionPoints>
        eval_coords{{
        {cmc_rbf(1.0),            cmc_rbf(3.0),             cmc_rbf(1.0),            cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0))},
        {cmc_rbf(1.0),            cmc_rbf(std::sqrt(5.0)),  cmc_rbf(std::sqrt(5.0)), cmc_rbf(3.0),             cmc_rbf(1.0)           },
        {cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0)), cmc_rbf(std::sqrt(10.0)), cmc_rbf(std::sqrt(2.0))},
        }};

    return eval_coords;
}

template <FloatType T>
inline std::array<T, kNumQuadControlPoints>
EvaluateQuadWeights(const std::array<std::array<T, kNumQuadControlPoints>, kNumQuadControlPoints>& M_inv, const std::array<T, kNumQuadControlPoints>& control_values)
{
    std::array<T, kNumQuadControlPoints> weights;

    /* Compute matrix vector multiplication in order to obtainthe weights */
    int idx = 0;
    for (auto row_iter = M_inv.begin(); row_iter != M_inv.end(); ++row_iter, ++idx)
    {
        weights[idx] = std::inner_product(row_iter->begin(), row_iter->end(), control_values.begin(), static_cast<T>(0.0));
    }

    return weights;
}

template <FloatType T>
inline std::vector<T>
EvaluateQuadPrediction(const std::array<T, kNumQuadControlPoints>& weights, const std::array<std::array<T, kNumQuadControlPoints>, kNumQuadPredictionPoints>& eval_prediction_coords)
{
    std::vector<T> predictions(kNumQuadPredictionPoints);

    for (int idx{0}; idx < kNumQuadPredictionPoints; ++idx)
    {
        const T prediction = std::inner_product(weights.begin(), weights.end(), eval_prediction_coords[idx].begin(), static_cast<T>(0.0));

        predictions[idx] = prediction;        
    }

    return predictions;
}

template <FloatType T>
std::array<T, kNumQuadPredictionPoints>
PerformQuadPrediction(const std::array<T, kNumQuadControlPoints>& control_values)
{
    /* Get the inverse sytem matrix */
    constexpr std::array<std::array<T, kNumQuadControlPoints>, kNumQuadControlPoints> M_inv = GetQuadMatrixInverse<T>();

    /* Compute the weights for the quad prediction */
    const std::array<T, kNumQuadControlPoints> weights{
        M_inv[0][0] * control_values[0] + M_inv[0][1] * control_values[1] + M_inv[0][2] * control_values[2] + M_inv[0][3] * control_values[3] + M_inv[0][4] * control_values[4],
        M_inv[1][0] * control_values[0] + M_inv[1][1] * control_values[1] + M_inv[1][2] * control_values[2] + M_inv[1][3] * control_values[3] + M_inv[1][4] * control_values[4],
        M_inv[2][0] * control_values[0] + M_inv[2][1] * control_values[1] + M_inv[2][2] * control_values[2] + M_inv[2][3] * control_values[3] + M_inv[2][4] * control_values[4],
        M_inv[3][0] * control_values[0] + M_inv[3][1] * control_values[1] + M_inv[3][2] * control_values[2] + M_inv[3][3] * control_values[3] + M_inv[3][4] * control_values[4],
        M_inv[4][0] * control_values[0] + M_inv[4][1] * control_values[1] + M_inv[4][2] * control_values[2] + M_inv[4][3] * control_values[3] + M_inv[4][4] * control_values[4]
    };

    /* Get the evaluation coordinates */
    constexpr std::array<std::array<T, kNumQuadControlPoints>, kNumQuadPredictionPoints> eval_coords = GetQuadPredictionCoordsEvaluation<T>();

    /* Perform the prediction */
    const std::array<T, kNumQuadPredictionPoints> prediction{
        weights[0] * eval_coords[0][0] + weights[1] * eval_coords[0][1] + weights[2] * eval_coords[0][2] + weights[3] * eval_coords[0][3] + weights[4] * eval_coords[0][4],
        weights[0] * eval_coords[1][0] + weights[1] * eval_coords[1][1] + weights[2] * eval_coords[1][2] + weights[3] * eval_coords[1][3] + weights[4] * eval_coords[1][4],
        weights[0] * eval_coords[2][0] + weights[1] * eval_coords[2][1] + weights[2] * eval_coords[2][2] + weights[3] * eval_coords[2][3] + weights[4] * eval_coords[2][4]
    };

    return prediction;
}


/** END OF Quad Element Operators **/

/** Start Extended Quad Element Operators */

template <FloatType T>
constexpr inline 
std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadControlPoints>
GetExtendedQuadMatrixInverse()
{   
    constexpr long double a = cmc_rbf<T>(std::sqrt(static_cast<T>(8.0)));
    constexpr long double b = cmc_rbf<T>(2.0);
    constexpr long double c = cmc_rbf<T>(4.0);
    constexpr long double d = cmc_rbf<T>(std::sqrt(static_cast<T>(20.0)));
    constexpr long double e = cmc_rbf<T>(std::sqrt(static_cast<T>(32.0)));

    constexpr long double f = (-4.0*a*c - 2.0*a*e - 2.0*a + 4.0*b*b + 8.0*b*d - 2.0*c*c - c*e - 3.0*c + 4.0*d*d - e - 1.0)
    /
    (8.0*a*a*a + 4.0*a*a*c + 4.0*a*a - 16.0*a*b*b - 16.0*a*b*d - 4.0*a*c - 2.0*a*e - 2.0*a + 8.0*b*b*c + 4.0*b*b*e + 8.0*b*b + 8.0*b*d - 2.0*c*c - c*e - 3.0*c + 4.0*d*d - e - 1.0);

    constexpr long double g = (2.0*a*a + a*c + a - 2.0*b*b - 2.0*b*d)
    /
    (8.0*a*a*a + 4.0*a*a*c + 4.0*a*a - 16.0*a*b*b - 16.0*a*b*d - 4.0*a*c - 2.0*a*e - 2.0*a + 8.0*b*b*c + 4.0*b*b*e + 8.0*b*b + 8.0*b*d - 2.0*c*c - c*e - 3.0*c + 4.0*d*d - e - 1.0);

    constexpr long double h = (-2.0*a*b - 2.0*a*d + 2.0*b*c + b*e + b)
    /
    (8.0*a*a*a + 4.0*a*a*c + 4.0*a*a - 16.0*a*b*b - 16.0*a*b*d - 4.0*a*c - 2.0*a*e - 2.0*a + 8.0*b*b*c + 4.0*b*b*e + 8.0*b*b + 8.0*b*d - 2.0*c*c - c*e - 3.0*c + 4.0*d*d - e - 1.0);

    constexpr long double i = (-4.0*a*a*a*b*b
    + 8.0*a*a*a*b*d
    + 8.0*a*a*a*c*c
    - 2.0*a*a*a*c*e
    - 14.0*a*a*a*c
    - 4.0*a*a*a*d*d
    + 2.0*a*a*a*e
    + 6.0*a*a*a
    - 2.0*a*a*b*b*c
    - 2.0*a*a*b*b
    + 4.0*a*a*b*c*d
    + 4.0*a*a*b*d
    + 4.0*a*a*c*c*c
    - a*a*c*c*e
    - 3.0*a*a*c*c
    - 2.0*a*a*c*d*d
    - 4.0*a*a*c
    - 2.0*a*a*d*d
    + a*a*e
    + 3.0*a*a
    + 8.0*a*b*b*b*b
    - 8.0*a*b*b*b*d
    - 16.0*a*b*b*c*c
    + 4.0*a*b*b*c*e
    + 28.0*a*b*b*c
    - 8.0*a*b*b*d*d
    - 2.0*a*b*b*e
    - 10.0*a*b*b
    - 16.0*a*b*c*c*d
    + 4.0*a*b*c*d*e
    + 28.0*a*b*c*d
    + 8.0*a*b*d*d*d
    - 8.0*a*b*d*e
    - 16.0*a*b*d
    - 4.0*a*c*c*c
    + 4.0*a*c*c
    + 2.0*a*c*e
    + 2.0*a*c
    + 2.0*a*d*d*e
    + 2.0*a*d*d
    - 2.0*a*e
    - 2.0*a
    - 4.0*b*b*b*b*e
    - 6.0*b*b*b*b
    + 8.0*b*b*b*d*e
    + 8.0*b*b*b*d
    + 8.0*b*b*c*c*c
    - 4.0*b*b*c*c
    - 4.0*b*b*c*e
    - 10.0*b*b*c
    - 4.0*b*b*d*d*e
    + 6.0*b*b*e
    + 8.0*b*b
    + 8.0*b*c*c*d
    - 4.0*b*c*d*e
    - 16.0*b*c*d
    + 4.0*b*d
    - 2.0*c*c*c*c
    + 4.0*c*c*d*d
    + c*c*e
    + 3.0*c*c
    - 6.0*c*d*d
    - 2.0*d*d*d*d
    + 2.0*d*d*e
    + 4.0*d*d
    - e
    - 1.0)
    /
    (32.0*a*a*a*b*b*c
    - 16.0*a*a*a*b*b*e
    - 16.0*a*a*a*b*b
    - 64.0*a*a*a*b*c*d
    + 32.0*a*a*a*b*d*e
    + 32.0*a*a*a*b*d
    - 16.0*a*a*a*c*c*e
    + 16.0*a*a*a*c*c
    + 32.0*a*a*a*c*d*d
    + 8.0*a*a*a*c*e*e
    + 16.0*a*a*a*c*e
    - 24.0*a*a*a*c
    - 16.0*a*a*a*d*d*e
    - 16.0*a*a*a*d*d
    - 8.0*a*a*a*e*e
    + 8.0*a*a*a
    + 16.0*a*a*b*b*c*c
    - 8.0*a*a*b*b*c*e
    + 8.0*a*a*b*b*c
    - 8.0*a*a*b*b*e
    - 8.0*a*a*b*b
    - 32.0*a*a*b*c*c*d
    + 16.0*a*a*b*c*d*e
    - 16.0*a*a*b*c*d
    + 16.0*a*a*b*d*e
    + 16.0*a*a*b*d
    - 8.0*a*a*c*c*c*e
    + 8.0*a*a*c*c*c
    + 16.0*a*a*c*c*d*d
    + 4.0*a*a*c*c*e*e
    - 4.0*a*a*c*c
    - 8.0*a*a*c*d*d*e
    + 8.0*a*a*c*d*d
    + 8.0*a*a*c*e
    - 8.0*a*a*c
    - 8.0*a*a*d*d*e
    - 8.0*a*a*d*d
    - 4.0*a*a*e*e
    + 4.0*a*a
    - 64.0*a*b*b*b*b*c
    + 32.0*a*b*b*b*b*e
    + 32.0*a*b*b*b*b
    + 64.0*a*b*b*b*c*d
    - 32.0*a*b*b*b*d*e
    - 32.0*a*b*b*b*d
    + 32.0*a*b*b*c*c*e
    - 48.0*a*b*b*c*c
    + 64.0*a*b*b*c*d*d
    - 16.0*a*b*b*c*e*e
    - 32.0*a*b*b*c*e
    + 48.0*a*b*b*c
    - 32.0*a*b*b*d*d*e
    - 32.0*a*b*b*d*d
    + 20.0*a*b*b*e*e
    + 8.0*a*b*b*e
    - 12.0*a*b*b
    + 32.0*a*b*c*c*d*e
    - 64.0*a*b*c*d*d*d
    - 16.0*a*b*c*d*e*e
    - 32.0*a*b*c*d*e
    + 48.0*a*b*c*d
    + 32.0*a*b*d*d*d*e
    + 32.0*a*b*d*d*d
    + 8.0*a*b*d*e*e
    - 16.0*a*b*d*e
    - 24.0*a*b*d
    + 8.0*a*c*c*c*e
    - 8.0*a*c*c*c
    - 16.0*a*c*c*d*d
    - 8.0*a*c*c*e
    + 8.0*a*c*c
    - 2.0*a*c*e*e*e
    - 2.0*a*c*e*e
    + 2.0*a*c*e
    + 2.0*a*c
    + 4.0*a*d*d*e*e
    + 8.0*a*d*d*e
    + 4.0*a*d*d
    + 2.0*a*e*e*e
    + 2.0*a*e*e
    - 2.0*a*e
    - 2.0*a
    + 32.0*b*b*b*b*c*c
    + 16.0*b*b*b*b*c
    - 8.0*b*b*b*b*e*e
    - 24.0*b*b*b*b*e
    - 16.0*b*b*b*b
    - 64.0*b*b*b*c*c*d
    + 16.0*b*b*b*d*e*e
    + 32.0*b*b*b*d*e
    + 16.0*b*b*b*d
    - 16.0*b*b*c*c*c*e
    + 8.0*b*b*c*c*c
    + 32.0*b*b*c*c*d*d
    + 8.0*b*b*c*c*e
    - 16.0*b*b*c*c
    - 32.0*b*b*c*d*d
    + 4.0*b*b*c*e*e*e
    + 10.0*b*b*c*e*e
    + 8.0*b*b*c*e
    - 14.0*b*b*c
    - 8.0*b*b*d*d*e*e
    + 8.0*b*b*d*d
    - 4.0*b*b*e*e*e
    - 6.0*b*b*e*e
    + 8.0*b*b*e
    + 10.0*b*b
    + 16.0*b*c*c*c*d
    - 16.0*b*c*c*d*e
    + 32.0*b*c*c*d
    + 4.0*b*c*d*e*e
    + 8.0*b*c*d*e
    - 28.0*b*c*d
    - 12.0*b*d*e*e
    - 8.0*b*d*e
    + 4.0*b*d
    + 4.0*c*c*c*c*e
    - 4.0*c*c*c*c
    - 8.0*c*c*c*d*d
    - 8.0*c*c*d*d*e
    - c*c*e*e*e
    - c*c*e*e
    - 3.0*c*c*e
    + 5.0*c*c
    + 16.0*c*d*d*d*d
    + 6.0*c*d*d*e*e
    + 12.0*c*d*d*e
    - 10.0*c*d*d
    - 8.0*d*d*d*d*e
    - 8.0*d*d*d*d
    - 2.0*d*d*e*e
    + 4.0*d*d*e
    + 6.0*d*d
    + e*e*e
    + e*e
    - e
    - 1.0);

    constexpr long double j = (4.0*a*a*a*b - 4.0*a*a*a*d
    + 2.0*a*a*b*c + 2.0*a*a*b
    - 2.0*a*a*c*d - 2.0*a*a*d
    - 10.0*a*b*b*b + 4.0*a*b*b*d
    + a*b*c*e - 3.0*a*b*c
    + 6.0*a*b*d*d - 2.0*a*b*e
    + 2.0*a*c*d + a*d*e + a*d
    + 4.0*b*b*b*c + 2.0*b*b*b*e + 5.0*b*b*b
    - 4.0*b*b*c*d - 2.0*b*b*d*e - b*b*d
    - b*c*c - b*c*e - b*c
    - 3.0*b*d*d - b
    + c*c*d + 2.0*c*d
    - d*d*d + d*e)
    /
    (16.0*a*a*a*b*b - 32.0*a*a*a*b*d
    - 8.0*a*a*a*c*e + 8.0*a*a*a*c
    + 16.0*a*a*a*d*d + 8.0*a*a*a*e - 8.0*a*a*a
    + 8.0*a*a*b*b*c + 8.0*a*a*b*b
    - 16.0*a*a*b*c*d - 16.0*a*a*b*d
    - 4.0*a*a*c*c*e + 4.0*a*a*c*c
    + 8.0*a*a*c*d*d + 8.0*a*a*d*d
    + 4.0*a*a*e - 4.0*a*a
    - 32.0*a*b*b*b*b + 32.0*a*b*b*b*d
    + 16.0*a*b*b*c*e - 24.0*a*b*b*c
    + 32.0*a*b*b*d*d - 20.0*a*b*b*e
    + 12.0*a*b*b
    + 16.0*a*b*c*d*e - 32.0*a*b*d*d*d
    - 8.0*a*b*d*e + 24.0*a*b*d
    + 4.0*a*c*c*e - 4.0*a*c*c
    - 8.0*a*c*d*d + 2.0*a*c*e*e
    - 4.0*a*c*e + 2.0*a*c
    - 4.0*a*d*d*e - 4.0*a*d*d
    - 2.0*a*e*e + 2.0*a
    + 16.0*b*b*b*b*c + 8.0*b*b*b*b*e
    + 16.0*b*b*b*b
    - 32.0*b*b*b*c*d - 16.0*b*b*b*d*e
    - 16.0*b*b*b*d
    - 8.0*b*b*c*c*e + 4.0*b*b*c*c
    + 16.0*b*b*c*d*d
    - 4.0*b*b*c*e*e + 2.0*b*b*c*e
    - 6.0*b*b*c
    + 8.0*b*b*d*d*e - 8.0*b*b*d*d
    + 4.0*b*b*e*e + 2.0*b*b*e
    - 10.0*b*b
    + 8.0*b*c*c*d - 4.0*b*c*d*e
    + 20.0*b*c*d + 12.0*b*d*e
    - 4.0*b*d
    + 2.0*c*c*c*e - 2.0*c*c*c
    - 4.0*c*c*d*d + c*c*e*e
    - c*c - 6.0*c*d*d*e
    - 2.0*c*d*d - 2.0*c*e + 2.0*c
    + 8.0*d*d*d*d + 2.0*d*d*e
    - 6.0*d*d - e*e + 1.0);

    constexpr long double l = (-4.0*a*a*a*b + 4.0*a*a*a*d
    + 2.0*a*a*c*d - 2.0*a*a*b*c
    - 2.0*a*a*b + 2.0*a*a*d
    + 6.0*a*b*b*b + 4.0*a*b*b*d
    + a*b*c*e + a*b*c - 10.0*a*b*d*d
    + 2.0*a*b - 2.0*a*c*d - a*d*e - a*d
    - 4.0*b*b*b*c - 2.0*b*b*b*e - 3.0*b*b*b
    + 4.0*b*b*c*d + 2.0*b*b*d*e - b*b*d
    + b*c*c + 2.0*b*c + b*d*d + b*e
    - c*c*d - c*d*e - c*d
    + 3.0*d*d*d - d)
    /
    (16.0*a*a*a*b*b - 32.0*a*a*a*b*d
    - 8.0*a*a*a*c*e + 8.0*a*a*a*c
    + 16.0*a*a*a*d*d + 8.0*a*a*a*e - 8.0*a*a*a
    + 8.0*a*a*b*b*c + 8.0*a*a*b*b
    - 16.0*a*a*b*c*d - 16.0*a*a*b*d
    - 4.0*a*a*c*c*e + 4.0*a*a*c*c
    + 8.0*a*a*c*d*d + 8.0*a*a*d*d
    + 4.0*a*a*e - 4.0*a*a
    - 32.0*a*b*b*b*b + 32.0*a*b*b*b*d
    + 16.0*a*b*b*c*e - 24.0*a*b*b*c
    + 32.0*a*b*b*d*d - 20.0*a*b*b*e
    + 12.0*a*b*b
    + 16.0*a*b*c*d*e - 32.0*a*b*d*d*d
    - 8.0*a*b*d*e + 24.0*a*b*d
    + 4.0*a*c*c*e - 4.0*a*c*c
    - 8.0*a*c*d*d + 2.0*a*c*e*e
    - 4.0*a*c*e + 2.0*a*c
    - 4.0*a*d*d*e - 4.0*a*d*d
    - 2.0*a*e*e + 2.0*a
    + 16.0*b*b*b*b*c + 8.0*b*b*b*b*e
    + 16.0*b*b*b*b
    - 32.0*b*b*b*c*d - 16.0*b*b*b*d*e
    - 16.0*b*b*b*d
    - 8.0*b*b*c*c*e + 4.0*b*b*c*c
    + 16.0*b*b*c*d*d
    - 4.0*b*b*c*e*e + 2.0*b*b*c*e
    - 6.0*b*b*c
    + 8.0*b*b*d*d*e - 8.0*b*b*d*d
    + 4.0*b*b*e*e + 2.0*b*b*e
    - 10.0*b*b
    + 8.0*b*c*c*d - 4.0*b*c*d*e
    + 20.0*b*c*d + 12.0*b*d*e
    - 4.0*b*d
    + 2.0*c*c*c*e - 2.0*c*c*c
    - 4.0*c*c*d*d + c*c*e*e
    - c*c - 6.0*c*d*d*e
    - 2.0*c*d*d - 2.0*c*e + 2.0*c
    + 8.0*d*d*d*d + 2.0*d*d*e
    - 6.0*d*d - e*e + 1.0);

    constexpr long double m = (2.0*a*a*a + a*a*c + a*a
    - 4.0*a*b*b - 4.0*a*b*d - 2.0*a*c
    + 4.0*b*b*c + b*b + 2.0*b*d
    - c*c - c + d*d)
    /
    (16.0*a*a*a*c - 8.0*a*a*a*e - 8.0*a*a*a
    + 8.0*a*a*c*c - 4.0*a*a*c*e
    + 4.0*a*a*c - 4.0*a*a*e - 4.0*a*a
    - 32.0*a*b*b*c + 16.0*a*b*b*e
    + 16.0*a*b*b - 32.0*a*b*c*d
    + 16.0*a*b*d*e + 16.0*a*b*d
    - 8.0*a*c*c + 2.0*a*e*e + 4.0*a*e + 2.0*a
    + 16.0*b*b*c*c + 8.0*b*b*c
    - 4.0*b*b*e*e - 12.0*b*b*e - 8.0*b*b
    + 16.0*b*c*d - 8.0*b*d*e - 8.0*b*d
    - 4.0*c*c*c - 4.0*c*c + 8.0*c*d*d
    + c*e*e + 2.0*c*e + c
    - 4.0*d*d*e - 4.0*d*d
    + e*e + 2.0*e + 1.0);


    constexpr long double n = (-4.0*a*a*a*b*b
    + 8.0*a*a*a*b*d
    - 8.0*a*a*a*c*c
    + 6.0*a*a*a*c*e
    + 10.0*a*a*a*c
    - 4.0*a*a*a*d*d
    - 6.0*a*a*a*e
    - 2.0*a*a*a
    - 2.0*a*a*b*b*c
    - 2.0*a*a*b*b
    + 4.0*a*a*b*c*d
    + 4.0*a*a*b*d
    - 4.0*a*a*c*c*c
    + 3.0*a*a*c*c*e
    + a*a*c*c
    - 2.0*a*a*c*d*d
    + 4.0*a*a*c
    - 2.0*a*a*d*d
    - 3.0*a*a*e
    - a*a
    + 8.0*a*b*b*b*b
    - 8.0*a*b*b*b*d
    + 16.0*a*b*b*c*c
    - 12.0*a*b*b*c*e
    - 20.0*a*b*b*c
    - 8.0*a*b*b*d*d
    + 14.0*a*b*b*e
    + 6.0*a*b*b
    + 16.0*a*b*c*c*d
    - 12.0*a*b*c*d*e
    - 20.0*a*b*c*d
    + 8.0*a*b*d*d*d
    + 8.0*a*b*d*e
    + 4.0*a*c*c*c
    - 4.0*a*c*c
    - 2.0*a*c*e*e
    - 2.0*a*c*e
    + 2.0*a*d*d*e
    + 2.0*a*d*d
    + 2.0*a*e*e
    + 2.0*a*e
    - 4.0*b*b*b*b*e
    - 6.0*b*b*b*b
    + 8.0*b*b*b*d*e
    + 8.0*b*b*b*d
    - 8.0*b*b*c*c*c
    + 4.0*b*b*c*c
    + 4.0*b*b*c*e*e
    + 8.0*b*b*c*e
    + 6.0*b*b*c
    - 4.0*b*b*d*d*e
    - 4.0*b*b*e*e
    - 6.0*b*b*e
    - 8.0*b*c*c*d
    + 4.0*b*c*d*e
    + 8.0*b*c*d
    - 8.0*b*d*e
    - 4.0*b*d
    + 2.0*c*c*c*c
    - 4.0*c*c*d*d
    - c*c*e*e
    - c*c*e
    - 2.0*c*c
    + 4.0*c*d*d*e
    + 6.0*c*d*d
    - 2.0*d*d*d*d
    - 2.0*d*d*e
    + e*e
    + e)
    /
    (32.0*a*a*a*b*b*c
    - 16.0*a*a*a*b*b*e
    - 16.0*a*a*a*b*b
    - 64.0*a*a*a*b*c*d
    + 32.0*a*a*a*b*d*e
    + 32.0*a*a*a*b*d
    - 16.0*a*a*a*c*c*e
    + 16.0*a*a*a*c*c
    + 32.0*a*a*a*c*d*d
    + 8.0*a*a*a*c*e*e
    + 16.0*a*a*a*c*e
    - 24.0*a*a*a*c
    - 16.0*a*a*a*d*d*e
    - 16.0*a*a*a*d*d
    - 8.0*a*a*a*e*e
    + 8.0*a*a*a
    + 16.0*a*a*b*b*c*c
    - 8.0*a*a*b*b*c*e
    + 8.0*a*a*b*b*c
    - 8.0*a*a*b*b*e
    - 8.0*a*a*b*b
    - 32.0*a*a*b*c*c*d
    + 16.0*a*a*b*c*d*e
    - 16.0*a*a*b*c*d
    + 16.0*a*a*b*d*e
    + 16.0*a*a*b*d
    - 8.0*a*a*c*c*c*e
    + 8.0*a*a*c*c*c
    + 16.0*a*a*c*c*d*d
    + 4.0*a*a*c*c*e*e
    - 4.0*a*a*c*c
    - 8.0*a*a*c*d*d*e
    + 8.0*a*a*c*d*d
    + 8.0*a*a*c*e
    - 8.0*a*a*c
    - 8.0*a*a*d*d*e
    - 8.0*a*a*d*d
    - 4.0*a*a*e*e
    + 4.0*a*a
    - 64.0*a*b*b*b*b*c
    + 32.0*a*b*b*b*b*e
    + 32.0*a*b*b*b*b
    + 64.0*a*b*b*b*c*d
    - 32.0*a*b*b*b*d*e
    - 32.0*a*b*b*b*d
    + 32.0*a*b*b*c*c*e
    - 48.0*a*b*b*c*c
    + 64.0*a*b*b*c*d*d
    - 16.0*a*b*b*c*e*e
    - 32.0*a*b*b*c*e
    + 48.0*a*b*b*c
    - 32.0*a*b*b*d*d*e
    - 32.0*a*b*b*d*d
    + 20.0*a*b*b*e*e
    + 8.0*a*b*b*e
    - 12.0*a*b*b
    + 32.0*a*b*c*c*d*e
    - 64.0*a*b*c*d*d*d
    - 16.0*a*b*c*d*e*e
    - 32.0*a*b*c*d*e
    + 48.0*a*b*c*d
    + 32.0*a*b*d*d*d*e
    + 32.0*a*b*d*d*d
    + 8.0*a*b*d*e*e
    - 16.0*a*b*d*e
    - 24.0*a*b*d
    + 8.0*a*c*c*c*e
    - 8.0*a*c*c*c
    - 16.0*a*c*c*d*d
    - 8.0*a*c*c*e
    + 8.0*a*c*c
    - 2.0*a*c*e*e*e
    - 2.0*a*c*e*e
    + 2.0*a*c*e
    + 2.0*a*c
    + 4.0*a*d*d*e*e
    + 8.0*a*d*d*e
    + 4.0*a*d*d
    + 2.0*a*e*e*e
    + 2.0*a*e*e
    - 2.0*a*e
    - 2.0*a
    + 32.0*b*b*b*b*c*c
    + 16.0*b*b*b*b*c
    - 8.0*b*b*b*b*e*e
    - 24.0*b*b*b*b*e
    - 16.0*b*b*b*b
    - 64.0*b*b*b*c*c*d
    + 16.0*b*b*b*d*e*e
    + 32.0*b*b*b*d*e
    + 16.0*b*b*b*d
    - 16.0*b*b*c*c*c*e
    + 8.0*b*b*c*c*c
    + 32.0*b*b*c*c*d*d
    + 8.0*b*b*c*c*e
    - 16.0*b*b*c*c
    - 32.0*b*b*c*d*d
    + 4.0*b*b*c*e*e*e
    + 10.0*b*b*c*e*e
    + 8.0*b*b*c*e
    - 14.0*b*b*c
    - 8.0*b*b*d*d*e*e
    + 8.0*b*b*d*d
    - 4.0*b*b*e*e*e
    - 6.0*b*b*e*e
    + 8.0*b*b*e
    + 10.0*b*b
    + 16.0*b*c*c*c*d
    - 16.0*b*c*c*d*e
    + 32.0*b*c*c*d
    + 4.0*b*c*d*e*e
    + 8.0*b*c*d*e
    - 28.0*b*c*d
    - 12.0*b*d*e*e
    - 8.0*b*d*e
    + 4.0*b*d
    + 4.0*c*c*c*c*e
    - 4.0*c*c*c*c
    - 8.0*c*c*c*d*d
    - 8.0*c*c*d*d*e
    - c*c*e*e*e
    - c*c*e*e
    - 3.0*c*c*e
    + 5.0*c*c
    + 16.0*c*d*d*d*d
    + 6.0*c*d*d*e*e
    + 12.0*c*d*d*e
    - 10.0*c*d*d
    - 8.0*d*d*d*d*e
    - 8.0*d*d*d*d
    - 2.0*d*d*e*e
    + 4.0*d*d*e
    + 6.0*d*d
    + e*e*e
    + e*e
    - e
    - 1.0);


    constexpr long double p = (8.0*a*a*a*a*e - 8.0*a*a*a*a
    - 4.0*a*a*b*b*c - 16.0*a*a*b*b*e + 12.0*a*a*b*b
    + 8.0*a*a*b*c*d - 16.0*a*a*b*d*e + 24.0*a*a*b*d
    - 4.0*a*a*c*d*d - 8.0*a*a*c*e + 8.0*a*a*c
    - 4.0*a*a*d*d - 2.0*a*a*e*e - 4.0*a*a*e + 6.0*a*a
    + 8.0*a*b*b*b*b - 8.0*a*b*b*b*d
    + 12.0*a*b*b*c*e - 12.0*a*b*b*c
    - 8.0*a*b*b*d*d + 4.0*a*b*b*e*e
    + 16.0*a*b*b*e - 20.0*a*b*b
    + 4.0*a*b*c*d*e - 4.0*a*b*c*d
    + 8.0*a*b*d*d*d + 20.0*a*b*d*e - 20.0*a*b*d
    + 4.0*a*d*d*e - 4.0*a*d*d
    - 4.0*b*b*b*b*c - 2.0*b*b*b*b*e - 4.0*b*b*b*b
    + 8.0*b*b*b*c*d + 4.0*b*b*b*d*e + 4.0*b*b*b*d
    - 2.0*b*b*c*c*e + 4.0*b*b*c*c
    - 4.0*b*b*c*d*d - b*b*c*e*e
    - 6.0*b*b*c*e + 11.0*b*b*c
    - 2.0*b*b*d*d*e + 2.0*b*b*d*d
    - 3.0*b*b*e*e - 2.0*b*b*e + 7.0*b*b
    - 4.0*b*c*c*d - 4.0*b*c*d*e
    - 4.0*b*c*d - 8.0*b*d*e + 4.0*b*d
    + 2.0*c*c*d*d + 2.0*c*c*e - 2.0*c*c
    + 4.0*c*d*d + c*e*e + 2.0*c*e - 3.0*c
    - 2.0*d*d*d*d - 2.0*d*d*e + 4.0*d*d
    + e*e - 1.0)
    /
    (32.0*a*a*a*a*b*b - 64.0*a*a*a*a*b*d
    - 16.0*a*a*a*a*c*e + 16.0*a*a*a*a*c
    + 32.0*a*a*a*a*d*d + 16.0*a*a*a*a*e - 16.0*a*a*a*a
    - 64.0*a*a*b*b*b*b + 64.0*a*a*b*b*b*d
    - 8.0*a*a*b*b*c*c + 32.0*a*a*b*b*c*e
    - 64.0*a*a*b*b*c + 64.0*a*a*b*b*d*d
    - 40.0*a*a*b*b*e + 16.0*a*a*b*b
    + 16.0*a*a*b*c*c*d + 32.0*a*a*b*c*d*e
    + 32.0*a*a*b*c*d - 64.0*a*a*b*d*d*d
    - 16.0*a*a*b*d*e + 64.0*a*a*b*d
    + 4.0*a*a*c*c*c*e - 4.0*a*a*c*c*c
    - 8.0*a*a*c*c*d*d + 12.0*a*a*c*c*e
    - 12.0*a*a*c*c - 32.0*a*a*c*d*d
    + 4.0*a*a*c*e*e - 12.0*a*a*c*e
    + 8.0*a*a*c - 8.0*a*a*d*d*e
    - 16.0*a*a*d*d - 4.0*a*a*e*e
    - 4.0*a*a*e + 8.0*a*a
    + 64.0*a*b*b*b*b*c + 16.0*a*b*b*b*b*e
    + 64.0*a*b*b*b*b - 96.0*a*b*b*b*c*d
    - 32.0*a*b*b*b*d*e - 64.0*a*b*b*b*d
    - 32.0*a*b*b*c*c*e + 32.0*a*b*b*c*c
    - 8.0*a*b*b*c*e*e + 8.0*a*b*b*c*e
    + 16.0*a*b*b*d*d*e - 48.0*a*b*b*d*d
    + 8.0*a*b*b*e*e + 24.0*a*b*b*e
    - 32.0*a*b*b - 16.0*a*b*c*c*d*e
    + 16.0*a*b*c*c*d + 32.0*a*b*c*d*d*d
    - 16.0*a*b*c*d*e + 16.0*a*b*c*d
    + 32.0*a*b*d*d*d + 32.0*a*b*d*e
    - 32.0*a*b*d - 8.0*a*c*d*d*e
    + 8.0*a*c*d*d + 16.0*a*d*d*d*d
    + 8.0*a*d*d*e - 8.0*a*d*d
    - 16.0*b*b*b*b*c*c - 8.0*b*b*b*b*c*e
    - 32.0*b*b*b*b*c - 8.0*b*b*b*b*e
    - 16.0*b*b*b*b + 32.0*b*b*b*c*c*d
    + 16.0*b*b*b*c*d*e + 48.0*b*b*b*c*d
    + 16.0*b*b*b*d*e + 16.0*b*b*b*d
    + 8.0*b*b*c*c*c*e - 4.0*b*b*c*c*c
    - 16.0*b*b*c*c*d*d + 4.0*b*b*c*c*e*e
    + 6.0*b*b*c*c*e + 2.0*b*b*c*c
    - 8.0*b*b*c*d*d*e - 8.0*b*b*c*d*d
    - 4.0*b*b*c*e + 16.0*b*b*c
    - 8.0*b*b*d*d*e + 8.0*b*b*d*d
    - 4.0*b*b*e*e - 2.0*b*b*e + 10.0*b*b
    - 8.0*b*c*c*c*d + 4.0*b*c*c*d*e
    - 28.0*b*c*c*d - 8.0*b*c*d*e
    - 16.0*b*c*d - 12.0*b*d*e + 4.0*b*d
    - 2.0*c*c*c*c*e + 2.0*c*c*c*c
    + 4.0*c*c*c*d*d - c*c*c*e*e
    - 2.0*c*c*c*e + 3.0*c*c*c
    + 6.0*c*c*d*d*e + 6.0*c*c*d*d
    - c*c*e*e + 2.0*c*c*e - c*c
    - 8.0*c*d*d*d*d + 4.0*c*d*d*e
    + 8.0*c*d*d + c*e*e + 2.0*c*e
    - 3.0*c - 8.0*d*d*d*d
    - 2.0*d*d*e + 6.0*d*d + e*e - 1.0);


    constexpr long double q = (4.0*a*a*a
    - 4.0*a*b*b
    - 4.0*a*b*d
    - 2.0*a*c
    - a*e
    - a
    + 2.0*b*b*c
    + b*b*e
    + 2.0*b*b
    + 2.0*b*d
    + d*d)
    /
    (16.0*a*a*a*a
    - 32.0*a*a*b*b
    - 32.0*a*a*b*d
    - 4.0*a*a*c*c
    - 16.0*a*a*c
    - 4.0*a*a*e
    - 8.0*a*a
    + 32.0*a*b*b*c
    + 8.0*a*b*b*e
    + 32.0*a*b*b
    + 16.0*a*b*c*d
    + 32.0*a*b*d
    + 8.0*a*d*d
    - 8.0*b*b*c*c
    - 4.0*b*b*c*e
    - 16.0*b*b*c
    - 4.0*b*b*e
    - 8.0*b*b
    - 8.0*b*c*d
    - 8.0*b*d
    + 2.0*c*c*c
    + c*c*e
    + 5.0*c*c
    - 4.0*c*d*d
    + 2.0*c*e
    + 4.0*c
    - 4.0*d*d
    + e
    + 1.0);


    constexpr long double r = (-8.0*a*a*a*a*e
    + 8.0*a*a*a*a
    - 4.0*a*a*b*b*c
    + 16.0*a*a*b*b*e
    - 20.0*a*a*b*b
    + 8.0*a*a*b*c*d
    + 16.0*a*a*b*d*e
    - 8.0*a*a*b*d
    + 4.0*a*a*c*c*e
    - 4.0*a*a*c*c
    - 4.0*a*a*c*d*d
    + 8.0*a*a*c*e
    - 8.0*a*a*c
    - 4.0*a*a*d*d
    + 2.0*a*a*e*e
    - 2.0*a*a
    + 8.0*a*b*b*b*b
    - 8.0*a*b*b*b*d
    - 20.0*a*b*b*c*e
    + 20.0*a*b*b*c
    - 8.0*a*b*b*d*d
    - 4.0*a*b*b*e*e
    - 8.0*a*b*b*e
    + 12.0*a*b*b
    - 12.0*a*b*c*d*e
    + 12.0*a*b*c*d
    + 8.0*a*b*d*d*d
    - 12.0*a*b*d*e
    + 12.0*a*b*d
    - 4.0*a*d*d*e
    + 4.0*a*d*d
    - 4.0*b*b*b*b*c
    - 2.0*b*b*b*b*e
    - 4.0*b*b*b*b
    + 8.0*b*b*b*c*d
    + 4.0*b*b*b*d*e
    + 4.0*b*b*b*d
    + 6.0*b*b*c*c*e
    - 4.0*b*b*c*c
    - 4.0*b*b*c*d*d
    + 3.0*b*b*c*e*e
    + 6.0*b*b*c*e
    - 5.0*b*b*c
    - 2.0*b*b*d*d*e
    + 2.0*b*b*d*d
    + b*b*e*e
    + 2.0*b*b*e
    - b*b
    - 4.0*b*c*c*d
    + 4.0*b*c*d*e
    - 12.0*b*c*d
    - 4.0*b*d
    - 2.0*c*c*c*e
    + 2.0*c*c*c
    + 2.0*c*c*d*d
    - c*c*e*e
    - 2.0*c*c*e
    + 3.0*c*c
    + 4.0*c*d*d*e
    - c*e*e
    + c
    - 2.0*d*d*d*d
    + 2.0*d*d*e)
    /
    (32.0*a*a*a*a*b*b
    - 64.0*a*a*a*a*b*d
    - 16.0*a*a*a*a*c*e
    + 16.0*a*a*a*a*c
    + 32.0*a*a*a*a*d*d
    + 16.0*a*a*a*a*e
    - 16.0*a*a*a*a
    - 64.0*a*a*b*b*b*b
    + 64.0*a*a*b*b*b*d
    - 8.0*a*a*b*b*c*c
    + 32.0*a*a*b*b*c*e
    - 64.0*a*a*b*b*c
    + 64.0*a*a*b*b*d*d
    - 40.0*a*a*b*b*e
    + 16.0*a*a*b*b
    + 16.0*a*a*b*c*c*d
    + 32.0*a*a*b*c*d*e
    + 32.0*a*a*b*c*d
    - 64.0*a*a*b*d*d*d
    - 16.0*a*a*b*d*e
    + 64.0*a*a*b*d
    + 4.0*a*a*c*c*c*e
    - 4.0*a*a*c*c*c
    - 8.0*a*a*c*c*d*d
    + 12.0*a*a*c*c*e
    - 12.0*a*a*c*c
    - 32.0*a*a*c*d*d
    + 4.0*a*a*c*e*e
    - 12.0*a*a*c*e
    + 8.0*a*a*c
    - 8.0*a*a*d*d*e
    - 16.0*a*a*d*d
    - 4.0*a*a*e*e
    - 4.0*a*a*e
    + 8.0*a*a
    + 64.0*a*b*b*b*b*c
    + 16.0*a*b*b*b*b*e
    + 64.0*a*b*b*b*b
    - 96.0*a*b*b*b*c*d
    - 32.0*a*b*b*b*d*e
    - 64.0*a*b*b*b*d
    - 32.0*a*b*b*c*c*e
    + 32.0*a*b*b*c*c
    - 8.0*a*b*b*c*e*e
    + 8.0*a*b*b*c*e
    + 16.0*a*b*b*d*d*e
    - 48.0*a*b*b*d*d
    + 8.0*a*b*b*e*e
    + 24.0*a*b*b*e
    - 32.0*a*b*b
    - 16.0*a*b*c*c*d*e
    + 16.0*a*b*c*c*d
    + 32.0*a*b*c*d*d*d
    - 16.0*a*b*c*d*e
    + 16.0*a*b*c*d
    + 32.0*a*b*d*d*d
    + 32.0*a*b*d*e
    - 32.0*a*b*d
    - 8.0*a*c*d*d*e
    + 8.0*a*c*d*d
    + 16.0*a*d*d*d*d
    + 8.0*a*d*d*e
    - 8.0*a*d*d
    - 16.0*b*b*b*b*c*c
    - 8.0*b*b*b*b*c*e
    - 32.0*b*b*b*b*c
    - 8.0*b*b*b*b*e
    - 16.0*b*b*b*b
    + 32.0*b*b*b*c*c*d
    + 16.0*b*b*b*c*d*e
    + 48.0*b*b*b*c*d
    + 16.0*b*b*b*d*e
    + 16.0*b*b*b*d
    + 8.0*b*b*c*c*c*e
    - 4.0*b*b*c*c*c
    - 16.0*b*b*c*c*d*d
    + 4.0*b*b*c*c*e*e
    + 6.0*b*b*c*c*e
    + 2.0*b*b*c*c
    - 8.0*b*b*c*d*d*e
    - 8.0*b*b*c*d*d
    - 4.0*b*b*c*e
    + 16.0*b*b*c
    - 8.0*b*b*d*d*e
    + 8.0*b*b*d*d
    - 4.0*b*b*e*e
    - 2.0*b*b*e
    + 10.0*b*b
    - 8.0*b*c*c*c*d
    + 4.0*b*c*c*d*e
    - 28.0*b*c*c*d
    - 8.0*b*c*d*e
    - 16.0*b*c*d
    - 12.0*b*d*e
    + 4.0*b*d
    - 2.0*c*c*c*c*e
    + 2.0*c*c*c*c
    + 4.0*c*c*c*d*d
    - c*c*c*e*e
    - 2.0*c*c*c*e
    + 3.0*c*c*c
    + 6.0*c*c*d*d*e
    + 6.0*c*c*d*d
    - c*c*e*e
    + 2.0*c*c*e
    - c*c
    - 8.0*c*d*d*d*d
    + 4.0*c*d*d*e
    + 8.0*c*d*d
    + c*e*e
    + 2.0*c*e
    - 3.0*c
    - 8.0*d*d*d*d
    - 2.0*d*d*e
    + 6.0*d*d
    + e*e
    - 1.0);

    constexpr std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadControlPoints>
        M_inv{{
        {f,g,h,g,h,h,g,h,g},
        {g,i,j,m,j,l,m,l,n},
        {h,j,p,j,q,q,l,r,l},
        {g,m,j,i,l,j,n,l,m},
        {h,j,q,l,p,r,j,q,l},
        {h,l,q,j,r,p,l,q,j},
        {g,m,l,n,j,l,i,j,m},
        {h,l,r,l,q,q,j,p,j},
        {g,n,l,m,l,j,m,j,i}
        }};

    return M_inv;
}

template <FloatType T>
constexpr inline
std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadControlPoints>
GetExtendedQuadMatrix()
{
    constexpr T a = cmc_rbf<T>(std::sqrt(static_cast<T>(8.0)));
    constexpr T b = cmc_rbf<T>(2.0);
    constexpr T c = cmc_rbf<T>(static_cast<T>(4.0));
    constexpr T d = cmc_rbf<T>(std::sqrt(20.0));
    constexpr T e = cmc_rbf<T>(std::sqrt(32.0));

    constexpr std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadControlPoints>
        M{{
        {1.0,a,b,a,b,b,a,b,a}, 
        {a,1.0,b,c,b,d,c,d,e},
        {b,b,1.0,b,a,a,d,c,d},
        {a,c,b,1.0,d,b,e,d,c},
        {b,b,a,d,1.0,c,b,a,d},
        {b,d,a,b,c,1.0,d,a,b},
        {a,c,d,e,b,d,1.0,b,c},
        {b,d,c,d,a,a,b,1.0,b},
        {a,e,d,c,d,b,c,b,1.0}
        }};

    return M;
}

template <FloatType T>
constexpr inline
std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadPredictionPoints>
GetExtendedQuadPredictionCoordsEvaluation()
{
    constexpr std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadPredictionPoints>
        eval_coords{{
        {cmc_rbf<T>(1.0),            cmc_rbf<T>(std::sqrt(13.0)), cmc_rbf<T>(std::sqrt(5.0)), cmc_rbf<T>(std::sqrt(5.0)), cmc_rbf<T>(3.0),             cmc_rbf<T>(1.0),             cmc_rbf<T>(std::sqrt(13.0)), cmc_rbf<T>(std::sqrt(5.0)),  cmc_rbf<T>(std::sqrt(5.0))},
        {cmc_rbf<T>(1.0),            cmc_rbf<T>(std::sqrt(5.0)),  cmc_rbf<T>(1.0),            cmc_rbf<T>(std::sqrt(5.0)), cmc_rbf<T>(std::sqrt(5.0)),  cmc_rbf<T>(std::sqrt(5.0)),  cmc_rbf<T>(std::sqrt(13.0)), cmc_rbf<T>(3.0),             cmc_rbf<T>(std::sqrt(13.0))},
        {cmc_rbf<T>(std::sqrt(2.0)), cmc_rbf<T>(std::sqrt(10.0)), cmc_rbf<T>(std::sqrt(2.0)), cmc_rbf<T>(std::sqrt(2.0)), cmc_rbf<T>(std::sqrt(10.0)), cmc_rbf<T>(std::sqrt(2.0)),  cmc_rbf<T>(std::sqrt(18.0)), cmc_rbf<T>(std::sqrt(10.0)), cmc_rbf<T>(std::sqrt(10.0))}
        }};

    return eval_coords;
}

template <FloatType T>
inline std::array<T, kNumExtendedQuadControlPoints>
EvaluateExtendedQuadWeights(const std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadControlPoints>& M_inv, const std::array<T, kNumExtendedQuadControlPoints>& control_values)
{
    std::array<T, kNumExtendedQuadControlPoints> weights;

    /* Compute matrix vector multiplication in order to obtainthe weights */
    int idx = 0;
    for (auto row_iter = M_inv.begin(); row_iter != M_inv.end(); ++row_iter, ++idx)
    {
        weights[idx] = std::inner_product(row_iter->begin(), row_iter->end(), control_values.begin(), static_cast<T>(0.0));
    }

    return weights;
}

template <FloatType T>
inline std::vector<T>
EvaluateExtendedQuadPrediction(const std::array<T, kNumExtendedQuadControlPoints>& weights, const std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadPredictionPoints>& eval_prediction_coords)
{
    std::vector<T> predictions(kNumExtendedQuadPredictionPoints);

    for (int idx{0}; idx < kNumExtendedQuadPredictionPoints; ++idx)
    {
        const T prediction = std::inner_product(weights.begin(), weights.end(), eval_prediction_coords[idx].begin(), static_cast<T>(0.0));

        predictions[idx] = prediction;        
    }

    return predictions;
}

template <FloatType T>
std::array<T, kNumExtendedQuadPredictionPoints>
PerformExtendedQuadPrediction(const std::array<T, kNumExtendedQuadControlPoints>& control_values)
{
    /* Get the inverse sytem matrix */
    constexpr std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadControlPoints> M_inv = GetExtendedQuadMatrixInverse<T>();

    /* Compute the weights for the quad prediction */
    const std::array<T, kNumExtendedQuadControlPoints> weights{
        M_inv[0][0] * control_values[0] + M_inv[0][1] * control_values[1] + M_inv[0][2] * control_values[2] + M_inv[0][3] * control_values[3] + M_inv[0][4] * control_values[4] + M_inv[0][5] * control_values[5] + M_inv[0][6] * control_values[6] + M_inv[0][7] * control_values[7] + M_inv[0][8] * control_values[8],
        M_inv[1][0] * control_values[0] + M_inv[1][1] * control_values[1] + M_inv[1][2] * control_values[2] + M_inv[1][3] * control_values[3] + M_inv[1][4] * control_values[4] + M_inv[1][5] * control_values[5] + M_inv[1][6] * control_values[6] + M_inv[1][7] * control_values[7] + M_inv[1][8] * control_values[8],
        M_inv[2][0] * control_values[0] + M_inv[2][1] * control_values[1] + M_inv[2][2] * control_values[2] + M_inv[2][3] * control_values[3] + M_inv[2][4] * control_values[4] + M_inv[2][5] * control_values[5] + M_inv[2][6] * control_values[6] + M_inv[2][7] * control_values[7] + M_inv[2][8] * control_values[8],
        M_inv[3][0] * control_values[0] + M_inv[3][1] * control_values[1] + M_inv[3][2] * control_values[2] + M_inv[3][3] * control_values[3] + M_inv[3][4] * control_values[4] + M_inv[3][5] * control_values[5] + M_inv[3][6] * control_values[6] + M_inv[3][7] * control_values[7] + M_inv[3][8] * control_values[8],
        M_inv[4][0] * control_values[0] + M_inv[4][1] * control_values[1] + M_inv[4][2] * control_values[2] + M_inv[4][3] * control_values[3] + M_inv[4][4] * control_values[4] + M_inv[4][5] * control_values[5] + M_inv[4][6] * control_values[6] + M_inv[4][7] * control_values[7] + M_inv[4][8] * control_values[8],
        M_inv[5][0] * control_values[0] + M_inv[5][1] * control_values[1] + M_inv[5][2] * control_values[2] + M_inv[5][3] * control_values[3] + M_inv[5][4] * control_values[4] + M_inv[5][5] * control_values[5] + M_inv[5][6] * control_values[6] + M_inv[5][7] * control_values[7] + M_inv[5][8] * control_values[8],
        M_inv[6][0] * control_values[0] + M_inv[6][1] * control_values[1] + M_inv[6][2] * control_values[2] + M_inv[6][3] * control_values[3] + M_inv[6][4] * control_values[4] + M_inv[6][5] * control_values[5] + M_inv[6][6] * control_values[6] + M_inv[6][7] * control_values[7] + M_inv[6][8] * control_values[8],
        M_inv[7][0] * control_values[0] + M_inv[7][1] * control_values[1] + M_inv[7][2] * control_values[2] + M_inv[7][3] * control_values[3] + M_inv[7][4] * control_values[4] + M_inv[7][5] * control_values[5] + M_inv[7][6] * control_values[6] + M_inv[7][7] * control_values[7] + M_inv[7][8] * control_values[8],
        M_inv[8][0] * control_values[0] + M_inv[8][1] * control_values[1] + M_inv[8][2] * control_values[2] + M_inv[8][3] * control_values[3] + M_inv[8][4] * control_values[4] + M_inv[8][5] * control_values[5] + M_inv[8][6] * control_values[6] + M_inv[8][7] * control_values[7] + M_inv[8][8] * control_values[8]
    };

    /* Get the evaluation coordinates */
    constexpr std::array<std::array<T, kNumExtendedQuadControlPoints>, kNumExtendedQuadPredictionPoints> eval_coords = GetExtendedQuadPredictionCoordsEvaluation<T>();

    /* Perform the prediction */
    const std::array<T, kNumExtendedQuadPredictionPoints> prediction{
        weights[0] * eval_coords[0][0] + weights[1] * eval_coords[0][1] + weights[2] * eval_coords[0][2] + weights[3] * eval_coords[0][3] + weights[4] * eval_coords[0][4] + weights[5] * eval_coords[0][5] + weights[6] * eval_coords[0][6] + weights[7] * eval_coords[0][7] + weights[8] * eval_coords[0][8],
        weights[0] * eval_coords[1][0] + weights[1] * eval_coords[1][1] + weights[2] * eval_coords[1][2] + weights[3] * eval_coords[1][3] + weights[4] * eval_coords[1][4] + weights[5] * eval_coords[1][5] + weights[6] * eval_coords[1][6] + weights[7] * eval_coords[1][7] + weights[8] * eval_coords[1][8],
        weights[0] * eval_coords[2][0] + weights[1] * eval_coords[2][1] + weights[2] * eval_coords[2][2] + weights[3] * eval_coords[2][3] + weights[4] * eval_coords[2][4] + weights[5] * eval_coords[2][5] + weights[6] * eval_coords[2][6] + weights[7] * eval_coords[2][7] + weights[8] * eval_coords[2][8],
    };

    return prediction;
}

/** End Extended Quad Element Operators */

}

#endif /* !CMC_AMR_LOSSY_PAR_MULTI_RES_INTERPOLATION_UTIL_HXX */
