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

namespace cmc::par::lossy::rbf::util
{

template<typename T>
concept FloatType = (std::is_floating_point_v<T> && std::is_arithmetic_v<T>);

constexpr int kNumHexControlPoints = 7;
constexpr int kNumHexPredictionPoints = 7;

constexpr int kNumQuadControlPoints = 5;
constexpr int kNumQuadPredictionPoints = 3;

constexpr float eps_general = 0.07;
constexpr float eps_gauss = 0.07;
//constexpr float eps_gauss = 1.0;
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
    //static_assert(Eps<Float>() >= 0.02); //Advised is eps >= 0.1
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
    return cmc_rbf_gaussian_eval(dist);
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




}

#endif /* !CMC_AMR_LOSSY_PAR_MULTI_RES_INTERPOLATION_UTIL_HXX */
