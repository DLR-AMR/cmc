#include "lossy/cmc_sz_like_compression.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <tuple>
#include <vector>

namespace cmc
{
#if 0
constexpr uint8_t kNoEncoding = 0x00;//{0b00000000};
constexpr uint8_t kConstFit = 0x01;//{0b00000001};
constexpr uint8_t kLinFit = 0x02; //{0b00000010};
constexpr uint8_t kQuadFit = 0x03; //{0b00000011};

enum CurveFitting {Constant = 0, Linear, Quadratic};


//template<typename T>
//static
//std::tuple<T, T, T>
//CalculateCurveFittings(const std::vector<T>& last_decompressed_values)
//{
//    const T cf_const = last_decompressed_values[0];
//
//    const T cf_lin = 2 * last_decompressed_values[0] - last_decompressed_values[1];
//
//    const T cf_quad = 3 * last_decompressed_values[0] - 3 * last_decompressed_values[1] + last_decompressed_values[2];
//
//    return std::make_tuple(cf_const, cf_lin, cf_quad);
//}

template<typename T>
static
std::vector<T>
CalculateCurveFittings(const std::vector<T>& last_decompressed_values)
{
    return std::vector<T>{
        last_decompressed_values[0], //Constant Fitting
        2 * last_decompressed_values[0] - last_decompressed_values[1], //Linear Fitting
        3 * last_decompressed_values[0] - 3 * last_decompressed_values[1] + last_decompressed_values[2] //Quadratic Fitting
    };
}

template<class T>
void
SZCurveFittingCompressor::Compress()
{
    cmc_assert(compression_variable_.size() > 3);

    cmc_debug_msg("kNoEncoding = ", kNoEncoding, ", kConstFit = ", kConstFit, ", kLinFit = ", kLinFit, ", kQuadFit = ", kQuadFit);

    const int num_data_points = static_cast<int>(compression_variable_.size());

    cmc_debug_msg("Num data points: ", num_data_points);

    std::vector<uint8_t> encoded_byte_stream;
    encoded_curve_fitting_.reserve( num_data_points / 2 + 1);

    std::vector<T> not_encoded_values_;
    not_encoded_values_.reserve(num_data_points / 2 + 1);

    std::vector<double> ne_vals_remaining_error;
    ne_vals_remaining_error.reserve(num_data_points / 2 + 1);

    int current_byte_pos = 0;
    int current_byte_stream_index = 0;

    /* The offset for each new encoded flag */
    const int increment_cf_bits = 2;

    /* The first three values cannot be encoded */
    encoded_curve_fitting_[current_byte_stream_index] |= NoEncoding;
    current_byte_pos += increment_cf_bits;
    encoded_curve_fitting_[current_byte_stream_index] |= (NoEncoding << current_byte_pos);
    current_byte_pos += increment_cf_bits;
    encoded_curve_fitting_[current_byte_stream_index] |= (NoEncoding << current_byte_pos);

    /* Get the first (not encoded) values of the variable */
    std::vector<T> last_decompressed_values{compression_variable_[0], compression_variable_[1], compression_variable_[2]};

    /* Save the not encoded values sparately */
    not_encoded_values_.push_back(compression_variable_[0]);
    not_encoded_values_.push_back(compression_variable_[1]);
    not_encoded_values_.push_back(compression_variable_[2]);

    //TODO:
    //ne_vals_remaining_error.push_back(...);
    //ne_vals_remaining_error.push_back(...);
    //ne_vals_remaining_error.push_back(...);

    cmc_debug_msg("Start of CF");
    /* Try to encode all remaining data points */
    for (int index = 3; index < num_data_points; ++index)
    {
        /* Chck if we are able to emplace the next encoding in the current byte */
        if (current_byte_pos + increment_cf_bits >= CHAR_BIT)
        {
            current_byte_pos = 0;
            ++current_byte_stream_index;
            encoded_curve_fitting_.emplace_back(0);
        } else
        {
            current_byte_pos += increment_cf_bits;
        }

        /* Calculate all curve fittings */
        //auto [c_cf, l_cf, q_cf] = CalculateCurveFittings(last_decompressed_values);

        /* The value which trying to be fitted */
        //const T current_value = compression_variable_[index];

        /* Check if deviation is compliant */
        //const auto [c_cf_error_satisfied, c_cf_abs_deviation] = compression_variable_.CheckErrorBoundsForValue(c_cf, index);
        //const auto [l_cf_error_satisfied, l_cf_abs_deviation] = compression_variable_.CheckErrorBoundsForValue(l_cf, index);
        //const auto [q_cf_error_satisfied, q_cf_abs_deviation] = compression_variable_.CheckErrorBoundsForValue(q_cf, index);

        T min_val;
        double min_deviation;
        uint8_t bit_encoding;
        bool is_error_satisfied{false};

        /* Calculate all curve fittings */
        const std::vector<T> curve_fittings = CalculateCurveFittings(last_decompressed_values);

        /* Calculate all deviaitons and check whether they fullfill all given thresholds */
        const std::vector<ErrorCompliance> error_evaluations = compression_variable_.CheckErrorBoundsForValues(curve_fittings, index);

        /* Get the fitting which introduces the minimal deviation */
        if (error_evaluations[CurveFitting::Constant].max_introduced_error <=  error_evaluations[CurveFitting::Linear].max_introduced_error &&
            error_evaluations[CurveFitting::Constant].max_introduced_error <= error_evaluations[CurveFitting::Quadratic].max_introduced_error)
        {
            /* Const is minimum */
            min_val = curve_fittings[CurveFittings::Constant];
            bit_encoding = ConstFit;
            min_deviation = error_evaluations[CurveFitting::Constant].max_introduced_error;
            is_error_satisfied = error_evaluations[CurveFitting::Constant].is_error_threshold_satisfied;
        } else if (error_evaluations[CurveFitting::Linear].max_introduced_error <=  error_evaluations[CurveFitting::Constant].max_introduced_error &&
                   error_evaluations[CurveFitting::Linear].max_introduced_error <= error_evaluations[CurveFitting::Quadratic].max_introduced_error)
        {
            /* Lin is minimum */
            min_val = curve_fittings[CurveFittings::Linear];
            bit_encoding = LinFit;
            min_deviation = error_evaluations[CurveFitting::Linear].max_introduced_error;
            is_error_satisfied = error_evaluations[CurveFitting::Linear].is_error_threshold_satisfied;
        } else
        {
            /* Quad is minimum */
            min_val = curve_fittings[CurveFittings::Quadratic];
            bit_encoding = QuadFit;
            min_deviation = error_evaluations[CurveFitting::Quadratic].max_introduced_error;
            is_error_satisfied = error_evaluations[CurveFitting::Quadratic].is_error_threshold_satisfied;
        }

        /* Check if the error has been fulfilled */
        if (is_error_satisfied)
        {
            /* Push the encoding */
            encoded_curve_fitting_[current_byte_stream_index] |= (bit_encoding << current_byte_pos);
        } else
        {
            /* The encodings did not fullfill the error thresholds */
            encoded_curve_fitting_[current_byte_stream_index] |= (NoEncoding << current_byte_pos);
            not_encoded_values_.push_back(current_value);
            //TODO: implement with remaining deviations for further bit analysis/compression
            //not_encoded_vals_allowed_err.push_back(...);
        }
    }
}


#endif
}
