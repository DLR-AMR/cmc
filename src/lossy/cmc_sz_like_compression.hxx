#ifndef CMC_SZ_LIKE_COMPRESSION_HXX
#define CMC_SZ_LIKE_COMPRESSION_HXX

#include "utilities/cmc_utilities.hxx"
#include "t8code/cmc_t8_data_variables.hxx"

#include <vector>
#include <variant>

namespace cmc
{

//constexpr uint8_t kNoEncoding{0b00000000};//  = 0x00;//{0b00000000};
//constexpr uint8_t kConstFit{0b00000001};// = 0x01;//{0b00000001};
//constexpr uint8_t kLinFit{0b00000010};// = 0x02; //{0b00000010};
//constexpr uint8_t kQuadFit{0b00000011};// = 0x03; //{0b00000011};

constexpr uint8_t kNoEncoding = 0;
constexpr uint8_t kConstFit = 1;
constexpr uint8_t kLinFit = 2;
constexpr uint8_t kQuadFit = 3;

enum CurveFitting {Constant = 0, Linear, Quadratic};

template<class T>
class SZCurveFittingCompressor;

using CmcSZCurveFittingCompressor = std::variant<SZCurveFittingCompressor<int8_t>, SZCurveFittingCompressor<char>, SZCurveFittingCompressor<int16_t>,
                                                 SZCurveFittingCompressor<int32_t>, SZCurveFittingCompressor<float>, SZCurveFittingCompressor<double>,
                                                 SZCurveFittingCompressor<uint8_t>, SZCurveFittingCompressor<uint16_t>, SZCurveFittingCompressor<uint32_t>,
                                                 SZCurveFittingCompressor<int64_t>, SZCurveFittingCompressor<uint64_t>>;

template<class T>
class SZCurveFittingCompressor
{
public:
    SZCurveFittingCompressor() = delete;
    SZCurveFittingCompressor(const Variable<T>& comrpession_variable_ref)
    : compression_variable_{comrpession_variable_ref} {};

    void Compress();
private:
    const Variable<T>& compression_variable_;
    
    std::vector<uint8_t> encoded_curve_fitting_;
    std::vector<T> not_encoded_values_;
};

class SZCompressor
{
public:
    SZCompressor() = delete;
    SZCompressor(const std::vector<Var>& variables)
    : variables_{variables} {
        compression_data_.reserve(variables.size());
    };

    void Compress()
    {
        /* Setup the comrpession for each variable */
        for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
        {
            std::visit([&](const auto& var){
                compression_data_.emplace_back(var);
            }, var_iter->GetInternalVariable());
        }

        /* Compress each variable */
        for (auto sz_iter = compression_data_.begin(); sz_iter != compression_data_.end(); ++sz_iter)
        {
            std::visit([&](auto&& sz_cv_compressor){
                //CurveFittingCompression Compress;
                //Compress(sz_cv_compressor);
                sz_cv_compressor.Compress();
            }, *sz_iter);
        }

        is_compression_applied_ = true;
    }

private:
    const std::vector<Var>& variables_;

    std::vector<CmcSZCurveFittingCompressor> compression_data_;

    bool is_compression_applied_{false};
};



template<typename T>
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
SZCurveFittingCompressor<T>::Compress()
{
    cmc_assert(compression_variable_.size() > 3);

    cmc_debug_msg("kNoEncoding = ", unsigned(kNoEncoding), ", kConstFit = ", unsigned(kConstFit), ", kLinFit = ", unsigned(kLinFit), ", kQuadFit = ", unsigned(kQuadFit));

    const int num_data_points = static_cast<int>(compression_variable_.size());

    cmc_debug_msg("Num data points: ", num_data_points);

    //std::vector<uint8_t> encoded_curve_fitting_;
    encoded_curve_fitting_.clear();
    encoded_curve_fitting_.reserve(num_data_points / 2 + 1);

    std::vector<T> not_encoded_values_;
    not_encoded_values_.reserve(num_data_points / 2 + 1);

    std::vector<double> ne_vals_remaining_error;
    ne_vals_remaining_error.reserve(num_data_points / 2 + 1);

    cmc_assert(std::holds_alternative<T>(compression_variable_.GetMissingValue()));

    const T missing_value = std::get<T>(compression_variable_.GetMissingValue());
    bool was_last_value_missing_value{false};

    int current_byte_pos = 0;
    int current_byte_stream_index = 0;

    /* The offset for each new encoded flag */
    const int increment_cf_bits = 2;

    encoded_curve_fitting_.push_back(0);

    /* The first three values cannot be encoded */
    encoded_curve_fitting_[current_byte_stream_index] |= kNoEncoding;
    current_byte_pos += increment_cf_bits;
    encoded_curve_fitting_[current_byte_stream_index] |= (kNoEncoding << current_byte_pos);
    current_byte_pos += increment_cf_bits;
    encoded_curve_fitting_[current_byte_stream_index] |= (kNoEncoding << current_byte_pos);

    /* Get the first (not encoded) values of the variable */
    std::vector<T> last_decompressed_values{compression_variable_[0], compression_variable_[1], compression_variable_[2]};

    /* Save the not encoded values sparately */
    not_encoded_values_.push_back(compression_variable_[2]);
    not_encoded_values_.push_back(compression_variable_[1]);
    not_encoded_values_.push_back(compression_variable_[0]);

    //TODO:
    ne_vals_remaining_error.push_back(compression_variable_.GetRemainingMaxAllowedAbsoluteError(2));
    //ne_vals_remaining_error.push_back(...);
    //ne_vals_remaining_error.push_back(...);

    cmc_debug_msg("Start of CF");
    /* Try to encode all remaining data points */
    for (int index = 3; index < num_data_points; ++index)
    {
        //cmc_debug_msg("index ", index);
        /* Chck if we are able to emplace the next encoding in the current byte */
        if (current_byte_pos + increment_cf_bits >= CHAR_BIT)
        {
            //cmc_debug_msg("before new byte is set: ", unsigned(encoded_curve_fitting_[current_byte_stream_index]));
            current_byte_pos = 0;
            ++current_byte_stream_index;
            encoded_curve_fitting_.push_back(0);
        } else
        {
            current_byte_pos += increment_cf_bits;
        }

        /* Calculate all curve fittings */
        //auto [c_cf, l_cf, q_cf] = CalculateCurveFittings(last_decompressed_values);

        T min_val;
        double min_deviation;
        uint8_t bit_encoding;
        bool is_error_satisfied{false};

        /* The value which trying to be fitted */
        const T current_value = compression_variable_[index];

        /* Check for missing values */
        if (ApproxCompare(current_value, missing_value))
        {
            /* Check if it is a series of missing values */
            if (was_last_value_missing_value)
            {
                /* Push back a constant fitting and continue */
                min_deviation = 0.0;
                bit_encoding = kConstFit;
                is_error_satisfied = true;
            } else
            {
                /* Save the missing value without encoding */
                min_val = missing_value;
                is_error_satisfied = false;
                was_last_value_missing_value = true;
            }
        } else
        {
            if (was_last_value_missing_value)
            {
                was_last_value_missing_value = false;
            }

            /* Check if deviation is compliant */
            //const auto [c_cf_error_satisfied, c_cf_abs_deviation] = compression_variable_.CheckErrorBoundsForValue(c_cf, index);
            //const auto [l_cf_error_satisfied, l_cf_abs_deviation] = compression_variable_.CheckErrorBoundsForValue(l_cf, index);
            //const auto [q_cf_error_satisfied, q_cf_abs_deviation] = compression_variable_.CheckErrorBoundsForValue(q_cf, index);


            /* Calculate all curve fittings */
            const std::vector<T> curve_fittings = CalculateCurveFittings(last_decompressed_values);

            /* Calculate all deviaitons and check whether they fullfill all given thresholds */
            const std::vector<ErrorCompliance> error_evaluations = compression_variable_.CheckErrorBoundsForValues(curve_fittings, index);

            /* Get the fitting which introduces the minimal deviation */
            if (error_evaluations[CurveFitting::Constant].max_introduced_error <=  error_evaluations[CurveFitting::Linear].max_introduced_error &&
                error_evaluations[CurveFitting::Constant].max_introduced_error <= error_evaluations[CurveFitting::Quadratic].max_introduced_error)
            {
                /* Const is minimum */
                min_val = curve_fittings[CurveFitting::Constant];
                bit_encoding = kConstFit;
                min_deviation = error_evaluations[CurveFitting::Constant].max_introduced_error;
                is_error_satisfied = error_evaluations[CurveFitting::Constant].is_error_threshold_satisfied;
            } else if (error_evaluations[CurveFitting::Linear].max_introduced_error <=  error_evaluations[CurveFitting::Constant].max_introduced_error &&
                       error_evaluations[CurveFitting::Linear].max_introduced_error <= error_evaluations[CurveFitting::Quadratic].max_introduced_error)
            {
                /* Lin is minimum */
                min_val = curve_fittings[CurveFitting::Linear];
                bit_encoding = kLinFit;
                min_deviation = error_evaluations[CurveFitting::Linear].max_introduced_error;
                is_error_satisfied = error_evaluations[CurveFitting::Linear].is_error_threshold_satisfied;
            } else
            {
                /* Quad is minimum */
                min_val = curve_fittings[CurveFitting::Quadratic];
                bit_encoding = kQuadFit;
                min_deviation = error_evaluations[CurveFitting::Quadratic].max_introduced_error;
                is_error_satisfied = error_evaluations[CurveFitting::Quadratic].is_error_threshold_satisfied;
            }
        }

        /* Check if the error has been fullfilled */
        if (is_error_satisfied)
        {
            //cmc_debug_msg("is any error satisfied");
            /* Push the encoding */
            encoded_curve_fitting_[current_byte_stream_index] |= (bit_encoding << current_byte_pos);
            //cmc_debug_msg("as any error satisfied: ", unsigned(encoded_curve_fitting_[current_byte_stream_index]));
        } else
        {
            /* The encodings did not fullfill the error thresholds */
            encoded_curve_fitting_[current_byte_stream_index] |= (kNoEncoding << current_byte_pos);
            not_encoded_values_.push_back(current_value);
            //cmc_debug_msg("if no encoding: ", unsigned(encoded_curve_fitting_[current_byte_stream_index]));
            //TODO: implement with remaining deviations for further bit analysis/compression
            //not_encoded_vals_allowed_err.push_back(...);
        }

        /* Change last decompressed values */
        last_decompressed_values[2] = last_decompressed_values[1];
        last_decompressed_values[1] = last_decompressed_values[0];
        last_decompressed_values[0] = min_val;
    }

    for (int i = 0; i <  encoded_curve_fitting_.size(); ++i)
    {
        std::cout << unsigned(encoded_curve_fitting_[i]) << ", ";
    }
    cmc_debug_msg("Numelemnts: ", compression_variable_.size(), " Bytes for encoding: ", encoded_curve_fitting_.size(), " no encoded vals: ", not_encoded_values_.size());
    
}



}

#endif /* !CMC_SZ_LIKE_COMPRESSION_HXX */
