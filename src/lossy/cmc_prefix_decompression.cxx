#include "lossy/cmc_prefix_decompression.hxx"
#include "netcdf/cmc_nc_io_conventions.hxx"
#include "netcdf/cmc_nc_reader.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adapt_callbacks.h"

#include <variant>

namespace cmc
{

namespace prefix
{


[[maybe_unused]]
static bool
IsPrefixExtractionFormat(const std::vector<NcAttribute>& global_attributes)
{
    /* Check for a given attribute indicating the compression which has been used */
    auto comp_scheme_iter = FindAttribute(global_attributes, kCompressionSchemeAttrName);

    /* Check if the attribute has been found */
    if (comp_scheme_iter == global_attributes.end())
    {
        return false;
    }

    /* Get the attribute and check it's value */
    const NcAttribute& comp_scheme_attr = *comp_scheme_iter;

    /* Get the value of the compression scheme */
    const CompressionScheme scheme = static_cast<CompressionScheme>(std::get<CompressionSchemeType>(comp_scheme_attr.GetValue()));

    /* Check if the compression scheme equals the PrefixExtraction method */
    if (scheme == CompressionScheme::PrefixExtraction)
    {
        return true;
    } else
    {
        return false;
    }
}

static size_t
EstimateNumberCompressionVariables(const std::vector<NcVariable>& variable_hulls)
{
    //TODO: implement
    return 2;
}

static void
TrimAtSecondToLastUnderscore(std::string& name)
{
    /* Find the position of the last '_' corresponding to the appended time step*/
    const std::string::size_type last_pos = name.rfind('_');
    cmc_debug_msg("last pos is: ", last_pos);
    /* Find the second to last '_' corresponding to the appended split variable ID */
    const std::string::size_type pos = name.rfind('_', last_pos - 1);
    cmc_debug_msg("pos is: ", pos);
    /* Trim the name */
    name.erase(pos, name.size() - pos);
}

static
std::string
GetVariableBaseName(const NcVariable& variable)
{
    /* Retrieve the name of the variable */
    std::string base_name = variable.GetName();

    /* Trim the name to it's base name */
    TrimAtSecondToLastUnderscore(base_name);

    return base_name;
}

static
std::vector<CompressedVariableInfo>
GetCompressionVariableInfo(const std::vector<NcVariable>& variable_hulls)
{
    std::vector<CompressedVariableInfo> variables;
    variables.reserve(EstimateNumberCompressionVariables(variable_hulls));

    for (auto vh_iter = variable_hulls.begin(); vh_iter != variable_hulls.end(); ++vh_iter)
    {
        /* Get the attributes that belong to the variable */
        const std::vector<NcAttribute>& attributes = vh_iter->GetAttributes();

        /* Get the base_name of the variable without the appended information */
        const std::string basename = GetVariableBaseName(*vh_iter);

        /* Check if a variable info with the given base name already exists */
        auto var_info_iter = std::find_if(variables.begin(), variables.end(), [&basename](const CompressedVariableInfo& info){
            return info.base_name == basename;
        });

        if (var_info_iter != variables.end())
        {
            /* We have found a variable information with the given name.
             * Check whether this variable has a higher time step and split variable ID
             * and if so, update the numbers */

            /* Read in the global_context_information  */
            const auto global_context_information_iter = FindAttribute(attributes, "global_context");
            cmc_assert(global_context_information_iter != attributes.end());
            const int new_spitvar_id = std::get<int>(global_context_information_iter->GetValue()); 
            if (var_info_iter->num_split_vars < new_spitvar_id)
            {
                /* In case the new ID is larger than the previous IDs, we update it */
                var_info_iter->num_split_vars = new_spitvar_id;
            }
            
            /* Read in the corresponding time step */
            const auto time_step_iter = FindAttribute(attributes, "time_step");
            cmc_assert(time_step_iter != attributes.end());
            const int new_time_step = std::get<int>(time_step_iter->GetValue());
            if (var_info_iter->num_time_steps < new_time_step)
            {
                /* In case the new time step is larger than the previous time steps, we update it */
                var_info_iter->num_time_steps = new_time_step;
            }

        } else
        {
            /* We have not found a variable information with the given name and need to create one */
            variables.emplace_back(basename);
            
            /* Read in the id */
            const auto id_iter = FindAttribute(attributes, "id");
            cmc_assert(id_iter != attributes.end());
            variables.back().id = std::get<int>(id_iter->GetValue()); 

            /* Read in the initial refinement level  */
            const auto initial_ref_lvl_iter = FindAttribute(attributes, "initial_refinement_level");
            cmc_assert(initial_ref_lvl_iter != attributes.end());
            variables.back().intial_refinement_level = std::get<int>(initial_ref_lvl_iter->GetValue()); 

            /* Read in the initial data layout */
            const auto initial_layout_iter = FindAttribute(attributes, "initial_layout");
            cmc_assert(initial_layout_iter != attributes.end());
            variables.back().initial_layout = static_cast<DataLayout>(std::get<int>(initial_layout_iter->GetValue())); 

            /* Read in the pre-compression data layout */
            const auto pre_compression_layout_iter = FindAttribute(attributes, "pre_compression_layout");
            cmc_assert(pre_compression_layout_iter != attributes.end());
            variables.back().pre_compression_layout = static_cast<DataLayout>(std::get<int>(pre_compression_layout_iter->GetValue())); 

            /* Read in the global_context_information  */
            const auto global_context_information_iter = FindAttribute(attributes, "global_context");
            cmc_assert(global_context_information_iter != attributes.end());
            variables.back().global_context_information = std::get<int>(global_context_information_iter->GetValue()); 

            /* Read in the data type of the variable */
            const auto data_type_iter = FindAttribute(attributes, "data_type");
            cmc_assert(data_type_iter != attributes.end());
            variables.back().type = static_cast<CmcType>(std::get<int>(data_type_iter->GetValue()));

            /* Read in the missing value */
            const auto missing_value_iter = FindAttribute(attributes, "missing_value");
            cmc_assert(missing_value_iter != attributes.end());
            variables.back().missing_value = missing_value_iter->GetValue(); 

            GeoDomain reconstructed_domain;

            /* Read lon, lat, lev data and create a geo domain of it */
            const auto lon_iter = FindAttribute(attributes, "lon");
            if (lon_iter != attributes.end())
            {
                /* Set the found lon diemnsion length in the geo domain */
                reconstructed_domain.UpdateDimension(DimensionInterval(Dimension::Lon, 0, std::get<int>(lon_iter->GetValue()))); 
            }

            const auto lat_iter = FindAttribute(attributes, "lat");
            if (lat_iter != attributes.end())
            {
                /* Set the found lon diemnsion length in the geo domain */
                reconstructed_domain.UpdateDimension(DimensionInterval(Dimension::Lat, 0, std::get<int>(lat_iter->GetValue()))); 
            }

            const auto lev_iter = FindAttribute(attributes, "lon");
            if (lev_iter != attributes.end())
            {
                /* Set the found lon diemnsion length in the geo domain */
                reconstructed_domain.UpdateDimension(DimensionInterval(Dimension::Lev, 0, std::get<int>(lev_iter->GetValue()))); 
            }

            variables.back().domain = reconstructed_domain;

            /* Read in the corresponding time step */
            const auto time_step_iter = FindAttribute(attributes, "time_step");
            cmc_assert(time_step_iter != attributes.end());
            variables.back().num_time_steps = std::get<int>(time_step_iter->GetValue());
        }
    }

    /* Since, we have collected the the highest time step and the highest split var ID, we need to add +1 to it, because the time steps and IDs start at zero 
     * and therefore, do not refer to the actual number of time steps and split variables */
    for (auto var_iter = variables.begin(); var_iter != variables.end(); ++var_iter)
    {
        var_iter->num_split_vars += 1;
        var_iter->num_time_steps += 1;
    }

    return variables;
}

static
std::vector<CompressedVariableInfo>
EvaluateCompressionVariables(const std::vector<NcVariable>& variable_hulls)
{
    /* Get Information about the variables that have been compressed */
    std::vector<CompressedVariableInfo> variables = GetCompressionVariableInfo(variable_hulls);
    return variables;
}

static
std::vector<ByteVar>
SetUpByteVariableForDecompression(const std::vector<CompressedVariableInfo>& var_info)
{
    std::vector<ByteVar> byte_variables;
    byte_variables.reserve(var_info.size());

    for (auto var_iter = var_info.begin(); var_iter != var_info.end(); ++var_iter)
    {
        /* Create a new byte variable from the given information */
        byte_variables.emplace_back(var_iter->id, var_iter->type, var_iter->base_name, var_iter->domain, var_iter->initial_layout, var_iter->pre_compression_layout,
                                    var_iter->global_context_information, var_iter->missing_value);
        
        /* Create the base mesh for this variable */
        t8_forest_t base_mesh = ReconstructBaseMesh(GetDimensionalityOfDataLayout(var_iter->initial_layout), MPI_COMM_WORLD);

        /* Assign the base mesh to the byte variable */
        byte_variables.back().SetMesh(base_mesh);

        /* The base mesh consists of a single element in either way. In order to append the extracted prefixes,
         * we need to add an initial empty CompressionValue "representing the root element" */
        byte_variables.back().AssignCompressionValueForDecompressionStart();
    }

    return byte_variables;
}


void
Decompressor::Setup()
{
    /* Create a reader for the compressed file */
    NcReader reader(file_name_);

    /* Get the global attributes */
    std::vector<NcAttribute> global_atts = reader.ReadGlobalAttrtibutes();

    /* Ensure that the file obeys to the compression norm */
    if (IsPrefixExtractionFormat(global_atts) == false)
    {
        cmc_err_msg("The file does not hold prefix-comrpessed data. Therefore, a decompression cannot be appllied.");
    }

    /* Get the meta data of the variables in order to evaluate which variable corresponds to which feature */
    std::vector<NcVariable> variable_hulls = reader.ReadVariableMetaData();

    /* Check all variables and in particular retrieve the time steps and number of time steps encoded */
    variable_info_ = EvaluateCompressionVariables(variable_hulls);

    /* Set up the byte variables for the first time step */
    compression_variables_ = SetUpByteVariableForDecompression(variable_info_);
}


void
Decompressor::Decompress()
{

}

static
bool
IsThereACompressedVariableWithTheGivenName(const std::vector<CompressedVariableInfo>& var_infos, const std::string& variable_base_name)
{
    /* Try to find a variable with the given name */
    auto var_info_iter = std::find_if(var_infos.begin(), var_infos.end(), [&variable_base_name](const CompressedVariableInfo& var_info){
        return !var_info.base_name.compare(variable_base_name);
    });

    /* Exit, if the variable does not exist */
    if (var_info_iter == var_infos.end())
    {
        return false;
    }

    return true;
}

static
int
GetNumberOfTimeSteps(const std::vector<CompressedVariableInfo>& var_infos, const std::string& variable_base_name)
{
    /* Try to find a variable with the given name */
    auto var_info_iter = std::find_if(var_infos.begin(), var_infos.end(), [&variable_base_name](const CompressedVariableInfo& var_info){
        return !var_info.base_name.compare(variable_base_name);
    });

    cmc_assert(var_info_iter != var_infos.end());

    /* Return the number of time steps which are stored */
    if (var_info_iter != var_infos.end())
    {
        return var_info_iter->num_time_steps;
    } else
    {
        return 0;
    }
}

static
int
GetNumberOfSplitVariables(const std::vector<CompressedVariableInfo>& var_infos, const std::string& variable_base_name)
{
    /* Try to find a variable with the given name */
    auto var_info_iter = std::find_if(var_infos.begin(), var_infos.end(), [&variable_base_name](const CompressedVariableInfo& var_info){
        return !var_info.base_name.compare(variable_base_name);
    });

    cmc_assert(var_info_iter != var_infos.end());

    /* Return the number of time steps which are stored */
    if (var_info_iter != var_infos.end())
    {
        return var_info_iter->num_split_vars;
    } else
    {
        return 0;
    }
}

static
int 
GetInitialRefinementLevel(const std::vector<CompressedVariableInfo>& var_infos, const std::string& variable_base_name)
{
    /* Try to find a variable with the given name */
    auto var_info_iter = std::find_if(var_infos.begin(), var_infos.end(), [&variable_base_name](const CompressedVariableInfo& var_info){
        return !var_info.base_name.compare(variable_base_name);
    });

    cmc_assert(var_info_iter != var_infos.end());

    /* Return the number of time steps which are stored */
    if (var_info_iter != var_infos.end())
    {
        return var_info_iter->intial_refinement_level;
    } else
    {
        return 0;
    }
}

void
Decompressor::DecompressVariable(const std::string& variable_base_name)
{
    /* Check if a variable with the given name exists */
    if (!IsThereACompressedVariableWithTheGivenName(variable_info_, variable_base_name))
    {
        cmc_err_msg("A variable with the (base)name ", variable_base_name, " does not exist in the compressed file.");
    }

    /* Get the number of time steps which are encoded */
    const int num_time_steps = GetNumberOfTimeSteps(variable_info_, variable_base_name);

    /* Get the number of split variables */
    const int num_split_variables = GetNumberOfSplitVariables(variable_info_, variable_base_name);

    /* Get an iterator to the corresponding byte variable */
    auto var_iter = std::find_if(compression_variables_.begin(), compression_variables_.end(), [&variable_base_name](const ByteVar& byte_var){
        return !byte_var.GetName().compare(variable_base_name);
    });

    cmc_assert(var_iter != compression_variables_.end());

    /* Create a reader for the compressed file */
    NcReader reader(file_name_);

    /* Decompress the variable data by iterating over all time steps and all split variables */
    for (int time_step = 0; time_step < num_time_steps; ++time_step)
    {
        for (int split_id = 0; split_id < num_split_variables; ++split_id)
        {
            /* Create the name of the variable we want to decompress */
            const std::string nc_var_name = variable_base_name + "_" + std::to_string(split_id) + "_" + std::to_string(time_step);

            /* Read the (serialized) data of the variable from the file */
            std::vector<uint8_t> serialized_data = reader.ReadVariableData<uint8_t>(nc_var_name);
            cmc_debug_msg("Size of serialized data:", serialized_data.size());
            /* Get the initial refinement level of the data */
            const int initial_refinement_level = GetInitialRefinementLevel(variable_info_, variable_base_name);

            /* Create an adapt data class which handles the compressed byte stream */
            DecompressPrefixAdaptData adapt_data(*var_iter, serialized_data);

            /* Decode the levels of the prefix encodings */
            for (int lvl = 0; lvl < initial_refinement_level; ++lvl)
            {
                /* Initialize the decompression iteration by sufficiently allocating memory */
                adapt_data.InitializeDecompressionIteration();

                /* Decompress the next level */
                t8_forest_t decompressed_forest = t8_forest_new_adapt(adapt_data.GetCurrentMesh(), DecompressPrefixEncoding, 0, 0, static_cast<void*>(&adapt_data));
                
                /* Store the decomrpessed mesh */
                adapt_data.SetCurrentMesh(decompressed_forest);
                
                /* Finalize the decompression iteration */
                adapt_data.FinalizeDecompressionIteration();
            }

            /* Perform the suffix decoding */
            adapt_data.InitializeSuffixDecompression();

            /* Decompress the next level */
            t8_forest_t decompressed_forest = t8_forest_new_adapt(adapt_data.GetCurrentMesh(), DecompressSuffixEncoding, 0, 0, static_cast<void*>(&adapt_data));
            
            /* Store the decomrpessed mesh */
            adapt_data.SetCurrentMesh(decompressed_forest);

            /* Finalize the decompression iteration */
            adapt_data.FinalizeSuffixDecompression();

            //....

        }
    }
}


}

}
