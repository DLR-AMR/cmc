#include "t8code/cmc_t8_data.hxx"
#include "t8code/cmc_t8_mpi.hxx"
#include "t8code/cmc_t8_adapt_callbacks.h"
#include "utilities/cmc_variable_transformer.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "utilities/cmc_geo_utilities.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#endif

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h>
#include <t8_element_c_interface.h>
#include <t8_forest/t8_forest_iterate.h>
#include <t8_forest/t8_forest_vtk.h>
#include <t8_forest/t8_forest_partition.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.hxx>
#endif

#if defined(CMC_WITH_T8CODE) && defined(CMC_ENABLE_DEBUG)
#include <p4est.h>
#include <p8est.h>
#endif

#include <algorithm>
#include <numeric>
#include <iterator>
#include <string>

//TODO: remove
#ifndef CMC_WITH_NETCDF
#define CMC_WITH_NETCDF
#endif
#include <iostream>
#include <bitset>

namespace cmc 
{

void
AmrData::SplitVariables()
{
    if (!compression_settings_.AreThereVariablesToSplit())
    {
        return;
    }

    std::vector<InputVar> splitted_variables;

    const auto CapturekSplitAllVariables = kSplitAllVariables;

    if (auto split_all_vars_iter = std::find_if(compression_settings_.GetSplitVariablesBegin(), compression_settings_.GetSplitVariablesEnd(), [=](const SplitVariable& arg){return arg.GetVariableID() == CapturekSplitAllVariables;});
       split_all_vars_iter != compression_settings_.GetSplitVariablesEnd())
    {
        /* If all variables are ought to be split similarly */
        auto input_variable_iter = input_variables_.begin();

        std::vector<InputVar> single_split_var = SplitIntoSubVariables(*input_variable_iter, split_all_vars_iter->GetSplitDimension());

        /* Reserve memory for all splitted variables */
        splitted_variables.reserve(single_split_var.size() * input_variables_.size());

        /* Add the splitted variables to the new vector */
        splitted_variables.insert(splitted_variables.end(), std::make_move_iterator(single_split_var.begin()), std::make_move_iterator(single_split_var.end()));
        single_split_var.clear();

        for (auto iv_iter = input_variables_.begin() + 1; iv_iter != input_variables_.end(); ++iv_iter)
        {
            single_split_var = SplitIntoSubVariables(*iv_iter, split_all_vars_iter->GetSplitDimension());
            /* Add the splitted variables to the new vector */
            splitted_variables.insert(splitted_variables.end(), std::make_move_iterator(single_split_var.begin()), std::make_move_iterator(single_split_var.end()));
            single_split_var.clear();
        }
    } else
    {
        /* If only some variables will be split */
        for (auto input_variable_iter = input_variables_.begin(); input_variable_iter != input_variables_.end(); ++input_variable_iter)
        {
            auto split_var_iter = std::find_if(compression_settings_.GetSplitVariablesBegin(), compression_settings_.GetSplitVariablesEnd(), [=](const SplitVariable& arg){return arg.GetVariableID() == input_variable_iter->GetID();});
            
            if (split_var_iter != compression_settings_.GetSplitVariablesEnd())
            {
                /* In this case the variable will be split */
                std::vector<InputVar> single_split_var = SplitIntoSubVariables(*input_variable_iter, split_var_iter->GetSplitDimension());
                /* Add the splitted variables to the new vector */
                splitted_variables.insert(splitted_variables.end(), std::make_move_iterator(single_split_var.begin()), std::make_move_iterator(single_split_var.end()));
                single_split_var.clear();
            } else
            {
                /* In this case the variable will be not split and can be moved to the new vector */
                splitted_variables.emplace_back(*std::make_move_iterator(input_variable_iter));
            }
        }
    }

    input_variables_ = std::move(splitted_variables);
}


void
AmrData::BuildInitialMesh()
{
    cmc_assert(!input_variables_.empty());

    /* Since all varibales are defined on the same global domain, we are able to take the domain of the first variable */
    const GeoDomain& global_domain = input_variables_.front().GetGlobalDomain();

    #if 0
    const int dimensionality = DetermineDimensionalityOfTheData(global_domain);

    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(DimensionToElementClass(dimensionality), comm_, 0, 0, 1);

    const int initial_refinement_level = CalculateInitialRefinementLevel(global_domain);

    initial_mesh_.SetInitialRefinementLevel(initial_refinement_level);

    ValidateInitialRefinementLevelForDimensionality(initial_refinement_level, dimensionality);

    /* Construct a forest from the cmesh */
    t8_forest_t initial_forest;
    t8_forest_init(&initial_forest);
    t8_forest_set_cmesh(initial_forest, cmesh, comm_);
    t8_forest_set_scheme(initial_forest, t8_scheme_new_default_cxx());
    t8_forest_set_level(initial_forest, 0);
    t8_forest_commit(initial_forest);
    #endif

    /* The mesh layout is equal to the layout of the domains (the exact ordering of the dimensions is allowed to vary) */
    const DataLayout initial_mesh_layout = input_variables_.front().GetInitialDataLayout();
    #if 0
    AdaptDataInitialMesh adapt_data(global_domain, initial_refinement_level, initial_mesh_layout);

    /* Build the initial mesh via recursive refinement */
    initial_forest = t8_forest_new_adapt(initial_forest, RefineToInitialMesh, 1, 0, static_cast<void*>(&adapt_data)); 
    #endif
    auto [initial_forest, initial_refinement_level] = cmc::BuildInitialMesh(global_domain, initial_mesh_layout, comm_);

    initial_mesh_.SetMesh(initial_forest);
    initial_mesh_.SetInitialRefinementLevel(initial_refinement_level);
}

std::vector<IndexReduction>
AmrData::UpdateLinearIndicesToTheInitialMesh()
{
    const t8_locidx_t num_local_elements = initial_mesh_.GetNumberLocalElements();

    const int initial_refinement_level = initial_mesh_.GetInitialRefinementLevel();

    const t8_eclass_t eclass = t8_forest_get_eclass(initial_mesh_.GetMesh(), 0);

    /* Get the scheme of the forest's only tree */
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (initial_mesh_.GetMesh(), eclass);

    /* Get the element scheme */
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    const int num_children = ts_c->t8_element_num_children(t8_forest_get_element_in_tree(initial_mesh_.GetMesh(), 0, 0));
    
    MortonIndex linear_index_start_elem = GetMortonIndexOnLevel(t8_forest_get_element_in_tree(initial_mesh_.GetMesh(), 0, 0),
                                                   ts_c, t8_eclass_to_dimension[eclass], initial_refinement_level);
    
    std::vector<IndexReduction> index_correction;
    index_correction.emplace_back(linear_index_start_elem, linear_index_start_elem);

    MortonIndex skipped_indices = 0;

    bool coarse_element_streak = false;

    for (auto iter = 0; iter < num_local_elements; ++iter)
    {
        t8_element_t* elem = t8_forest_get_element_in_tree(initial_mesh_.GetMesh(), 0, iter);

        if (t8_element_level(ts, elem) != initial_refinement_level)
        {
            /* Get the number of uniform indices which were skipped by this element which lays outside of the domain */
            skipped_indices += std::pow(num_children, initial_refinement_level - t8_element_level(ts, elem)) - 1;
            coarse_element_streak = true;
        } else if (coarse_element_streak)
        {
            MortonIndex uniform_index_of_elem = GetMortonIndexOnLevel(elem, ts_c, t8_eclass_to_dimension[eclass], initial_refinement_level);
            index_correction.emplace_back(uniform_index_of_elem, skipped_indices);
            coarse_element_streak = false;
        }
    }

    return index_correction;
}


void
AmrData::DistributeDataOnInitialMesh()
{
    //TODO: Re-Integrate the parallelization
    #if 0
    cmc_assert(initial_mesh.IsValid());

    DataOffsets offsets = GatherGlobalDataOffsets(initial_mesh_, comm_);

    //TransformToMortonIndices sollte hier hin 

    //GroupVariablesForDistributing();//TODO: implement

    //Call (friend) function supplying the send list for each variable
    int comm_size = 1;

    #ifdef CMC_ENABLE_MPI
    int ret_val = MPI_Comm_size(comm_, &comm_size);
    MPICheckError(ret_val);
    #endif

    std::vector<VariableSendMessage> send_messages;
    send_messages.reserve(comm_size * input_variables_.size());

    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        std::vector<VariableSendMessage> var_send_messages = GatherDataToBeDistributed(*input_var_iter);
        send_messages.insert(std::back_insert_iterator(send_messages), std::make_move_iterator(var_send_messages.begin()), std::make_move_iterator(var_send_messages.end()));
        var_send_messages.clear();
    }

    /* Send the data for all variables */
    for (auto send_iter = send_messages.begin(); send_iter != send_messages.end(); ++send_iter)
    {
        send_iter->SendMessage();
    }

    /* Move the data which remains local */
    //....

    MPIBarrier(...);

    /* Receive the data */
    //....


    //if MPI
    GatherGlobalDataOffsets();

    ComputeMortonIndices();

    // if MPI
    CommunicateNonCompliantData();

    //SortReceivedDataLocally();
    #endif


    // Serial code below!!!

    const MortonIndex start_offset = 0;

    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        input_var_iter->TransformCoordinatesToLinearIndices();
    }

    /* Update the uniform morton indices to linear element IDs corresponding to the initial mesh */
    UpdateLinearIndices local_indices_update(UpdateLinearIndicesToTheInitialMesh());

    std::vector<InputVar> sorted_input_variables;
    sorted_input_variables.reserve(input_variables_.size());

    //Create new input variable with the same type
    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++ input_var_iter)
    {
        sorted_input_variables.push_back(MetaCopy(*input_var_iter));
        sorted_input_variables.back().SetUpFilledVariable(initial_mesh_.GetNumberLocalElements(), input_var_iter->GetMissingValue());
        sorted_input_variables.back().AssignDataAtLinearIndices(*input_var_iter, start_offset, local_indices_update);
    }

    /* Swap the sorted input variables with the initial input variables */
    std::swap(input_variables_, sorted_input_variables);


    //TODO: remove below:

    //std::vector<float> converted_data;
    //const int num_elems = 700*300;
    //converted_data.reserve(num_elems);
    //const float add_offset = 268.252690054217;
    //const float scale_factor = 0.00151686422239082;
    //const short missing_value = -32767;
    //const CmcInputVariable& invar = sorted_input_variables[0].GetInternalVariant();
    //const InputVariable<short>& ac_var = std::get<InputVariable<short>>(invar);

    //for (int j = 0; j < num_elems; ++j)
    //{
    //    if (ac_var[j] == missing_value)
    //    {
    //        converted_data.push_back(static_cast<float>(missing_value));
    //    } else
    //    {
    //        converted_data.push_back(scale_factor * static_cast<float>(ac_var[j]) + add_offset);
    //    }
    //}

    //cmc_assert(sizeof(float) == sizeof(uint32_t));

    //cmc_debug_msg("Binary presentation of all floats (in SFC order)");
    //for (auto ii = converted_data.begin(); ii != converted_data.end(); ++ii)
    //{
    //    float val = *ii;

    //    //float val = *ii;
    //    //char* ptr = std::reinterpret_cast<char*>(&val);
    //    //for (int j = 0; j < sizeof(float); ++j)
    //    //{
    //    //    std::cout << 
    //    //}
    //    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&val)) <<" fuer " << *ii << std::endl;
    //}

    //std::exit(1);
}

void
AmrData::ApplyScalingAndOffset()
{
    for (auto iter = input_variables_.begin(); iter != input_variables_.end(); ++iter)
    {
        iter->ApplyScalingAndOffset();
    }
}

std::vector<ByteVar>
AmrData::GetByteVariablesForCompression()
{
    cmc_assert(!variables_.empty());
    cmc_assert(initial_mesh_.IsValid());
    
    std::vector<ByteVar> byte_variables;
    byte_variables.reserve(variables_.size());

    /* Create a callable transforming the input variables to  compression variables */
    TransformerCompressionToByteVariable transform_to_byte_variable;

    /* Transform each compression variable to a byte variable */
    for (std::vector<cmc::Var>::iterator var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        byte_variables.emplace_back(transform_to_byte_variable(*var_iter));
    }

    /* Delete the input variables (since they are no longer needed) */
    variables_.clear();

    return byte_variables;
}

void
AmrData::SetupVariablesForCompression()
{
    cmc_assert(initial_mesh_.IsValid());

    /* Create a callable transforming the input variables to  compression variables */
    TransformerInputToCompressionVariable transform_to_compression_variable;

    /* Allocate memory for the compression variables */
    variables_.reserve(input_variables_.size());

    /* Transform each InputVariable to a compression variable */
    for (std::vector<cmc::InputVar>::iterator input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        variables_.emplace_back(transform_to_compression_variable(*input_var_iter));
    }

    /* Delete the input variables (since they are no longer needed) */
    input_variables_.clear();

    /* Check whether general compresison criteria have been supplied */
    const auto CapturekErrorCriterionHoldsForAllVariables = kErrorCriterionHoldsForAllVariables;
    auto general_specifications_iter = std::find_if(compression_settings_.GetSpecificationsBegin(), compression_settings_.GetSpecificationsEnd(), [=](const CompressionSpecifications& arg){return arg.GetVariableID() == CapturekErrorCriterionHoldsForAllVariables;});
    const bool are_general_specifications_given = (general_specifications_iter != compression_settings_.GetSpecificationsEnd() ? true : false);

    /* Compelete the setup for each variable with it's comrpession criteria, interpolation function and the inaccuracy container */
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        //TODO: Maybe set general settings for all variables 

        /* Set the compression criteria */
        if (auto var_specifications_iter = std::find_if(compression_settings_.GetSpecificationsBegin(), compression_settings_.GetSpecificationsEnd(), [&](const CompressionSpecifications& arg){return arg.GetVariableID() == var_iter->GetID();});
            var_specifications_iter != compression_settings_.GetSpecificationsEnd())
        {
            var_iter->SetUpCompressionCriteria(*var_specifications_iter);
        } else if (are_general_specifications_given)
        {
            var_iter->SetUpCompressionCriteria(*general_specifications_iter);
        } else 
        {
            cmc_err_msg("There are no compression settings for the variable ('", var_iter->GetName(), "') given.");
        }
        
        /* Set the inaccuracy container up */
        var_iter->SetUpInaccuracyStorage(initial_mesh_.GetNumberLocalElements());

        /* Set the intial mesh */
        var_iter->SetAmrMesh(initial_mesh_);
    }
}

size_t AmrData::GetNumberOfInputVariables() const
{
    return input_variables_.size();
}

size_t AmrData::GetNumberOfCompressionVariables() const
{
    return variables_.size();
}

AdaptData
AmrData::CreateAdaptationData(const CoarseningSample& adaptation_sample, const CompressionMode mode)
{
    return AdaptData(compression_settings_, variables_, adaptation_sample, mode);
}

std::vector<CoarseningSample>
AmrData::RetrieveMeshesToBeCoarsened(const CompressionMode compression_mode) const
{
    cmc_assert(!variables_.empty());

    std::vector<CoarseningSample> mesh_samples_to_be_coarsened;

    switch (compression_mode)
    {
        case CompressionMode::OneForAll:
            /* In case of an One For All compression, all variables are defined and coarsend on the same mesh */
            mesh_samples_to_be_coarsened.emplace_back(CoarseningSample(kMeshCorrespondsToAllVariables));
        break;
        case CompressionMode::OneForOne:
            /* In case of a One For One compression, every variable is coarsened on it's own mesh */
            mesh_samples_to_be_coarsened.reserve(variables_.size());
            for (auto iter = variables_.begin(); iter != variables_.end(); ++iter)
            {
                /* We push back the 'id' of each variable (i.e. 0, 1, 2,...,#num_variables - 1), since each variable has it's own mesh which needs to be coarsened */
                mesh_samples_to_be_coarsened.emplace_back(CoarseningSample(std::distance(variables_.begin(), iter)));
            }
        break;
        default:
            cmc_err_msg("The compression mode has not been specified.");
    }

    return mesh_samples_to_be_coarsened;
}

bool
AmrData::CheckConsistencyOfInputVariables() const
{
    if (!input_variables_.empty())
    {
        /* Get the layout and domain of the first variable */
        const DataLayout& first_layout = input_variables_.front().GetInitialDataLayout();
        const GeoDomain& first_domain = input_variables_.front().GetGlobalDomain();
    
        for (auto input_var_iter = std::next(input_variables_.begin()); input_var_iter != input_variables_.end(); ++input_var_iter)
        {
            if (first_layout != input_var_iter->GetInitialDataLayout())
                return false;
            if (!CompareGeoDomains(first_domain, input_var_iter->GetGlobalDomain()))
                return false;
        }
    }

    return true;
}

bool
AmrData::IsValidForCompression() const
{
    /* Checkt he variables's data and the initial mesh */
    const t8_locidx_t num_elements = initial_mesh_.GetNumberLocalElements();

    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        const cmc::AmrMesh& var_mesh = var_iter->GetAmrMesh();
        if (var_mesh.GetNumberLocalElements() != num_elements)
            return false;

        if (!var_iter->IsValidForCompression())
            return false;
    }

    return true;
}

int
AmrData::GetMpiSize() const
{
    int mpi_size, err;
    err = MPI_Comm_size(comm_, &mpi_size);
    MPICheckError(err);
    return mpi_size;
}

int
AmrData::GetMpiRank() const
{
    int mpi_rank, err;
    err = MPI_Comm_rank(comm_, &mpi_rank);
    MPICheckError(err);
    return mpi_rank;
}

void
AmrData::CompressByAdaptiveCoarsening(const CompressionMode compression_mode)
{
    #ifdef CMC_WITH_T8CODE
    /* Retrieve all forests which are going to be coarsened. In case of a 'One For All'-compression we will just
     * receive a single forest on which all variables are defined. However, in case of a 'One for One'-compression,
     * we will receive multiple forests which need to be coarsened, in particular a single forest for each variable */
    std::vector<CoarseningSample> meshes_to_be_coarsened = RetrieveMeshesToBeCoarsened(compression_mode);

    for (auto current_sample = meshes_to_be_coarsened.begin(); current_sample != meshes_to_be_coarsened.end(); ++current_sample)
    {
        /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
        AdaptData adapt_data = CreateAdaptationData(*current_sample, compression_mode);

        while (adapt_data.IsCompressionProgressing())
        {
            adapt_data.InitializeCompressionIteration();

            t8_forest_t previous_forest = adapt_data.GetCurrentMesh();

            /* Keep the 'previous forest' after the adaptation step */
            t8_forest_ref(previous_forest);

            /* Perform a coarsening iteration */
            t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, adapt_data.GetAdaptationFunction(), 0, 0, static_cast<void*>(&adapt_data));

            /* Interpolate the data from the previous forest to the (new) coarsened forest */
            adapt_data.UpdateCompressionData();

            /* Repartition the mesh as well as the data and receive the newly partitioned forest */
            adapted_forest = adapt_data.RepartitionData(adapted_forest);

            /* Free the former forest */
            t8_forest_unref(&previous_forest);

            adapt_data.SetCurrentMesh(adapted_forest);

            adapt_data.FinalizeCompressionIteration();

        }
    }

    #endif
}


void 
AmrData::DecompressToInitialRefinementLevel(const bool restrict_to_global_domain)
{
    /* For the decompression, we are interpolating the data back onto the initial refinement level.
     * Currently, a constant interpolation is chosen.
     */

    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {   
        const t8_forest_t compressed_forest = var_iter->GetAmrMesh().GetMesh();

        const int initial_refinement_level = var_iter->GetAmrMesh().GetInitialRefinementLevel();

        const DataLayout initial_data_layout = var_iter->GetInitialDataLayout();

        const GeoDomain& var_domain = var_iter->GetGlobalDomain();

        const t8_eclass_t eclass = t8_forest_get_eclass(initial_mesh_.GetMesh(), 0);
        t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (initial_mesh_.GetMesh(), eclass);
        t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);
        const int num_children = ts_c->t8_element_num_children(t8_forest_get_element_in_tree(initial_mesh_.GetMesh(), 0, 0));

        /* An estimation for the memory allocation (we estimate 10% more elements if dummy elements are decompressed as well) */
        const size_t memory_estimation = var_domain.GetNumberReferenceCoordsCovered() + (restrict_to_global_domain ? 0 : 0.1 * var_domain.GetNumberReferenceCoordsCovered());
        var_iter->AllocateForDecompression(memory_estimation);

        const t8_locidx_t num_local_elements = t8_forest_get_local_num_elements(compressed_forest);

        for (t8_locidx_t elem_id = 0; elem_id < num_local_elements; ++elem_id)
        {
            t8_element_t* element = t8_forest_get_element_in_tree(compressed_forest, 0, elem_id);

            /* In case that the dummy elements are excluded from the decompression, we need to ensure that only the amount of elements within the domain 
             * are inserted in the decompressed vector (especially for relatively coarse elements). Therefore, we need a correction.
             */
            const int num_copies = DetermineNumberOfDecompressedElements(restrict_to_global_domain, element, ts, num_children, var_domain, initial_refinement_level, initial_data_layout);

            var_iter->DecompressElementConstantly(elem_id, num_copies);
        }

        var_iter->UpdateDecompressionData();
    }
}

void
AmrData::TransferToDecompressionData(OutputVar& output_variable, const Var& compression_variable, const GeoDomain& resulting_domain)
{
    const t8_forest_t compressed_forest = compression_variable.GetAmrMesh().GetMesh();

    const int initial_refinement_level = compression_variable.GetAmrMesh().GetInitialRefinementLevel();

    const DataLayout initial_data_layout = compression_variable.GetInitialDataLayout();

    const GeoDomain& var_domain = compression_variable.GetGlobalDomain();

    const t8_eclass_t eclass = t8_forest_get_eclass(initial_mesh_.GetMesh(), 0);
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (initial_mesh_.GetMesh(), eclass);

    const t8_locidx_t num_local_elements = t8_forest_get_local_num_elements(compressed_forest);

    for (t8_locidx_t elem_id = 0; elem_id < num_local_elements; ++elem_id)
    {
        t8_element_t* element = t8_forest_get_element_in_tree(compressed_forest, 0, elem_id);

        const Hyperslab element_hyperslab = DetermineHyperslabOfDecompressedElements(element, ts, var_domain, initial_refinement_level, initial_data_layout);

        output_variable.AssignDataToHyperslab(compression_variable.GetValueAt(elem_id), element_hyperslab, resulting_domain);
    }
}

static
int CheckNumberOfCorrespondingVariables(const std::vector<Var>& variables, const int var_id)
{
    int corresponding_variables = 0;
    for (auto var_iter = variables.begin(); var_iter != variables.end(); ++var_iter)
    { 
        if (var_iter->GetID() == var_id)
            ++corresponding_variables;
    } 
    return corresponding_variables;
}

OutputVar
AmrData::DecompressVariable(const int variable_id)
{
    const int mpi_size = GetMpiSize();
    cmc_assert (mpi_size <= 1); //For parallel this does not work (and the end we also assign the global domain)

    const int corresponding_variables = CheckNumberOfCorrespondingVariables(variables_, variable_id);

    auto first_var_iter = std::find_if(variables_.begin(), variables_.end(), [&variable_id](const Var& var){return var.GetID() == variable_id;});
    
    if (first_var_iter == variables_.end())
    {
        cmc_err_msg("There is no variable with the specified ID (", variable_id, ")");
        return OutputVar(CMC_ERR);
    }

    const DataLayout resulting_layout = (first_var_iter->GetPreCompressionLayout() != DataLayout::LayoutUndefined ? first_var_iter->GetPreCompressionLayout() : first_var_iter->GetInitialDataLayout());
    const GeoDomain resulting_domain = first_var_iter->HasSplitDimension() != Dimension::DimensionUndefined ?  ExtendGeoDomain(first_var_iter->GetGlobalDomain(), DimensionInterval(first_var_iter->HasSplitDimension(), 0, corresponding_variables)) : first_var_iter->GetGlobalDomain();
    const int memory_estimation = resulting_domain.GetNumberReferenceCoordsCovered() / mpi_size;

    OutputVar output_variable(variable_id, first_var_iter->GetType(), memory_estimation, first_var_iter->GetName(), resulting_layout, resulting_domain, first_var_iter->GetMissingValue());

    for (auto var_iter = first_var_iter; var_iter != variables_.end(); ++var_iter)
    {
        if (var_iter->GetID() == variable_id)
        {
            TransferToDecompressionData(output_variable, *var_iter, resulting_domain);
        }
    }

    /* Set the global domain as hyperslab */
    output_variable.PushBackHyperslab(TransformGeoDomainToHyperslab(resulting_domain));
    
    return output_variable;
}



//TODO:: in serial allocate a vector holding all possible values of the domain. Get the current compressed element, transfer its anchor coordinates
// to the initial refinement level, obtain the hyperslab this element covers, set the values at the right position (corresponding to the the precompression layout of the data)
// in the output variable vector.
//In parallel it will be much harder, since the the data on each rank will not form a nice hyperslab. (by coarsening the local elements as much as possible, we could obtain, the coarsest possible hyperslab level)
// Moreover, if the variable has been split, each layer may have a different partition. Therefore the data needs be parititon uniformly/consistently for all those (slice-) variables


/* This a raw decompression, just decompression each variable leaving the data within Morton order */
std::vector<OutputVar>
AmrData::SeizeRawDecompressedVariable()
{
    std::vector<OutputVar> output_variables;
    output_variables.reserve(variables_.size());

    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        /* If the variable has not been split */
        TransformerCompressionToOutputVariable transform_to_output_variable;
        output_variables.push_back(std::move(transform_to_output_variable(*var_iter)));
    }

    return output_variables;
}

static
std::vector<std::vector<uint8_t>>
DetermineForestBits(t8_forest_t forest)
{
    t8_forest_ref(forest);

    std::vector<std::vector<uint8_t>> serialized_forest_refinements;

    t8_forest_t forest_adapt = nullptr;

    while (t8_forest_get_local_num_elements(forest) > 1)
    {
        RefinementBits adapt_data(t8_forest_get_local_num_elements(forest));

        forest_adapt = t8_forest_new_adapt(forest, FindRefinementBits, 0, 0, static_cast<void*>(&adapt_data));

        forest = forest_adapt;

        serialized_forest_refinements.push_back(std::move(adapt_data.refinement_indicator));
    }

    if (forest_adapt != nullptr)
    {
        t8_forest_unref(&forest_adapt);
    }

    const size_t level_iterations = serialized_forest_refinements.size();

    size_t num_bytes = 0;
    for (size_t i = 0; i < level_iterations; ++i)
    {
        num_bytes += serialized_forest_refinements[i].size();
    }

    return serialized_forest_refinements;
}



static int
CmcTypeToNcType(const CmcType type)
{
    #ifdef CMC_WITH_NETCDF
    switch (type)
    {
        case CmcType::Int8_t:
            return NC_BYTE;
        break;
        case CmcType::Char:
            return NC_CHAR;
        break;
        case CmcType::Int16_t:
            return NC_SHORT;
        break;
        case CmcType::Int32_t:
            return NC_INT;
        break;
        case CmcType::Float:
            return NC_FLOAT;
        break;
        case CmcType::Double:
            return NC_DOUBLE;
        break;
        case CmcType::Uint8_t:
            return NC_UBYTE;
        break;
        case CmcType::Uint16_t:
            return NC_USHORT;
        break;
        case CmcType::Uint32_t:
            return NC_UINT;
        break;
        case CmcType::Int64_t:
            return NC_INT64;
        break;
        case CmcType::Uint64_t:
            return NC_UINT64;
        break;
        default:
            cmc_err_msg("The supplied CmcType is not convertible to a NC-Type.");
            return CMC_ERR;
    }

    #else
    cmc_err_msg("CMC is not cofigured with netCDF.");
    return CMC_ERR;
    #endif
}


void AmrData::WriteCompressedData(const std::string& file_name) const
{
    /* Create the netCDF file to which the cmompressed data will be written */
    cmc_assert(GetMpiSize() <= 1);

    std::vector<std::vector<uint8_t>>  variables_mesh_refinements;
    variables_mesh_refinements.reserve(variables_.size());

    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        /* Obtain the refinement bits for the compressed data */
        std::vector<std::vector<uint8_t>> levelwise_refinements = DetermineForestBits(var_iter->GetAmrMesh().GetMesh());

        const int num_bytes = std::accumulate(levelwise_refinements.begin(), levelwise_refinements.end(), 0, [](const int& count, const std::vector<uint8_t>& refinement_bits){
            return count + refinement_bits.size();
        });

        variables_mesh_refinements.emplace_back();
        variables_mesh_refinements.back().reserve(num_bytes);

        for (auto lr_iter = levelwise_refinements.rbegin(); lr_iter != levelwise_refinements.rend(); ++lr_iter)
        {
            std::copy(lr_iter->begin(), lr_iter->end(), std::back_inserter(variables_mesh_refinements.back()));
        }
    }

    /* Create an array holding the IDs of the variables and dimensions */
    int variable_ids[variables_.size()];
    int variable_dim_ids[variables_.size()];
    int mesh_ids[variables_.size()];
    int mesh_dim_ids[variables_.size()];

    /* Create a new netCDF File */
    int ncid;
    /* We need to create a CDF-5 file in order to use unisgned data types (needed for the encoded mesh refinements) */
    int err = nc__create(file_name.c_str(), NC_CLOBBER|NC_CDF5, NC_SIZEHINT_DEFAULT, NULL, &ncid);
    NcCheckError(err);

    const std::string var_dim_name{"ce"}; // number of 'compressed elements'
    const std::string mesh_dim_name{"rf"}; // number of encoded 'refinements'

    /* Write a dimension for each variable and it's refinement bits */
    int index = 0;
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter, ++index)
    {
        const std::string var_dimension_name = var_dim_name + std::to_string(index);
        const std::string mesh_dimension_name = mesh_dim_name + std::to_string(index);

        const size_t var_dim_length = (var_iter->size() != 0 ? var_iter->size() : 1);

        /* Define a dimension for the comrpessed elements */
        err = nc_def_dim(ncid, var_dimension_name.c_str(), var_dim_length, variable_dim_ids + index);
        NcCheckError(err);

        const size_t mesh_dim_legnth = (variables_mesh_refinements[index].size() != 0 ? variables_mesh_refinements[index].size() : 1);

        /* Define a dimension for the encoded mesh refinements */
        err = nc_def_dim(ncid, mesh_dimension_name.c_str(), mesh_dim_legnth, mesh_dim_ids + index);
        NcCheckError(err);

        /* Define the variable holding the compressed data */
        err = nc_def_var(ncid, var_iter->GetName().c_str(), CmcTypeToNcType(var_iter->GetType()), 1, &(variable_dim_ids[index]), variable_ids + index);
        NcCheckError(err);

        /** Define attributes for the variable **/

        /* Set the cmc internal ID */
        const int id_ = var_iter->GetID();
        err = nc_put_att(ncid, variable_ids[index], "id", NC_INT, 1, &id_);
        NcCheckError(err);

        /* Set the data layout */
        const int var_layout = static_cast<int>(var_iter->GetInitialDataLayout());
        err = nc_put_att(ncid, variable_ids[index], "layout", NC_INT, 1, &var_layout);
        NcCheckError(err);

        /* Store the dimension legnths of the data */
        const GeoDomain& var_domain = var_iter->GetGlobalDomain();

        if (const int lon = var_domain.GetDimensionLength(Dimension::Lon);
            lon > 1)
        {
            err = nc_put_att(ncid, variable_ids[index], "lon", NC_INT, 1, &lon);
            NcCheckError(err);
        }
        if (const int lat = var_domain.GetDimensionLength(Dimension::Lat);
            lat > 1)
        {
            err = nc_put_att(ncid, variable_ids[index], "lat", NC_INT, 1, &lat);
            NcCheckError(err);
        }
        if (const int lev = var_domain.GetDimensionLength(Dimension::Lev);
            lev > 1)
        {
            err = nc_put_att(ncid, variable_ids[index], "lev", NC_INT, 1, &lev);
            NcCheckError(err);
        }
        if (const int time = var_domain.GetDimensionLength(Dimension::Time);
            time > 1)
        {
            err = nc_put_att(ncid, variable_ids[index], "t", NC_INT, 1, &time);
            NcCheckError(err);
        }

        /* The missing value attribute has to be set */
        var_iter->SetMissingValueInNCFile(ncid, variable_ids[index], CmcTypeToNcType(var_iter->GetType()));
        
        /* Set the global context information */
        if (const int gci = var_iter->GetGlobalContextInformation();
            gci != kGlobalContextInformationNotGiven)
        {
            err = nc_put_att(ncid, variable_ids[index], "gci", NC_INT, 1, &gci);
            NcCheckError(err);
        }
        if (const DataLayout pc_layout = var_iter->GetPreCompressionLayout();
            pc_layout != DataLayout::LayoutUndefined)
        {
            const int pre_compression_layout = static_cast<int>(pc_layout);
            //err = nc_put_att(ncid, variable_ids[index], "pcl", NC_INT, 1, &pre_compression_layout);
            err = nc_put_att(ncid, variable_ids[index], "p", NC_INT, 1, &pre_compression_layout);
            NcCheckError(err);
        }

        /* Define the variable holding the mesh refinement of the compressed forest */
        const std::string mesh_name = "m" + std::to_string(index);
        err = nc_def_var(ncid, mesh_name.c_str(), NC_UBYTE, 1, &(mesh_dim_ids[index]), mesh_ids + index);
        NcCheckError(err);
    }

    /* All dimensions and variables have been defined */
    /* Therefore, we are leaving the 'define-mode' and switch to the data mode */
    err = nc_enddef(ncid);
    NcCheckError(err);

    int var_index = 0;
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter, ++var_index)
    {
        /* Write the variable's data */
        var_iter->SetVariableDataInNCFile(ncid, variable_ids[var_index]);

        /* Write the variable's mesh refinement encoding */
        /* We need to check whether the mesh variable holds any refinements (if not we need to set a placeholder value. Otherwise, we would have had
         * to define an unlimited dimension, which is not the intention if there are no values present)
         */
        if (variables_mesh_refinements[var_index].size() > 0)
        {
            err = nc_put_var(ncid, mesh_ids[var_index], static_cast<void*>(variables_mesh_refinements[var_index].data()));
            NcCheckError(err);
        } else
        {
            uint8_t placeholder_refinement{0};
            err = nc_put_var(ncid, mesh_ids[var_index], static_cast<void*>(&placeholder_refinement));
            NcCheckError(err);
        }
    }

    /* All data has been written. Therefore, the file may be closed */
    err = nc_close(ncid);
    NcCheckError(err);
}


void
AmrData::WriteVTKFilePerVariable(const std::string& file_name) const
{
    //TODO: This does not work in all cases (for the decompression it needs to be updated)

    if (variables_.empty())
    {
        cmc_warn_msg("There are no variables present. Therefore no .vtu-files were written.");
        return;
    }

    /* Create a new vtk field holding the element data arrays */
    t8_vtk_data_field_t *vtk_data = new t8_vtk_data_field_t[1];

    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {

        /* Set the type of the data and pointer to the data */
        std::ignore = snprintf(vtk_data[0].description, 128, "%s", var_iter->GetName().c_str());
        
        vtk_data[0].type = T8_VTK_SCALAR;

        std::vector<double> converted_data = var_iter->GetDataAsDoubleVector();

        vtk_data[0].data = converted_data.data();

        const int vtk_err = t8_forest_vtk_write_file(var_iter->GetAmrMesh().GetMesh(), file_name.c_str(), 1, 1, 1, 1, 0, 1, vtk_data);
        
        if (vtk_err == 0)
            cmc_err_msg("An error occrued during the creation of the t8code-forest vtk file.");
    }
    
    delete[] vtk_data;
}


#if 0





struct LevelwisePrefixes
{
    LevelwisePrefixes(const int num_elements)
    {
        prefix_indicator.reserve(num_elements / 8 + 1);
        prefixes.reserve (num_elements / 4);
    }

    std::vector<uint8_t> prefix_indicator;
    std::vector<uint8_t> prefixes;
};





//vielleciht in Variabte<T> als Funktion
std::vector<
AmrData::GatherPrefixData()
{
    t8_forest_ref(forest);

    std::vector<std::vector<uint8_t>> levelwise_prefix_indicator;
    std::vector<std::vector<uint8_t>> levelwise_prefixes;

    t8_forest_t forest_adapt = nullptr;

    std::vector<Prefix<>> prefix_data = GetVariableDataAsPrefixVector();

    while (t8_forest_get_local_num_elements(forest) > 1)
    {
        LevelwisePrefixes adapt_data(t8_forest_get_local_num_elements(forest));

        forest_adapt = t8_forest_new_adapt(forest, FindRefinementBits, 0, 0, static_cast<void*>(&adapt_data));

        forest = forest_adapt;

        levelwise_prefix_indicator.push_back(std::move(adapt_data.prefix_indicator));
        levelwise_prefixes.push_back(std::move(adapt_data.levelwise_prefixes));
    }

    if (forest_adapt != nullptr)
    {
        t8_forest_unref(&forest_adapt);
    }

    const size_t level_iterations = serialized_forest_refinements.size();

    size_t num_bytes = 0;
    for (size_t i = 0; i < level_iterations; ++i)
    {
        num_bytes += serialized_forest_refinements[i].size();
    }

    return serialized_forest_refinements;
}

#endif



}

