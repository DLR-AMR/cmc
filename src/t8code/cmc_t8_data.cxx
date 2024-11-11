#include "t8code/cmc_t8_data.hxx"
#include "t8code/cmc_t8_mpi.hxx"
#include "t8code/cmc_t8_adapt_callbacks.hxx"
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
#include <t8_forest/t8_forest.h>
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
AmrData::SetCompressionSettings(const CompressionSettings& settings)
{
    compression_settings_ = settings;
}

void
AmrData::SetCompressionSettings(CompressionSettings&& settings)
{
    compression_settings_ = std::move(settings);
}

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

    /* The mesh layout is equal to the layout of the domains (the exact ordering of the dimensions is allowed to vary) */
    const DataLayout initial_mesh_layout = input_variables_.front().GetInitialDataLayout();

    auto [initial_forest, initial_refinement_level, dimensionality] = cmc::BuildInitialMesh(global_domain, initial_mesh_layout, comm_);

    initial_mesh_.SetMesh(initial_forest);
    initial_mesh_.SetInitialRefinementLevel(initial_refinement_level);
    initial_mesh_.SetDimensionality(dimensionality);
    
    /* In most cases, there will be dummy elements present within the mesh */
    initial_mesh_.IndicateWhetherDummyElementsArePresent(true);
}

void
AmrData::SetInitialMesh(const AmrMesh& mesh)
{
    cmc_assert(mesh.AreDummyElementsPresent() == false);

    /* Set the initial mesh */
    initial_mesh_ = mesh;

    /* Since we receive a mesh, we expect that there are no dummy elements present */
    initial_mesh_.IndicateWhetherDummyElementsArePresent(false);
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

    const t8_locidx_t first_ltree_id = 0;

    const int num_children = ts_c->t8_element_num_children(t8_forest_get_element_in_tree(initial_mesh_.GetMesh(), first_ltree_id, 0));

    const MortonIndex linear_index_start_elem = GetMortonIndexOnLevel(t8_forest_get_element_in_tree(initial_mesh_.GetMesh(), first_ltree_id, 0),
                                                   ts_c, t8_eclass_to_dimension[eclass], initial_refinement_level);
    
    /* Locally, this update function reduces the global_index to zero for the first local element */
    std::vector<IndexReduction> index_correction;
    index_correction.emplace_back(linear_index_start_elem, linear_index_start_elem);

    /* All elements that are not holding data need to be accumulated in order to be subtracted additionally for the local index correction */
    MortonIndex skipped_indices = linear_index_start_elem;

    bool coarse_element_streak = false;

    /* Iterate through all local elements and find how the global Morton indices need to be adjusted in order to comply the local data ordering */
    for (auto iter = 0; iter < num_local_elements; ++iter)
    {
        const t8_element_t* elem = t8_forest_get_element_in_tree(initial_mesh_.GetMesh(), 0, iter);

        if (t8_element_level(ts, elem) != initial_refinement_level)
        {
            /* Get the number of uniform indices which were skipped by this element which lays outside of the domain */
            skipped_indices += std::pow(num_children, initial_refinement_level - t8_element_level(ts, elem)) - 1;
            coarse_element_streak = true;
        } else if (coarse_element_streak)
        {
            /* We accumulate the amount of skipped indices (with regard to the initial refinement level) and store the offset 
             * once we have reached again an element on the initial refinement level */
            const MortonIndex uniform_index_of_elem = GetMortonIndexOnLevel(elem, ts_c, t8_eclass_to_dimension[eclass], initial_refinement_level);
            index_correction.emplace_back(uniform_index_of_elem, skipped_indices);
            coarse_element_streak = false;
        }
    }

    /* Return the correction scheme for the global indices */
    return index_correction;
}

void
AmrData::CreateInternIDsForRedistribution()
{
    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        if (input_var_iter->GetGlobalContextInformation() != kNoGlobalContext)
        {
            /* If the variable is part of a higher dimensional one and resembles a subband of the variable */
            const int intern_id = CreateSubbandID(input_var_iter->GetID(), input_var_iter->GetGlobalContextInformation());
            input_var_iter->SetInternID(intern_id);
        } else
        {
            /* If the variable has not been split */
            const int intern_id = input_var_iter->GetID();
            input_var_iter->SetInternID(intern_id);
        }
    }
}


/* TODO: Maybe pass output variables as in/out reference in order to be sure, that the send messages are not deallocated */
[[nodiscard]]
std::pair<std::vector<VariableSendMessage>, std::vector<MPI_Request>>
AmrData::SendInitialData()
{
    int comm_size{1};

    const int ret_val = MPI_Comm_size(comm_, &comm_size);
    MPICheckError(ret_val);

    /* Gather the partitioning of the initial mesh */
    const DataOffsets offsets = GatherGlobalDataOffsets(initial_mesh_, comm_);
    cmc_debug_msg("Offsets inquired");

    /* Inquire all data that has to be sent */
    std::vector<VariableSendMessage> send_messages;
    send_messages.reserve(comm_size * input_variables_.size());

    /* Gather the data thas has to be communicated for all variables */
    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        GatherDistributionData(*input_var_iter, offsets, send_messages);
    }

    cmc_debug_msg("Size of send_messages: ", send_messages.size());


    /* A vector collecting all send requests */
    std::vector<MPI_Request> send_requests;
    send_requests.reserve(2 * send_messages.size());

    /* Send all messages */
    for (auto sm_iter = send_messages.begin(); sm_iter != send_messages.end(); ++sm_iter)
    {
        /* Send the message and receive the returned requests */
        auto [req_morton_ids, req_data] = sm_iter->Send(comm_);

        /* Store the requests */
        send_requests.push_back(std::move(req_morton_ids));
        send_requests.push_back(std::move(req_data));
    }

    return std::make_pair(std::move(send_messages), std::move(send_requests));
}

std::vector<VariableRecvMessage>
AmrData::ReceiveInitialData()
{
    /* Receive all messages */
    bool are_messages_incoming{true};

    std::vector<VariableRecvMessage> recv_messages;
    recv_messages.reserve(input_variables_.size() * 2);

    /* Message receiving loop */
    while(are_messages_incoming)
    {
        int message_flag{0};
        MPI_Status probe_status;

        /* Check if there are any messages incoming */
        int err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm_, &message_flag, &probe_status);
        MPICheckError(err);

        /* Get the message if there is one */
        if (message_flag != 0)
        {
            /* Get the necessary information for retrieving the message */
            const int source = probe_status.MPI_SOURCE;
            const int tag = probe_status.MPI_TAG;
            const int variable_id = GetIDFromTag(tag);
            const CmcType var_type = GetDataTypeFromVariableViaInternID(input_variables_, variable_id);
            const MPI_Datatype data_type = ConvertCmcTypeToMPIType(var_type);

            /* Get the number of elements from this message */
            int num_elems{0};
            if (IsMessageADataMessage(tag))
            {
                /* If the data contains actual data values */
                err = MPI_Get_count(&probe_status, data_type, &num_elems);
                MPICheckError(err);
            } else
            {
                /* If the data contains Morton indices */
                err = MPI_Get_count(&probe_status, MPI_MORTON_INDEX_T, &num_elems);
                MPICheckError(err);
            }

            cmc_debug_msg("In Receive Initial Data: Id is: ", variable_id);
            /** We only do receive two messages from each process that sends data, since it is collected prior to the communication ,
              * Therefore, we need to check if we have already received elements from the rank this messag originated from. */
            std::vector<cmc::VariableRecvMessage>::iterator msg_iter = std::find_if(recv_messages.begin(), recv_messages.end(),
                                                                                    [&source, &variable_id](auto& msg){return (source == msg.GetSendingRank() && variable_id == msg.GetVariableID());});
            
            /* We receive (potentially) two messages from each process for each variable */

            /* Check if there has been a message found coming from this rank */
            if (msg_iter == recv_messages.end())
            {
                /* We have not yet received a message from this rank */
                /* Allocate a variable message */
                recv_messages.emplace_back(var_type, CreateVariableMessage(source, variable_id, var_type, num_elems));

                /* Assign the newly created message to the iterator */
                msg_iter = std::prev(recv_messages.end());
            }

            /* Check if it is a Morton indices message or the actual data message */
            if (IsMessageADataMessage(tag))
            {
                cmc_debug_msg("Trying to Received ", num_elems, " data values from rank ", source, " for variable ", variable_id);
                void* dataptr = msg_iter->GetInitialDataPtr();
                /* Receive the actual data from this process */
                err = MPI_Recv(dataptr, num_elems, data_type, source, tag, comm_, MPI_STATUS_IGNORE);
                MPICheckError(err); 

                cmc_debug_msg("Received ", num_elems, " data values from rank ", source, " for variable ", variable_id);
            } else
            {
                cmc_debug_msg("Trying to Received ", num_elems, " Morton indices from rank ", source, " for variable ", variable_id);
                void* miptr = msg_iter->GetInitialMortonIndicesPtr();
                /* Otherwise, we receive the sent Morton indices */
                err = MPI_Recv(miptr, num_elems, MPI_MORTON_INDEX_T, source, tag, comm_, MPI_STATUS_IGNORE);
                MPICheckError(err);

                cmc_debug_msg("Received ", num_elems, " Morton indices from rank ", source, " for variable ", variable_id);
            } 
        } else
        {
            /* If all messages have been processed, we break the receiving loop */
            are_messages_incoming = false;
        }
    }

    return recv_messages;
}

void
AmrData::SortInitialDataIntoVariables(const std::vector<VariableRecvMessage>& messages)
{
    /* Update the uniform morton indices to linear element IDs corresponding to the initial mesh */
    UpdateLinearIndices local_indices_update(UpdateLinearIndicesToTheInitialMesh());

    /* Setup the variable with the right amount of global elements */
    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++ input_var_iter)
    {
        input_var_iter->SetUpFilledVariable(initial_mesh_.GetNumberLocalElements(), input_var_iter->GetMissingValue());
    }
    cmc_debug_msg("\n\n\nSizeof messages: ", messages.size(), "\n");
    /* Iterate over all messages and assign their data to the correct variables */
    for (auto msg_iter = messages.begin(); msg_iter != messages.end(); ++msg_iter)
    {
        /* Find the variable to which the current message is assigned */
        auto var_iter = std::find_if(input_variables_.begin(), input_variables_.end(), [&msg_iter](auto& var){return var.GetInternID() == msg_iter->GetVariableID();});
        cmc_debug_msg("\n\n\nSort in intial data\n\nName: ", var_iter->GetName(), ", id: ", var_iter->GetID(), ", global context: ", var_iter->GetGlobalContextInformation());
        /* All processes should have prior knowledge to all variables. Therefore, a variable with the ID from the message has to be found locally */
        if (var_iter == input_variables_.end())
        {
            cmc_err_msg("A message has been received, but the variable ID does not correspond to any (local) variable.");
        }

        /* Assign the data from the messag at the right position within the variable */
        var_iter->AssignDataAtLinearIndices(*msg_iter, local_indices_update);
    }
}

/* Sorting the initial data only locally */
void
AmrData::SortLocalDataOnInitialMesh()
{
    /* Update the uniform morton indices to linear element IDs corresponding to the initial mesh */
    UpdateLinearIndices local_indices_update(UpdateLinearIndicesToTheInitialMesh());

    std::vector<InputVar> sorted_input_variables;
    sorted_input_variables.reserve(input_variables_.size());

    /* Create new variables and assign the sorted data to them */
    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        sorted_input_variables.emplace_back(MetaCopy(*input_var_iter));
        sorted_input_variables.back().SetUpFilledVariable(initial_mesh_.GetNumberLocalElements(), input_var_iter->GetMissingValue());
        sorted_input_variables.back().AssignDataAtLinearIndices(*input_var_iter, local_indices_update);
    }

    /* Swap the sorted input variables with the initial input variables */
    std::swap(input_variables_, sorted_input_variables);
}

void
AmrData::DistributeDataOnInitialMesh()
{
    #ifdef CMC_ENABLE_MPI
    cmc_debug_msg("in DistributeDataOnInitialMesh\n\n");
    cmc_assert(initial_mesh_.IsValid());

    /* We create intern IDs that are uniquely identifiable for the parallel communication */
    CreateInternIDsForRedistribution();

    cmc_debug_msg("Size of input variables: ", input_variables_.size());

    /* Transform all coordinates to Morton indices */
    for (auto input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        input_var_iter->TransformCoordinatesToLinearIndices();
    }

    /* Gather and send all the data that has to be communicated. 
     * The send messages has to be returned in order to no get deallocated before 
     * the actual sending has happend. Therefore, we keep them here 'alive' until 
     * all send_requests have been completed */
    auto [send_messages, send_requests] = SendInitialData();
    cmc_debug_msg("send_messages.size() = ", send_messages.size());
    /* Wait until all messages have been staged */
    int err = MPI_Barrier(comm_);
    MPICheckError(err);

    /* After all messages have been staged, we will receive them */
    const std::vector<VariableRecvMessage> received_messages = ReceiveInitialData();
    cmc_debug_msg("received_messages.size() = ", received_messages.size());
    /* Sort the data accordingly to the Morton indices */
    SortInitialDataIntoVariables(received_messages);

    /* Wait until any send messages have been completed */
    err = MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUS_IGNORE);
    MPICheckError(err);

    #else
    /* Call a locally sorting function which setups the data compliant to the Morton order */
    cmc_assert(initial_mesh_.IsValid());
    SortLocalDataOnInitialMesh();
    #endif
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
AmrData::TransformInputToCompressionVariables()
{
    cmc_assert(initial_mesh_.IsValid());

    /* Check if the input variables have already been transformed */
    if (input_variables_.empty()) {return;}

    /* Create a callable transforming the input variables to  compression variables */
    TransformerInputToCompressionVariable transform_to_compression_variable;

    /* Allocate memory for the compression variables */
    variables_.reserve(input_variables_.size());

    /* Transform each InputVariable to a compression variable */
    for (std::vector<cmc::InputVar>::iterator input_var_iter = input_variables_.begin(); input_var_iter != input_variables_.end(); ++input_var_iter)
    {
        /* Store the transformed compression variable */
        variables_.emplace_back(transform_to_compression_variable(*input_var_iter));

        /* Set the intial mesh */
        variables_.back().SetAmrMesh(initial_mesh_);
    }

    /* Delete the input variables (since they are no longer needed) */
    input_variables_.clear();

    /* Allocate memory for the indication bits resulting from the adaptive caorsening procedure */
    ac_indications_.reserve(variables_.size());
}

void
AmrData::SetupVariablesForCompression()
{
    TransformInputToCompressionVariables();

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
            for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
            {
                /* We push back the 'id' of each variable (i.e. 0, 1, 2,...,#num_variables - 1), since each variable has it's own mesh which needs to be coarsened */
                //mesh_samples_to_be_coarsened.emplace_back(CoarseningSample(std::distance(variables_.begin(), iter)));
                mesh_samples_to_be_coarsened.emplace_back(CoarseningSample(var_iter->GetInternalID()));
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
    /* Check the variables's data and the initial mesh */
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
            cmc_debug_msg("Befroe adapt compression");
            /* Perform a coarsening iteration */
            t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, adapt_data.GetAdaptationFunction(), 0, 0, static_cast<void*>(&adapt_data));
            cmc_debug_msg("Befroe update compression");
            /* Interpolate the data from the previous forest to the (new) coarsened forest */
            adapt_data.UpdateCompressionData();
            cmc_debug_msg("Befroe repartition");
            /* Repartition the mesh as well as the data and receive the newly partitioned forest */
            adapted_forest = adapt_data.RepartitionData(adapted_forest);
            cmc_debug_msg("after repartition");
            /* Free the former forest */
            t8_forest_unref(&previous_forest);
            cmc_debug_msg("Befroe set partitioned mesh");
            adapt_data.SetCurrentMesh(adapted_forest);

            adapt_data.FinalizeCompressionIteration();
        }

        ac_indications_.emplace_back(adapt_data.TransferIndicationBits());
        }

    #endif
}

[[nodiscard]] std::vector<AdaptiveCoarseningIndications>
AmrData::TransferIndicationBits()
{
    std::vector<AdaptiveCoarseningIndications> temp_ = std::move(ac_indications_);
    ac_indications_.clear();
    return temp_;
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
            const t8_element_t* element = t8_forest_get_element_in_tree(compressed_forest, 0, elem_id);

            /* In case that the dummy elements are excluded from the decompression, we need to ensure that only the amount of elements within the domain 
             * are inserted in the decompressed vector (especially for relatively coarse elements). Therefore, we need a correction.
             */
            const int num_copies = DetermineNumberOfDecompressedElements(restrict_to_global_domain, element, ts, num_children, var_domain, initial_refinement_level, initial_data_layout);

            var_iter->DecompressElementConstantly(elem_id, num_copies);
        }

        var_iter->UpdateDecompressionData();
    }
    #if 0
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        std::vector<double> decompressed_data = var_iter->GetDataAsDoubleVector();
        cmc_debug_msg("Size of decompressed data vector: ", decompressed_data.size());
        std::vector<short> finit_data(decompressed_data.size(), 0);
        std::string fname = "era5_initial_t_sfc_ordered_" + std::to_string(var_iter->GetGlobalContextInformation());
        const int num_elems_initial_mesh = 1440*721;
        FILE* file = fopen(fname.c_str(), "rb");
        fread(finit_data.data(), sizeof(short), num_elems_initial_mesh, file);
        fclose(file);
        cmc_debug_msg("Size decompressed data var: ", decompressed_data.size(), " und Soll_wert: 1440*721 = ", 1440*721);
        int larger_errs = 0;
        double max_err = 0.0;
        double dmax_err = 0.0;
        double scale_factor = 0.00200295932540368;
        double add_offset = 250.108564255201;
        for  (size_t init_acces = 0; init_acces < decompressed_data.size(); ++init_acces)
        {

            double ferr = std::abs(static_cast<double>(finit_data[init_acces]) - decompressed_data[init_acces]);
            double dferr = std::abs((scale_factor * static_cast<double>(finit_data[init_acces]) + add_offset) - (scale_factor * decompressed_data[init_acces] + add_offset));

            #if 0
            if (ferr > 30000.0)
            {
                ++larger_errs;
            }
            else 
            #endif

            if (ferr > max_err)
            {
                max_err = ferr;
            }
            if (dferr > dmax_err)
            {
                dmax_err = dferr;
            }
        }

        cmc_debug_msg("For variable ", var_iter->GetGlobalContextInformation(), " the CHAR max error is: ", max_err);
        cmc_debug_msg("For variable ", var_iter->GetGlobalContextInformation(), " the DOUBLE max error is: ", dmax_err);
        //cmc_debug_msg("It had larger erroros: ", larger_errs);
    }
    #endif

    #if 0
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        std::vector<double> decompressed_data = var_iter->GetDataAsDoubleVector();
        cmc_debug_msg("Size of decompressed data vector: ", decompressed_data.size());

        std::vector<float> finit_data(decompressed_data.size(), 0);
        std::string fname = "era5_initial_t_sfc_ordered_" + std::to_string(var_iter->GetGlobalContextInformation());
        const int num_elems_initial_mesh = 1440*721;
        FILE* file = fopen(fname.c_str(), "rb");
        fread(finit_data.data(), sizeof(float), num_elems_initial_mesh, file);
        fclose(file);

        double max_err = 0.0;
        //const double scale = 9.25404335362387E-08;
        //const double offset = 0.00790779508439177;

        
        for (size_t i = 0; i < decompressed_data.size(); ++i)
        {
            double fdata = static_cast<double>(finit_data[i]);

            float rel_err = std::abs(fdata - decompressed_data[i]) / std::abs(fdata);
            //if (rel_err > 2.0)
            //{
            //    cmc_debug_msg("Error für Index ", i, " greater than limit: ", rel_err);
            //}
            if (rel_err > max_err)
            {
                max_err = rel_err;
            }
        }
        cmc_debug_msg("For variable ", var_iter->GetGlobalContextInformation(), " the relative max error is: ", max_err);
    }
    #endif
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
        const t8_element_t* element = t8_forest_get_element_in_tree(compressed_forest, 0, elem_id);

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

    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter)
    {
        cmc_debug_msg("Id of variable: ", var_iter->GetID());
    }
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
    int i = 0;
    for (auto var_iter = variables_.begin(); var_iter != variables_.end(); ++var_iter, ++i)
    {

        /* Set the type of the data and pointer to the data */
        std::ignore = snprintf(vtk_data[0].description, 128, "%s_%d", var_iter->GetName().c_str(), i);
        
        vtk_data[0].type = T8_VTK_SCALAR;

        std::vector<double> converted_data = var_iter->GetDataAsDoubleVector();

        vtk_data[0].data = converted_data.data();

        const int vtk_err = t8_forest_vtk_write_file(var_iter->GetAmrMesh().GetMesh(), (file_name + "_" + var_iter->GetName()).c_str(), 0, 1, 0, 0, 0, 1, vtk_data);
        
        if (vtk_err == 0)
            cmc_err_msg("An error occrued during the creation of the t8code-forest vtk file.");
    }
    
    delete[] vtk_data;
}

}

