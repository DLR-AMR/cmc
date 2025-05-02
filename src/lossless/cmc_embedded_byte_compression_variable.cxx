#include "lossless/cmc_embedded_byte_compression_variable.hxx"
#include "utilities/cmc_geo_domain.hxx"


#include <array>
#include <tuple>
#include <vector>
#include <utility>

namespace cmc::lossless
{

struct AdaptDataInitialEmbeddedMesh
{
public:
    AdaptDataInitialEmbeddedMesh() = delete;
    AdaptDataInitialEmbeddedMesh(const GeoDomain& domain, const int initial_refinement_lvl, const DataLayout layout)
    : global_domain{domain}, initial_refinement_level{initial_refinement_lvl}, initial_layout{layout}{};

    const GeoDomain& global_domain;
    const int initial_refinement_level;
    const DataLayout initial_layout;
};

static t8_locidx_t
RefineToInitialEmbeddedMesh (t8_forest_t forest,
                             [[maybe_unused]] t8_forest_t forest_from,
                             [[maybe_unused]] t8_locidx_t which_tree,
                             t8_locidx_t lelement_id,
                             t8_eclass_scheme_c * ts,
                             [[maybe_unused]] const int is_family,
                             const int num_elements,
                             t8_element_t * elements[])
{
    AdaptDataInitialEmbeddedMesh* adapt_data = static_cast<AdaptDataInitialEmbeddedMesh*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);
    
    /* Check if the element is already on the initial refinement level (or if it is still coarser) */
    if (t8_element_level(ts, elements[0]) >= adapt_data->initial_refinement_level)
    {
        /* If the element's level is already on the initial refinement level the refinement process stops */
        return 0;
    }

    /* If the element is inside the global domain, it will be refined until the intial refinement level is reached */
    if (IsMeshElementWithinGlobalDomain(elements[0], ts, adapt_data->global_domain, adapt_data->initial_refinement_level, adapt_data->initial_layout))
    {
        return 1;
    } else
    {
        return 0;
    }
}

static t8_eclass_t
DimensionToElementClass(const int dimensionality)
{
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    return (dimensionality == 3 ? T8_ECLASS_HEX : T8_ECLASS_QUAD);
}

static std::array<int, 3>
GetElementAnchorOfElement(const t8_element_t* element, t8_eclass_scheme_c* ts)
{
    std::array<int, 3> element_anchor;

    /* Get the element scheme */
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    /* Receive the integer anchor coordinates of the element */
    ts_c->t8_element_anchor (element, element_anchor.data());

    return element_anchor;
}

static int
CalculateInitialRefinementLevel(const GeoDomain& global_domain)
{
    const size_t max_elem_per_direction = global_domain.GetLargestDimensionLength();
    /* Calculate the induced initial refinement level needed in order to build an enclosing mesh */
    return static_cast<int>(std::ceil(std::log2(max_elem_per_direction) + std::numeric_limits<double>::epsilon()));
}

static std::tuple<t8_forest_t, int, int>
BuildInitialEmbeddedMesh(const GeoDomain& domain, const DataLayout initial_layout, MPI_Comm comm)
{
    /* Get the dimensionality of the domain on which the variable is defined */
    const int dimensionality = domain.GetDimensionality();

    /* Create the cmesh; either a quad or hex tree based on the dimension */
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(DimensionToElementClass(dimensionality), comm, 0, 0, 1);

    /* Determine the initial refinement level for embedding the data */
    const int initial_refinement_level = CalculateInitialRefinementLevel(domain);
    ValidateInitialRefinementLevelForDimensionality(initial_refinement_level, dimensionality);

    /* Construct a forest from the cmesh */
    t8_forest_t initial_forest;
    t8_forest_init(&initial_forest);
    t8_forest_set_cmesh(initial_forest, cmesh, comm);
    t8_forest_set_scheme(initial_forest, t8_scheme_new_default_cxx());
    t8_forest_set_level(initial_forest, 0);
    t8_forest_commit(initial_forest);

    AdaptDataInitialEmbeddedMesh adapt_data(domain, initial_refinement_level, initial_layout);

    t8_forest_t adapted_forest;
    for (int adaptation_steps = 0; adaptation_steps <= initial_refinement_level; ++adaptation_steps)
    {
        t8_forest_init(&adapted_forest);
        t8_forest_set_adapt(adapted_forest, initial_forest, RefineToInitialMesh, 0);
        const int set_partition_for_coarsening = 0; //TODO change to one later
        t8_forest_set_partition(adapted_forest, NULL, set_partition_for_coarsening);
        t8_forest_set_user_data(adapted_forest, static_cast<void*>(&adapt_data));
        t8_forest_commit(adapted_forest);

        initial_forest = adapted_forest;
    }

    return std::make_tuple(initial_forest, initial_refinement_level, dimensionality); 
}

template <typename T>
AmrMesh
AbstractEmbeddedByteCompressionVariable<T>::BuildInitialMesh(const input::Variable<T>& input_variable)
{
    cmc_assert(input_variable.IsValid());

    /* Since all varibales are defined on the same global domain, we are able to take the domain of the first variable */
    const GeoDomain& global_domain = input_variable.GetGlobalDomain();

    /* The mesh layout is equal to the layout of the domains (the exact ordering of the dimensions is allowed to vary) */
    const DataLayout initial_mesh_layout = input_variable.GetInitialDataLayout();

    /* Build the actual embedded mesh based on the given features */
    auto [initial_forest, initial_refinement_level, dimensionality] = BuildInitialMesh(global_domain, initial_mesh_layout, GetMPIComm());

    return AmrMesh(initial_forest, initial_refinement_level, dimensionality);
}


/////////////////////////////////
/////////////////////////////////


template <typename T>
static std::vector<IndexReduction>
AbstractEmbeddedByteCompressionVariable<T>::UpdateLinearIndicesToTheInitialMesh()
{
    const t8_locidx_t num_local_elements = mesh_.GetNumberLocalElements();

    const int initial_refinement_level = mesh_.GetInitialRefinementLevel();

    const t8_eclass_t eclass = t8_forest_get_eclass(mesh_.GetMesh(), 0);

    /* Get the scheme of the forest's only tree */
    t8_eclass_scheme_c* ts =  t8_forest_get_eclass_scheme (mesh_.GetMesh(), eclass);

    /* Get the element scheme */
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    const t8_locidx_t first_ltree_id = 0;

    const int num_children = ts_c->t8_element_num_children(t8_forest_get_element_in_tree(mesh_.GetMesh(), first_ltree_id, 0));

    const MortonIndex linear_index_start_elem = GetMortonIndexOnLevel(t8_forest_get_element_in_tree(mesh_.GetMesh(), first_ltree_id, 0),
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
        const t8_element_t* elem = t8_forest_get_element_in_tree(mesh_.GetMesh(), 0, iter);

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

template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SortInitialDataIntoVariables(input::Variable<T>& input_variable, const std::vector<VariableRecvMessage>& messages)
{
    /* Update the uniform morton indices to linear element IDs corresponding to the initial mesh */
    UpdateLinearIndices local_indices_update(UpdateLinearIndicesToTheInitialMesh());

    /* Setup the variable with the right amount of global elements */
    input_variable.SetUpFilledVariable(initial_mesh_.GetNumberLocalElements(), input_var_iter->GetMissingValue());

    /* Iterate over all messages and assign their data to the correct variables */
    for (auto msg_iter = messages.begin(); msg_iter != messages.end(); ++msg_iter)
    {
        cmc_assert(input_variable.GetInternID() == msg_iter->GetVariableID());
        if (input_variable.GetInternID() == msg_iter->GetVariableID())
        {
            /* Assign the data from the messag at the right position within the variable */
            input_variable.AssignDataAtLinearIndices(*msg_iter, local_indices_update);
        } else
        {
            cmc_warn_msg("A message has been received, but the variable ID does not correspond to any (local) variable.");
        }
    }
}


template <typename T>
[[nodiscard]]
std::pair<std::vector<VariableSendMessage>, std::vector<MPI_Request>>
AbstractEmbeddedByteCompressionVariable<T>::SendInitialData(input::Variable<T>& input_variable)
{
    cmc_assert(input_variable.GetActiveDataFormat());

    int comm_size{1};

    const int ret_val = MPI_Comm_size(GetMPIComm(), &comm_size);
    MPICheckError(ret_val);

    /* Gather the partitioning of the initial mesh */
    const DataOffsets offsets = GatherGlobalDataOffsets(mesh_, GetMPIComm());


    /* Inquire all data that has to be sent */
    std::vector<VariableSendMessage> send_messages;

    input_variable.GatherDistributionData(offsets, send_messages);

    /* A vector collecting all send requests */
    std::vector<MPI_Request> send_requests;
    send_requests.reserve(2 * send_messages.size());

    /* Send all messages */
    for (auto sm_iter = send_messages.begin(); sm_iter != send_messages.end(); ++sm_iter)
    {
        /* Send the message and receive the returned requests */
        auto [req_morton_ids, req_data] = sm_iter->Send(GetMPIComm());

        /* Store the requests */
        send_requests.push_back(std::move(req_morton_ids));
        send_requests.push_back(std::move(req_data));
    }

    return std::make_pair(std::move(send_messages), std::move(send_requests));
}

template <typename T>
std::vector<VariableRecvMessage>
AbstractEmbeddedByteCompressionVariable<T>::ReceiveInitialData(input::Variable<T>& input_variable)
{
    /* Receive all messages */
    bool are_messages_incoming{true};

    std::vector<VariableRecvMessage> recv_messages;

    /* Message receiving loop */
    while(are_messages_incoming)
    {
        int message_flag{0};
        MPI_Status probe_status;

        /* Check if there are any messages incoming */
        int err = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, GetMPIComm(), &message_flag, &probe_status);
        MPICheckError(err);

        /* Get the message if there is one */
        if (message_flag != 0)
        {
            /* Get the necessary information for retrieving the message */
            const int source = probe_status.MPI_SOURCE;
            const int tag = probe_status.MPI_TAG;
            const int variable_id = GetIDFromTag(tag);
            const CmcType var_type = input_variable.GetType();
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
                void* dataptr = msg_iter->GetInitialDataPtr();
                /* Receive the actual data from this process */
                err = MPI_Recv(dataptr, num_elems, data_type, source, tag, GetMPIComm(), MPI_STATUS_IGNORE);
                MPICheckError(err); 
            } else
            {
                void* miptr = msg_iter->GetInitialMortonIndicesPtr();
                /* Otherwise, we receive the sent Morton indices */
                err = MPI_Recv(miptr, num_elems, MPI_MORTON_INDEX_T, source, tag, GetMPIComm(), MPI_STATUS_IGNORE);
                MPICheckError(err);
            } 
        } else
        {
            /* If all messages have been processed, we break the receiving loop */
            are_messages_incoming = false;
        }
    }

    return recv_messages;
}

/* Sorting the initial data only locally */
template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::SortLocalDataOnInitialMesh(input::Variable<T>& input_variable)
{
    /* Update the uniform morton indices to linear element IDs corresponding to the initial mesh */
    UpdateLinearIndices local_indices_update(UpdateLinearIndicesToTheInitialMesh());

    input::Variable<T> sorted_input_variable = input::HollowCopy(input_variable);

    sorted_input_variable.SetUpFilledVariable(mesh_.GetNumberLocalElements(), input_variable.GetMissingValue());
    sorted_input_variable.AssignDataAtLinearIndices(input_variable, local_indices_update);

    /* Swap the sorted input variables with the initial input variables */
    std::swap(input_variable, sorted_input_variable);
}


template <typename T>
void
AbstractEmbeddedByteCompressionVariable<T>::DistributeDataOnInitialMesh(input::Variable<T>& input_variable)
{
    #ifdef CMC_ENABLE_MPI

    cmc_assert(mesh_.IsValid());

    /* Gather and send all the data that has to be communicated. 
     * The send messages has to be returned in order to no get deallocated before 
     * the actual sending has happend. Therefore, we keep them here 'alive' until 
     * all send_requests have been completed */
    auto [send_messages, send_requests] = SendInitialData(input_variable);

    /* Wait until all messages have been staged */
    int err = MPI_Barrier(comm_);
    MPICheckError(err);

    /* After all messages have been staged, we will receive them */
    const std::vector<VariableRecvMessage> received_messages = ReceiveInitialData(input_variable);

    /* Sort the data accordingly to the Morton indices */
    SortInitialDataIntoVariables(input_variable, received_messages);

    /* Wait until any send messages have been completed */
    err = MPI_Waitall(send_requests.size(), send_requests.data(), MPI_STATUS_IGNORE);
    MPICheckError(err);

    #else
    /* Call a locally sorting function which setups the data compliant to the Morton order */
    cmc_assert(mesh.IsValid());
    SortLocalDataOnInitialMesh();
    #endif
}

}
