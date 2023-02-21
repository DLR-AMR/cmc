#include "cmc_t8_mpi.h"
#include "utilities/cmc_constants_definitions.h"
#include "mpi/cmc_mpi.h"

enum DATA_DISTRIBUTION {DISTRIBUTION_UNDEFINED = 0, CMC_SERIAL_DATA, CMC_LINEARIZED, CMC_BLOCKED, CMC_ZCURVE};

struct cmc_mpi_reading_data_distribution
{
    cmc_mpi_data_distribution(){};
    cmc_mpi_data_distribution(std::vector<uint32_t>&& _start_values, std::vector<uint32_t>&& _count_values, DATA_DISTRIBUTION _distribution)
    : start_values{_start_values}, count_values{_count_values}, distribution{_distribution}{};
    ~cmc_mpi_data_distribution(){};

    std::vector<uint32_t> start_values;
    std::vector<uint32_t> count_values;
    std::vector<std::vector<uint32_t>> id_blocks;

    DATA_DISTRIBUTION distribution{DISTRIBUTION_UNDEFINED};
};

struct cmc_coordinate
{
    std::tuple<uint32_t, uint32_t, uint32_t>> cartesian_coordinate;
    uint64_t morton_index;
    uint64_t linear_id;
}
struct cmc_coordinates
{
    std::vector<std::tuple<uint32_t, uint32_t, uint32_t>> cartesian_coordinates;
    std::vector<uint64_t> morton_indices;
    std::vector<uint32_t> linear_ids;
};

struct cmc_mpi_data_distribution_info
{
    int initial_refinement_level{0};
    int data_dimensionality{0};
    uint32_t local_mesh_elements{0};
    uint64_t global_mesh_elements{0};
    std::vector<uint32_t> start_values;
    std::vector<uint32_t> count_values;
    std::vector<cmc_coordinates> offsets;
};

struct cmc_mpi_transfer_data
{
    int partner_rank{-1};
    size_t num_elements{0};
    std::vector<var_dynamic_array_t*> data;
    std::vector<cmc_coordinates> coordinates;
    DATA_DISTRIBUTION distribution{DISTRIBUTION_UNDEFINED};
};


static int
cmc_t8_elem_in_storage_domain(const int dim, const int initial_refinement_lvl, const t8_element_t* element, t8_eclass_scheme_c* ts,
                              const int x_min, const int y_min, const int z_min, const int x_max, const int y_max, const int z_max)
{
    #if defind CMC_WITH_T8CODE && defined CMC_ENABLE_MPI

    int element_anchor[3];
    t8_default_scheme_common_c* ts_c = static_cast<t8_default_scheme_common_c*>(ts);

    /* Maximum refinement level depending on the dimension of the data */
    const int element_anchor_max_lvl{(dim == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL)};
    /* Receive the integer anchor coodinates of the element */ 
    ts_c->t8_element_anchor (element, element_anchor);

    /* Transform coordinates into the range of the initial-refinement-level coordinate values */
    for (int index{0}; index < dim; ++index)
    {
        element_anchor[index] >>= (element_anchor_max_lvl - initial_refinement_lvl);
    }

    if (dim == 2)
    {
        if (element_anchor[0] >= x_min && element_anchor[0] < x_max &&
            element_anchor[1] >= y_min && element_anchor[1] < y_max)
        {
            /* The 2D element is inside the "lat x lon" data */
            return 1;
        }
    }
    else
    {
        if (element_anchor[0] >= x_min && element_anchor[0] < x_max &&
            element_anchor[1] >= y_min && element_anchor[1] < y_max &&
            element_anchor[2] >= z_min && element_anchor[2] < z_max)
        {
            /* The 3D is inside the "lat x lon x lev" data */
            return 1;
        }
    }
    return 0;
    #else
    return CMC_ERR;
    #endif
}

/* Transform a cartesian coordinate to a morton index */
static uint64_t
cmc_linear_id_to_zcurve_id(const int dim, uint64_t x, uint64_t y, uint64_t z)
{
    cmc_assert(dim == 2 || dim == 3);

    /* Prepare the coordinates for interleaving */
    if (dim == 2)
    {
        /* 2D case */
        x = (x | (x << 16)) & 0x0000ffff0000ffff;
        x = (x | (x << 8)) & 0x00ff00ff00ff00ff;
        x = (x | (x << 4)) & 0x0f0f0f0f0f0f0f0f;
        x = (x | (x << 2)) & 0x3333333333333333;
        x = (x | (x << 1)) & 0x5555555555555555;

        y = (y | (y << 16)) & 0x0000ffff0000ffff;
        y = (y | (y << 8)) & 0x00ff00ff00ff00ff;
        y = (y | (y << 4)) & 0x0f0f0f0f0f0f0f0f;
        y = (y | (y << 2)) & 0x3333333333333333;
        y = (y | (y << 1)) & 0x5555555555555555;

        /* Interleave the two integers */
        return (x | (y << 1));
    } else
    {
        /* 3D case */
        x &= 0x1fffff;
        x = (x | x << 32) & 0x1f00000000ffff;
        x = (x | x << 16) & 0x1f0000ff0000ff;
        x = (x | x << 8) & 0x100f00f00f00f00f;
        x = (x | x << 4) & 0x10c30c30c30c30c3;
        x = (x | x << 2) & 0x1249249249249249;

        y &= 0x1fffff;
        y = (y | y << 32) & 0x1f00000000ffff;
        y = (y | y << 16) & 0x1f0000ff0000ff;
        y = (y | y << 8) & 0x100f00f00f00f00f;
        y = (y | y << 4) & 0x10c30c30c30c30c3;
        y = (y | y << 2) & 0x1249249249249249;

        z &= 0x1fffff;
        z = (z | z << 32) & 0x1f00000000ffff;
        z = (z | z << 16) & 0x1f0000ff0000ff;
        z = (z | z << 8) & 0x100f00f00f00f00f;
        z = (z | z << 4) & 0x10c30c30c30c30c3;
        z = (z | z << 2) & 0x1249249249249249;

        /* Interleave the three integers */
        return (x | (y << 1) | (z << 2));
    }

}

/* Examine to which rank the supplied morton index belongs */
static int
get_rank_assignment_zcurve(const std::vector<uint64_t>& offsets, const uint64_t& morton_index)
{
    
    for (int rank_id{0}; rank_id < offsets.size() -1; ++rank_id)
    {
        if (offsets[rank_id +1] > morton_index)
        {
            return rank_id;
        } else
        {
            continue;
        }
    }

    return offsets.size() -1;
}


static
cmc_mpi_reading_data_distribution
cmc_mpi_calculate_reading_data_distribution_blocked(const MPI_Comm comm, const cmc_mpi_data_distribution_info& info)
{
    int err, rank, size;

    /* Get the size of the communicator */
    err = MPI_Comm_size(t8_data.comm, &comm_size);
    cmc_mpi_check_err(err);
    
    /* Get the (local) rank id */
    err = MPI_Comm_rank(t8_data.comm, &rank);
    cmc_mpi_check_err(err);

    /* Variables for calculating the distribution */
    uint32_t leftover_dim_points{0}, evenly_split_dim_points{0};
    uint32_t global_offset{0};
    /* Flag indicating that one diemnsion was already not split (in oder to cover the whole domain) */
    bool non_splitted_dim_applied{false};

    std::vector<uint32_t> start_values;
    std::vector<uint32_t> count_values;

    /* Iterate over all dimensions and determine a blocked distribution */
    for (auto iter{info.count_values.begin()}; iter != info.count_values.end(); ++iter)
    {
        /* If a dimension is not considered in the data hyperslab, just skip it */
        if (*iter <= 1)
        {
            start_values.push_back(0);
            count_values.puch_back(1);
            continue;
        } else
        {
            //Currently, there is only one option to calculate a blocked distribution, therefore there are some assertion 
            cmc_assert(*iter >= static_cast<size_t>(comm_size));
            if (non_splitted_dim_applied)
            {
                /* If the dimension is considered, calculate an offset for this dimension */
                //The assertion above indicates that at least one dimension point is assigned to each rank
                /* The leftover dimension points will be distributed */
                evenly_split_dim_points = static_cast<uint32_t>((*iter) / comm_size);
                /* Calculate the offset for the start vector */
                global_offset = rank * evenly_split_dim_points;
                /* Calculate the points which cannot be evenly split between all processses */
                leftover_dim_points = (*iter) % comm_size;
                /* Adjust the offset for the leftover points */
                /* Rank 0 starts at position zero nevertheless, therefore, the offset cannot change */
                if (rank != 0)
                {
                    if (rank < static_cast<int>(leftover_dim_points))
                    {
                        /* The ranks (starting from the lowest to the highest id) obtain another additional dimension point when the diemnsion cannot be split evenly */
                        ++(evenly_split_dim_points);
                        /* Adjust the global offset */
                        global_offset += leftover_dim_points - rank;
                    } else
                    {
                        /* Adjust the global offset */
                        global_offset += leftover_dim_points;
                    }
                }

                /* Save the start and count values for this dimension */
                start_values.push_back(info.start_values[dims] + global_offset);
                count_values-push_back(evenly_split_dim_points);
                nc_data->vars[current_var_id]->start_ptr[dims] += global_offset;
                nc_data->vars[current_var_id]->count_ptr[dims] = evenly_split_dim_points;
            }
            else
            {
                /* Just copy the whole dimension length */
                start_values.push_back(info.start_values[dims]);
                count_values-push_back(info.count_values[dims]);
                /* The whole dimensions are already saved by default in start and count pointers */
                non_splitted_dim_applied = true;
            }
        }
    }

    /* Return a struct dscribing the data which should be read in on each rank */
    return cmc_mpi_reading_data_distribution(std::move(start_values), std::move(count_values), DATA_DISTRIBUTION::CMC_BLOCKED);
}

cmc_mpi_reading_data_distribution
cmc_mpi_calculate_reading_data_distribution(const MPI_Comm comm, const cmc_mpi_data_distribution_info& info, const DATA_DISTRIBUTION preferred_data_distribution)
{
    int err, rank, size;

    /* Get the size of the communicator */
    err = MPI_Comm_size(t8_data.comm, &comm_size);
    cmc_mpi_check_err(err);
    
    /* Get the (local) rank id */
    err = MPI_Comm_rank(t8_data.comm, &rank);
    cmc_mpi_check_err(err);

    switch(preferred_data_distribution)
    {
        case CMC_SERIAL_DATA:

        break;
        case CMC_LINEARIZED:

        break;
        case CMC_BLOCKED:
            return cmc_mpi_calculate_reading_data_distribution_blocked(comm, info);
        break;
        case CMC_ZCURVE:

        break;
        default:
            cmc_err_msg("An unknown data distribution was supplied.");
    }

}


#if 0
struct cmc_mpi_data_distribution_info
{
    int initial_refinement_level{0}; //!< Initial uniform refinement level of a mesh wich embeds the geo-spatial variable data
    int data_dimensionality{0}; //!< The dimensionality of the concerned data (only 2D or 3D data is currently possible)
    uint32_t local_mesh_elements{0}; //!< The local number of mesh elements
    uint64_t global_mesh_elements{0}; //!< The global number of mesh elements
    std::array<uint64_t, 3> start_values{0,0,0}; //!< Global start values of the local data of a 2D or 3D variable
    std::array<uint64_t, 3> count_values{0,0,0}; //!< Local count values of the local data of a 2D or 3D variable
    std::vector<cmc_coordinates> offsets;
    t8_forest_t forest{nullptr}; //!< The forest which describes the morton/zcurve order
    uint64_t linear_zcurve_offset{0}; //!< Global linear offset of a uniform refinement (corresponding to the @var initial_refinement_level)
};

struct cmc_mpi_transfer_data
{
    int partner_rank{-1};
    size_t num_elements{0};
    std::vector<var_dynamic_array_t*> data;
    std::vector<cmc_coordinates> coordinates;
    DATA_DISTRIBUTION distribution{DISTRIBUTION_UNDEFINED};
};
#endif
cmc_mpi_transfer_data
cmc_mpi_calculate_transfer_data_blocked_to_zcurve(const MPI_Comm comm, const cmc_mpi_data_distribution_info& info)
{
    cmc_assert(info.data_dimensionality == 2 || info.data_dimensionality == 3);
    int err, rank, size;

    /* Get the size of the communicator */
    err = MPI_Comm_size(comm, &comm_size);
    cmc_mpi_check_err(err);
    
    /* Get the (local) rank id */
    err = MPI_Comm_rank(comm, &rank);
    cmc_mpi_check_err(err);

    //TODO: Eigentlich müssen die counts nicht nit kommuniziert werden
    /* We also do need all start and count ptrs in order to know from which rank we will receive data */
    std::array<uint32_t, 7> all_offsets;
    /* The first entry is the offset of the first local element */
    all_offsets[0] = elem_offset;
    /* Followed by start_pointers and lasty three count_pointers */
    /* Start pointers for the first and second dimension */
    all_offsets[1] = static_cast<uint32_t>(info.start_ptr[0]);
    all_offsets[2] = static_cast<uint32_t>(info.start_ptr[1]);
    /* If 3D data is considered */
    all_offsets[3] = (info.data_dimensionality == 2) ? 0 : static_cast<uint32_t>(info.start_ptr[2]);
    /* Count pointers for the first and second dimension */
    all_offsets[4] = static_cast<uint32_t>(info.count_ptr[0]);
    all_offsets[5] = static_cast<uint32_t>(info.count_ptr[1]);
    /* If 3D data is considered */
    all_offsets[6] = (info.data_dimensionality == 2) ? 0 : static_cast<uint32_t>(info.count_ptr[2]);
    /* Disribute the offsets between all ranks */
    uint64_t* offsets = (uint64_t*) malloc(sizeof(uint64_t) * 7 * comm_size);
    err = MPI_Allgather(all_offsets.data(), 7, MPI_UINT32_T, offsets, 7, MPI_UINT32_T, comm);
    cmc_mpi_check_err(err);

    /* Seperate the element offsets */
    std::vector<uint64_t> elem_offsets;
    elem_offsets.reserve(comm_size);
    for (int i{0}; i < comm_size; ++i)
    {
        elem_offsets.push_back(static_cast<uint64_t>(offsets[i * 7]));
    }

    //delete later
    if (rank == 0)
    {
        std::cout << "offsets in rank 0 are:\n";
        for (int i{0}; i < comm_size; ++i)
        {
            std::cout << elem_offsets[i] << ", ";
        }
    }

    /** Calculate the data which will be send **/
    /* Vector storing all 'send information' */
    std::vector<cmc_t8_mpi_send_data> send_list;
    send_list.reserve(comm_size);
    uint64_t morton_index{0};
    int recv_rank{-1};

    /* A size_hint for allocation for members of 'cmc_t8_mpi_send_data' */
    const size_t size_hint{global_num_elems / comm_size +1};

    /* This variable resembles the linear storage id of the current data point */
    size_t linear_count_of_data_entries{0};
    /* This flag indicates whether the rank which will receive the data is already in the 'send_list' */
    bool rank_already_receives_data{false};
    int receiving_rank{-1};

    /**************************************************/
    /** Examine to which processes data will be sent **/
    std::cout << "rank " << rank << " has " << info.count_ptr[0] * info.count_ptr[1] << "data points" << std::endl;
    /* Check if we got a 2D or 3D variable */
    if (info.data_dimensionality == 2)
    {
        /* If the variable is 2D */
        /* Iterate over the coordinates the local data is defined on and determine to which rank the data needs to be transferred */
        //TODO: currently, this depends on the start and count ptr of the first variable, this works in case all vairables are defined on the same domain, but this information about start and count should be stored in 'general place'
        std::cout << "Rank " << rank << ", hat start_ptr: " << info.start_ptr[0] << " und " << info.start_ptr[1] <<  "; und count_ptr: " << info.count_ptr[0] << " und " << info.count_ptr[1] <<  std::endl;
        for (uint64_t dim1{info.start_ptr[0]}; dim1 < info.start_ptr[0] + info.count_ptr[0]; ++dim1)
        {
            for (uint64_t dim2{info.start_ptr[1]}; dim2 < info.start_ptr[1] + info.count_ptr[1]; ++dim2)
            {
                for (uint64_t dim3{info.start_ptr[2]}; dim3 < info.start_ptr[2] + info.count_ptr[2]; ++dim3)
                {
                    /* Calculate the morton index */
                    morton_index = cmc_linear_id_to_zcurve_id(info.data_dimensionality, dim1, dim2, dim3);
                    /* Check to which rank the element has to be transferred to */
                    recv_rank = get_rank_assignment_zcurve(elem_offsets, morton_index);
                    /* Check if the receiving rank already is in the 'send list'. If so, assign the new element, oterwise a new receiver will be added */
                    for (size_t receiver{0}; receiver < send_list.size(); ++receiver)
                    {
                        if (send_list[receiver].partner_rank == recv_rank)
                        {
                            receiving_rank = receiver;
                            rank_already_receives_data = true;
                        }
                    }
                    /* If the rank does not yet receive data, add it to the list */
                    if (rank_already_receives_data == false)
                    {
                        /* Set the correct rank id */
                        receiving_rank = send_list.size();
                        /* If the receiving rank is not yet in the 'send list', it will be added */
                        send_list.emplace_back(cmc_t8_mpi_send_data(recv_rank, t8_data.vars.size()));
                        /* Reserve memory for the morton indices */
                        send_list[receiving_rank].morton_indices.reserve(size_hint);
                        /* Allocate the dynamic arrays holding the data */
                        for (size_t var_id{0}; var_id < t8_data.vars.size(); ++var_id)
                        {
                            send_list[receiving_rank].data.push_back(new var_dynamic_array_t(size_hint, t8_data.vars[var_id]->get_type()));
                        }
                    }

                    /* Add the morton index */
                    send_list[receiving_rank].morton_indices.push_back(morton_index);
                    /* Add the data of each variable */
                    for (size_t var_id{0}; var_id < t8_data.vars.size(); ++var_id)
                    {
                        send_list[receiving_rank].data[var_id]->push_back(t8_data.vars[var_id]->var->data->operator[](linear_count_of_data_entries));
                    }

                    /* Reset the flag */
                    rank_already_receives_data = false;
                    /* Increment the pointer */
                    ++linear_count_of_data_entries;
                }
            }
        }
    }

    //delete later
    #if 1
    if (rank >= 0)
    {
        /* The sender list has been filled */
        std::cout << "Size of send list: " << send_list.size() << std::endl;
        std::cout << "Rank " << rank << " hat folgende Send Liste:" << "\n";
        for (size_t s{0}; s < send_list.size(); ++s)
        {
            std::cout << "Receiving rank: " << send_list[s].partner_rank << " -> bekommt " << send_list[s].data[0]->size() << " Elemente.";
            std::cout << "\n";
            std::cout << "Size morton indices: " << send_list[s].morton_indices.size() << "\n";
        }
        std::cout <<std::endl;
    }
    #endif


    /********************************************************/
    /** Examine from which processes data will be received **/
    /* List of receiving data */
    std::vector<cmc_t8_mpi_recv_data> recv_list;
    recv_list.reserve(comm_size);

    /* Cartesian coordinate */
    std::tuple<uint64_t,uint64_t,uint64_t> linear_id;

    /* Variables for a mesh element and its coordinates */
    t8_element_t *elem;
    //int x_coord{0}, y_coord{0}, z_coord{0};
    //int vertex_coords[3];

    /* offset variable for accesssing the all_offsets vector */
    size_t offset_geo_domains = 0;

    /* Flag indicating whether the a new sender will be added to the receivng list */
    bool sender_already_in_recv_list{false};


    /* Check from which rank data will be received */
    for (size_t elem_id{0}; elem_id < local_num_elems; ++elem_id)
    {
        /* Get the element to the corresponding element id */
        elem = t8_forest_get_element_in_tree(forest, 0, elem_id);

        /* Check which rank needs this data */
        for (int rank_id{0}; rank_id < comm_size; ++rank_id)
        {
            /* Skip the entries of all fromer ranks and the local element offset */
            offset_geo_domains = rank_id * 7 + 1;
            if (cmc_t8_elem_in_storage_domain(info.dim, info.initial_refinement_lvl, elem, scheme_eclass,
                                              all_offsets[offset_geo_domains], all_offsets[offset_geo_domains + 1], all_offsets[offset_geo_domains + 2],
                                              all_offsets[offset_geo_domains] + all_offsets[offset_geo_domains + 3],
                                              all_offsets[offset_geo_domains + 1] + all_offsets[offset_geo_domains + 4],
                                              all_offsets[offset_geo_domains + 2] + all_offsets[offset_geo_domains + 5]) != 0)
            {
               
                /* If the element is inside of the domain of the rank 'rank_id', we need to receive this data */
                /* Check if the receiving rank already is in the 'send list'. If so, assign the new element, oterwise a new receiver will be added */
                for (size_t sender{0}; sender < recv_list.size(); ++sender)
                {
                    if (recv_list[sender].partner_rank == rank_id)
                    {
                        sender_already_in_recv_list = true;
                        ++(recv_list[sender].num_receiving_elems);
                        break;
                    }
                }
                /* If the rank is not already in the receiving_list, we add it and the increment the count for this rank */
                if (sender_already_in_recv_list == false)
                {
                    recv_list.emplace_back(cmc_t8_mpi_recv_data(rank_id, t8_data.vars.size()));
                    ++(recv_list.back().num_receiving_elems);
                }
                break;
            }
        }

        /* Reset the flag */
        sender_already_in_recv_list = false;
    }

}

cmc_mpi_calculate_transfer_data(const MPI_Comm comm, const DATA_DISTRIBUTION from_distribution, const DATA_DISTRIBUTION to_distribution)
{
    switch(to_distribution)
    {
        case CMC_SERIAL_DATA:

        break;
        case CMC_LINEARIZED:

        break;
        case CMC_BLOCKED:
            switch (from_distribution)
            {
                case CMC_BLOCKED:
                    return cmc_mpi_calculate_transfer_data_blocked_to_zcurve();
                break;
                default:
                    cmc_err_msg("This data redistribution is not implemented yet.");
            }
        break;
        case CMC_ZCURVE:

        break;
        default:
            cmc_err_msg("An unknown data distribution was supplied.");
    }
}