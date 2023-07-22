#include "cmc_t8code_data.hxx"
#include "mpi/cmc_mpi.h"
#include "utilities/cmc_geo_util.h"
#include "utilities/cmc_container.h"
#include "utilities/cmc_log_functions.h"

cmc_type 
cmc_t8_var::get_type() const
{
    cmc_assert(var != nullptr && var->data != nullptr);
    return var->data->get_data_type();
}

size_t 
cmc_t8_var::get_data_size() const
{
    cmc_assert(var != nullptr && var->data != nullptr);
    return cmc_type_to_bytes[var->data->get_data_type()];
}

int
cmc_t8_var::get_dimensionality() const
{
    return var->num_dimensions;
}

void*
cmc_t8_var::get_initial_data_ptr() const
{
    cmc_assert(var != nullptr && var->data != nullptr);
    return var->data->get_initial_data_ptr();
}

void*
cmc_t8_var::get_initial_data_new_ptr() const
{
    cmc_assert(var != nullptr && var->data != nullptr);
    return var->data_new->get_initial_data_ptr();
}

DATA_LAYOUT
cmc_t8_var::get_data_layout() const
{
    cmc_assert(var != nullptr);
    return var->data_layout;
}

size_t
cmc_t8_geo_data::get_coord_length(const CMC_COORD_IDS cmc_coord_id) const
{
    cmc_assert(coordinates->coords.size() >= CMC_NUM_COORD_IDS);
    return (coordinates->coords.operator[](cmc_coord_id)).size();
}

/* Write the forest as well as the variables in a vtk-file (depending on the compression mode (if supplied) one file with all variables or one file per variable is written) */
//TODO: This function needs to be rewritten in a bit more efficient way 
void
cmc_t8_write_forest_all_vars(cmc_t8_data_t t8_data, const char* file_prefix)
{
    #ifdef CMC_WITH_T8CODE
    /* Arrays holding the ids of variables of the given data ordering */
    std::vector<int> var_ids_zcurve;
    var_ids_zcurve.reserve(t8_data->vars.size());

    /* Count the variables which cannot be written in a vtk file */
    size_t num_vars_undefined_data_scheme{0};

    /* Iterate over all variables and check if their data_scheme is compliant */
    for (size_t id{0}; id < t8_data->vars.size(); ++id)
    {
        if (t8_data->vars[id]->var->data_scheme == CMC_DATA_ORDERING_SCHEME::CMC_GEO_DATA_Z_CURVE)
        {
            var_ids_zcurve.push_back(id);
        }
        else
        {
            ++num_vars_undefined_data_scheme;
        }
    }

    /* If all data schemes are undefined, just write out the forest without any additional data */
    if (num_vars_undefined_data_scheme == t8_data->vars.size())
    {
        /* Check if an initial forest is supplied in t8_data */
        if (t8_data->initial_forest == nullptr)
        {
            cmc_warn_msg("There is no initial_forest in t8_data and there are none variables with a compliant ordering scheme.");
            cmc_warn_msg("No vtk-file was written.");
            return;
        }

        /* Write out the forest in a vtk file */
        int vtk_err{t8_forest_vtk_write_file (t8_data->initial_forest, file_prefix, 1, 1, 1, 1, 0, 0, NULL)};
        if (vtk_err == 0)
        {
            cmc_err_msg("An error occrued during the creation of the t8code-forest vtk-file.");
        }
    }
    /* In case there is some z-curve ordered data is present, write out the data with the forest itself */ 
    else
    {
        /* If all variables are defined on the same forest with an 'One For All'-compression mode */
        if (t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
        {
            /* Get the amount of variables to write out */
            const size_t num_data_vars{var_ids_zcurve.size()};
            
            /* Flag indicating if non-double variables have to be transformed */
            std::vector<bool> calculate_data;
            calculate_data.reserve(num_data_vars);

            /* Create a new vtk field holding the element data arrays */
            t8_vtk_data_field_t *vtk_data = new t8_vtk_data_field_t[num_data_vars];

            /* Create an data array holding the element values */
            double **data_array = T8_ALLOC(double*, num_data_vars);
            
            /* Get the number of elements of the forest */
            const t8_locidx_t num_elems{t8_forest_get_local_num_elements(t8_data->assets->forest)};

            /* Standard offset and scaling_factor */
            double scale_factor{1.0};
            double add_offset{0.0};

            int vtk_var_id{0};

            /* Write out variables in zcurve order */
            for (auto iter{var_ids_zcurve.begin()}; iter != var_ids_zcurve.end(); ++iter, ++vtk_var_id)
            {
                /* Get scaling and offset attributes from the varibale, in order to check if the data has to be scaled before it is written */
                scale_factor = cmc_get_universal_data<double>(t8_data->vars[*iter]->var->scale_factor);
                add_offset = cmc_get_universal_data<double>(t8_data->vars[*iter]->var->add_offset);

                /* Checl if the data for this variable has to be calculated first */
                calculate_data.push_back(!((t8_data->vars[*iter]->get_type() == cmc_type::CMC_DOUBLE) && (t8_data->vars[*iter]->var->applied_offset_scaling || (cmc_approx_cmp(add_offset, 0.0) && cmc_approx_cmp(scale_factor, 1.0)))));
                if (calculate_data.back())
                {
                    /* If the variable is not of type double */
                    /* Allocate an array for the current variable's element data */ 
                    data_array[vtk_var_id] = T8_ALLOC(double, num_elems);

                    /* Check if the meta data of the variable holds offset and scaling attributes and if so, if they already have been applied */
                    if (t8_data->vars[*iter]->var->applied_offset_scaling || (cmc_approx_cmp(add_offset, 0.0) && cmc_approx_cmp(scale_factor, 1.0)))
                    {
                        /* Assign data WITHOUT applying value offsets and scalings */
                        for (int elem_id{0}; elem_id < num_elems; ++elem_id)
                        {
                            data_array[vtk_var_id][elem_id] = cmc_get_universal_data<double>(t8_data->vars[*iter]->var->data->operator[](elem_id));
                        }
                    } else
                    {
                        /* Assign data WITH applying value offsets and scalings */        
                        for (int elem_id{0}; elem_id < num_elems; ++elem_id)
                        {
                            data_array[vtk_var_id][elem_id] = add_offset + scale_factor * cmc_get_universal_data<double>(t8_data->vars[*iter]->var->data->operator[](elem_id));
                        }
                    }
                    /* Set the type of the data and pointer to the data */
                    snprintf(vtk_data[vtk_var_id].description, 50, "%s", (t8_data->vars[*iter]->var->name).c_str());
                    vtk_data[vtk_var_id].type = T8_VTK_SCALAR;
                    vtk_data[vtk_var_id].data = static_cast<double*>(data_array[vtk_var_id]);
                } else
                {
                    /* If the variable is already in z-curve order and of type double */
                    /* Set the type of the data and pointer to the data */
                    snprintf(vtk_data[vtk_var_id].description, 50, "%s", (t8_data->vars[*iter]->var->name).c_str());
                    vtk_data[vtk_var_id].type = T8_VTK_SCALAR;
                    vtk_data[vtk_var_id].data = static_cast<double*>(t8_data->vars[*iter]->get_initial_data_ptr());

                }
            }

            /* Write out the vtk file of the forest with the additional variable data */
            int vtk_err{t8_forest_vtk_write_file(t8_data->assets->forest, file_prefix, 1, 1, 1, 1, 0, num_data_vars, vtk_data)};
            if (vtk_err == 0)
            {
                cmc_err_msg("An error occrued during the creation of the t8code-forest vtk file.");
            }

            /* Free allocated memory */
            delete[] vtk_data;
            for(size_t i{0}; i < num_data_vars; ++i)
            {
                if (calculate_data[i])
                {
                    T8_FREE(data_array[i]);
                }
            }
            T8_FREE(data_array);
            cmc_debug_msg("A VTK-File of all variables has been written (Compression-Mode: 'One for All').");
        }
        /* If each variable is defined on its own forest */
        else {
            /* Number of elements within the forest */
            t8_locidx_t num_elems{0};

            /* String holding the fileprefix */
            char _file_prefix[100];

            /* VTK Error check */
            int vtk_err{0};

            /* Create a new vtk field holding the element data arrays */
            t8_vtk_data_field_t *vtk_data = new t8_vtk_data_field_t[1];

            /* Create an data array holding the element values */
            double* data_array{nullptr};

            /* Standard offset and scaling_factor */
            double scale_factor{1.0};
            double add_offset{0.0};

            /* A flag indicating whther the data of the variable has to be calculated first */
            bool calculate_data{true};

            /* Write out variables in zcurve order */
            for (auto iter{var_ids_zcurve.begin()}; iter != var_ids_zcurve.end(); ++iter)
            {
                /* Get the offset and scaling attribute from the varibale */
                scale_factor = cmc_get_universal_data<double>(t8_data->vars[*iter]->var->scale_factor);
                add_offset = cmc_get_universal_data<double>(t8_data->vars[*iter]->var->add_offset);

                /** If it is a double variable and no offset and scaling needs to be applied, we can directly copy the data from the variable.
                 *  Otherwise, we have to allocate a new array ad save the data as double data and/or apply the scaling and offset to the data
                **/
                calculate_data = !((t8_data->vars[*iter]->get_type() == cmc_type::CMC_DOUBLE) && (t8_data->vars[*iter]->var->applied_offset_scaling || (cmc_approx_cmp(add_offset, 0.0) && cmc_approx_cmp(scale_factor, 1.0))));
                if (calculate_data)
                {
                    /* Get the number of forest elements */
                    num_elems = t8_forest_get_local_num_elements(t8_data->vars[*iter]->assets->forest);

                    /* Allocate an array for the current variable's element data */ 
                    data_array = T8_ALLOC(double, num_elems);
                    
                    /* Check if the meta data of the variable holds offset and scaling attributes and if so, if they already have been applied */
                    if (t8_data->vars[*iter]->var->applied_offset_scaling || (cmc_approx_cmp(add_offset, 0.0) && cmc_approx_cmp(scale_factor, 1.0)))
                    {
                        /* Assign data WITHOUT applying value offsets and scalings */
                        for (int elem_id{0}; elem_id < num_elems; ++elem_id)
                        {
                            data_array[elem_id] = cmc_get_universal_data<double>(t8_data->vars[*iter]->var->data->operator[](elem_id));
                        }
                    } else
                    {
                        /* Assign data WITH applying value offsets and scalings */        
                        for (int elem_id{0}; elem_id < num_elems; ++elem_id)
                        {
                            data_array[elem_id] = add_offset + scale_factor * cmc_get_universal_data<double>(t8_data->vars[*iter]->var->data->operator[](elem_id));
                        }
                    }

                    /* Set the type of the data and pointer to the data */
                    snprintf(vtk_data[0].description, 50, "%s", (t8_data->vars[*iter]->var->name).c_str());
                    vtk_data[0].type = T8_VTK_SCALAR;
                    vtk_data[0].data = static_cast<double*>(data_array);

                } else
                {
                    cmc_assert(t8_data->vars[*iter]->var->data->size() == static_cast<size_t>(t8_forest_get_local_num_elements(t8_data->vars[*iter]->assets->forest)));
                    /* Set the type of the data and pointer to the data */
                    snprintf(vtk_data[0].description, 50, "%s", (t8_data->vars[*iter]->var->name).c_str());
                    //Quick fix for visualization
                    //snprintf(vtk_data[0].description, 3, "%s", "O3");
                    vtk_data[0].type = T8_VTK_SCALAR;
                    vtk_data[0].data = static_cast<double*>(t8_data->vars[*iter]->var->data->get_initial_data_ptr());
                    #if 0
                    std::cout << "VarID: " << *iter << ", Forest Num elements: " << t8_forest_get_global_num_elements(t8_data->vars[*iter]->assets->forest) << ", Num data array: " << t8_data->vars[*iter]->var->data->size() << std::endl;
                    std::cout << "Type: "  << t8_data->vars[*iter]->get_type() << std::endl;
                    if (*iter == 0)
                    {
                        for (int i{0}; i < t8_forest_get_global_num_elements(t8_data->vars[*iter]->assets->forest); ++i)
                        {
                            std::cout << (vtk_data[0].data)[i];
                        }
                        std::cout << std::endl;
                    }
                    #endif
                }

                /* Assign a file prefix */
                vtk_err = snprintf(_file_prefix, strlen(file_prefix) + 2 + strlen((t8_data->vars[*iter]->var->name).c_str()), "%s_%s", file_prefix, (t8_data->vars[*iter]->var->name).c_str());
                //Quick fix for visualization
                //vtk_err = snprintf(_file_prefix, strlen(file_prefix) + 10 + strlen((t8_data->vars[*iter]->var->name).c_str()), "%s_%s%d", file_prefix,  "_O3_", *iter);
                
                /* Write out the variable data with the forest on which its data is defined */
                vtk_err = t8_forest_vtk_write_file(t8_data->vars[*iter]->assets->forest, _file_prefix, 1, 1, 1, 1, 0, 1, vtk_data);
                if (vtk_err == 0)
                {
                    cmc_err_msg("An error occrued during the creation of the t8code-forest vtk file.");
                }
                
                /* Free allocated memory if it was allocated */
                if (calculate_data)
                {
                    T8_FREE(data_array);
                }
            }

            /* Free the vtk data field */
            delete[] vtk_data;

            cmc_debug_msg("A VTK-File of all variables has been written (Compression-Mode: 'One for One').");
        }
    }
    #endif
}
