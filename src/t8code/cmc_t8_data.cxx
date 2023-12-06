#include "t8code/cmc_t8_data.hxx"
#include <limits.h>

namespace cmc 
{

inline static
int CalculateInitialRefinementLevel()
{

}

static
int DetermineDimensionalityOfTheData()
{

}

inline static
t8_eclass_t DimensionToElementClass(const int dimensionality)
{
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    return (dimensionality == 3 ? T8_ECLASS_HEX : T8_ECLASS_QUAD);
}

inline static
void ValidateInitialRefinementLevelForDimensionality(const int initial_refinement_level, const int dimensionality)
{
    #ifdef CMC_ENABLE_DEBUG
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    
    const int maximum_possible_refinement_level = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);
    
    if (initial_refinement_level < 0 ||
        initial_refinement_level > maximum_possible_refinement_level)
    {
        cmc_err_msg("The corresponding refinement level is not within the range of an computationally posiible refinement level.");
    }

    #endif
}

void
AmrData::BuildInitialEmbeddedMesh()
{
    const int dimensionality = DetermineDimensionalityOfTheData();

    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(DimensionToElementClass(dimensionality), comm, 0, 0, 1);

    const int initial_refinement_level = CalculateInitialRefinementLevel();

    initial_mesh.SetInitialRefinementLevel(initial_refinement_level);

    ValidateInitialRefinementLevelForDimensionality(initial_refinement_level, dimensionality);
    
    t8_forest_t initial_forest = t8_forest_new_uniform(cmesh, t8_scheme_new_default_cxx(), initial_refinement_lvl, 0, comm);

    initial_forest = CoarsenInitialEmebddedMeshToGlobalDomain(initial_forest, );

    initial_mesh.SetMesh(initial_forest);
}

void
AmrData::DistributeDataOnInitialMesh()
{
    cmc_assert(initial_mesh.IsValid());

    //if MPI
    GatherGlobalDataOffsets();

    ComputeMortonIndices();

    // if MPI
    CommunicateNonCompliantData();

    SortReceivedDataLocally();
}


inline AdaptData
AmrData::CreateAdaptationData(const AmrMesh& adaptation_sample) const
{
    return AdaptData(compression_settings, adaptation_sample.variable_id, variables);
}

inline std::vector<CoarseningSample>
AmrData::RetrieveMeshesToBeCoarsened(const CompressionMode compression_mode) const
{
    cmc_assert(!variables->empty());

    std::vector<CoarseningSample> mesh_samples_to_be_coarsened;

    switch (compression_mode)
    {
        case CompressionMode::OneForAll:
            /* In case of an One For All compression, all variables are defined and coarsend on the same mesh */
            mesh_samples_to_be_coarsened.emplace_back(CoarseningSample(kMeshCorrespondsToAllVariables));
        break;
        case CompressionMode::OneForOne:
            /* In case of a One For One compression, every variable is coarsened on it's own mesh */
            mesh_samples_to_be_coarsened.reserve(variables->size());
            for (auto iter = variables->begin(); iter != variables->end(); ++iter)
            {
                /* We push back the 'id' of each variable (i.e. 0, 1, 2,...,#num_variables - 1), since each variable has it's own mesh which needs to be coarsened */
                mesh_samples_to_be_coarsened.emplace_back(CoarseningSample(std::distance(variables->begin(), iter)));
            }
        break;
        default:
            cmc_err_msg("The compression mode has not been specified.");
    }

    return mesh_samples_to_be_coarsened;
}

void
AmrData::CompressByAdaptiveCoarsening(const CompressionMode compression_mode)
{
    #ifdef CMC_WITH_T8CODE
    //If Debug
    //Save initial data

    /* Retrieve all forests which are going to be coarsened. In case of a 'One For All'-compression we will just
     * receive a single forest on which all variables are defined. However, in case of a 'One for One'-compression,
     * we will receive multiple forests which need to be coarsened, in particular a single forest for each variable */
    std::vector<CoarseningSample> meshes_to_be_coarsened = AmrData.RetrieveMeshesToBeCoarsened(compression_mode);

    for (auto current_sample = meshes_to_be_coarsened.begin(); current_sample != meshes_to_be_coarsened.end(); ++current_sample)
    {
        /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
        AdaptData adapt_data = CreateAdaptationData(*current_sample);

        /* Get the adaptation function for the coarsening step (which is based on the compression settings) */
        t8_forest_adapt_t adaptive_coarsening_function = adapt_data.GetAdaptationFunction();

        /* Iterate until none further coarsening is possible */
        while (adapt_data.IsCompressionProgressing())
        {
            /* Get the forest before the coarsening step */
            t8_forest_t previous_forest = adapt_data.GetCurrentMesh();

            /* Keep the 'previous forest' after the adaptation step */
            t8_forest_ref(previous_forest);

            /* Perform a coarsening iteration */
            t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, adaptive_coarsening_function, 0, 0, static_cast<void*>(&adapt_data));

            /* Interpolate the data from the previous forest to the (new) coarsened forest */
            adapt_data.InterpolateData(adapted_forest);
        
            /* Repartition the mesh as well as the data */
            adapted_forest = adapt_data.RepartitionData(adapted_forest);

            /* Free the former forest */
            t8_forest_unref(&previous_forest);

            /* Store the (new) coarsened forest */
            adapt_data.SetCurrentMesh(adapted_forest);
        }
    }

    #endif
}

}

