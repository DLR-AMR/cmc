#include "lossy/cmc_prefix_lossy_compression.hxx"
#include "t8code/cmc_t8_adapt_callbacks.h"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_span.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#endif

#include <vector>

namespace cmc
{

int
PrefixCompressionData::GetMpiSize() const
{
    int comm_size{1};
    int err = MPI_Comm_size(comm_, &comm_size);
    MPICheckError(err);
    return comm_size;
}

void
PrefixCompressionData::Setup(const bool with_default_lossy_amr_compression)
{
    compression_data_->SplitVariables();

    compression_data_->CheckConsistencyOfInputVariables();

    compression_data_->BuildInitialMesh();

    compression_data_->DistributeDataOnInitialMesh();

    compression_data_->ApplyScalingAndOffset();

    compression_data_->SetupVariablesForCompression();

    if (!with_default_lossy_amr_compression)
    {
        /* In case a default lossy compression needs to be perfomed first, we have to postpone the transfomration to byte variables */
        compression_variables_ = compression_data_->GetByteVariablesForCompression();
    } else
    {
        perform_default_lossy_compression_ = true;
    }
}

void
PrefixCompressionData::Compress()
{
    /* Check if the default lossy compression is about to be applied */
    if (perform_default_lossy_compression_)
    {
        /* Perform the default lossy compression */
        compression_data_->CompressByAdaptiveCoarsening(CompressionMode::OneForOne);

        /* Transform the data to byte variables */
        compression_variables_ = compression_data_->GetByteVariablesForCompression();

        /* Store the initial compressed mesh */
        for (auto cv_iter = compression_variables_.begin(); cv_iter != compression_variables_.end(); ++cv_iter)
        {
            cv_iter->StoreInitialMesh();
        }
    }

    /* Perform trail truncation until the error trhesholds are exhausted */
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        var_iter->PerformTailTruncation();
    }

    /* Afterwards, we create prefixes in the tree hierachy */
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
        PrefixAdaptDataEGU adapt_data{*var_iter};

        while (adapt_data.IsCompressionProgressing())
        {
            adapt_data.InitializeCompressionIteration();

            t8_forest_t previous_forest = adapt_data.GetCurrentMesh();

            /* Perform a coarsening iteration */
            t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, FindPrefixBits, 0, 0, static_cast<void*>(&adapt_data));

            adapt_data.SetCurrentMesh(adapted_forest);

            adapt_data.FinalizeCompressionIteration();
        }
    }
}

//New version for now
template<int N>
void
PrefixCompressionData::DecompressVariableEGU(std::vector<PrefixDecompressionData>& decompression_data) const
{
    /* Build the initial mesh, corresponding to the data */
    for (auto var_iter = decompression_data.begin(); var_iter != decompression_data.end(); ++var_iter)
    {
        /* Get a view on the prefix indication bytes */
        const int prefix_indication_byte_offset = 0;
        VectorView<uint8_t> prefix_indications(var_iter->serialized_variable.data() + prefix_indication_byte_offset, var_iter->num_bytes_prefix_indications);
        
        /* Get a view on the prefix lengths */
        const int prefix_lengths_byte_offset = prefix_indication_byte_offset + var_iter->num_bytes_prefix_indications;
        VectorView<uint8_t> prefix_lengths(var_iter->serialized_variable.data() + prefix_lengths_byte_offset, var_iter->num_bytes_prefix_lengths);

        /* Get a view on the prefix encodings */
        const int prefix_encodings_byte_offset = prefix_lengths_byte_offset + var_iter->num_bytes_prefix_lengths;
        VectorView<uint8_t> prefix_encodings(var_iter->serialized_variable.data() + prefix_encodings_byte_offset, var_iter->num_bytes_prefix_encodings);

        /* Get a view on the runlength encoding of the suffix indicator */
        const int rl_suffix_lengths_byte_offset = prefix_encodings_byte_offset + var_iter->num_bytes_prefix_encodings;
        VectorView<uint8_t> rl_suffix_indication(var_iter->serialized_variable.data() + rl_suffix_lengths_byte_offset, var_iter->num_bytes_run_length_suffix_indicator);

        /* Get a view on the suffix lengths */
        const int suffix_lengths_byte_offset = rl_suffix_lengths_byte_offset + var_iter->num_bytes_run_length_suffix_indicator;
        VectorView<uint8_t> suffix_lengths(var_iter->serialized_variable.data() + suffix_lengths_byte_offset, var_iter->num_bytes_suffix_lengths);

        /* Get a view on the suffix encodings */
        const int suffix_encodings_byte_offset = suffix_lengths_byte_offset + var_iter->num_bytes_suffix_lengths;
        VectorView<uint8_t> suffix_encodings(var_iter->serialized_variable.data() + suffix_encodings_byte_offset, var_iter->num_bytes_suffix_encodings);
        
        cmc_debug_msg("All vector views have been generated");

        /* Reconstruct the base mesh */
        t8_forest_t mesh = ReconstructBaseMesh(var_iter->domain.GetDimensionality(), comm_);

        cmc_debug_msg("Initial mesh was reconstructed with num elements: ", t8_forest_get_local_num_elements(mesh));

        DecompressionPrefixAdaptDataEGU<N> adapt_data(prefix_indications, prefix_lengths, prefix_encodings);
        
        cmc_debug_msg("Size prefix_indications: ", prefix_indications.size(), " size prefix_lengths: ", prefix_lengths.size(), " size prefix_encodings: ", prefix_encodings.size());
        adapt_data.domain = var_iter->domain;
        adapt_data.current_layout = var_iter->current_layout;

        cmc_debug_msg("Prefix decompression started");
        
        /* Decode and apply the prefixes to the data */
        while (adapt_data.IsDeCompressionProgressing())
        {
            cmc_debug_msg("Iteration: ", adapt_data.count_adaptation_step_);
            adapt_data.InitializeCompressionIteration(t8_forest_get_local_num_elements(mesh));

            mesh = t8_forest_new_adapt(mesh, DecodePrefixEGU<N>, 0, 0, static_cast<void*>(&adapt_data));

            adapt_data.FinalizeCompressionIteration();

        }
        cmc_debug_msg("Prefix decompression finished");
        
        /* Decode and apply the suffixes to the data */
        //TODO: Test direct indication bit
        //const auto [suffix_indication_bits, si_num_bits] = DecodeRunLengthEncoding(rl_suffix_indication);
        //VectorView<uint8_t> suffix_indication(suffix_indication_bits.data(), suffix_indication_bits.size());
        //cmc_debug_msg("Size suffix indication: ", suffix_indication_bits.size(), " und num_bits: ", si_num_bits);
        //cmc_debug_msg("Size of reconstructed variable before suffix addition: ", adapt_data.reconstructed_variable.size());
        //DecompressionPrefixAdaptDataEGU<N> suffix_adapt_data(suffix_indication, suffix_lengths, suffix_encodings);

        cmc_debug_msg("Size suffix indication: ", rl_suffix_indication.size(), " und num_bits: ", rl_suffix_indication.size() * 8);
        cmc_debug_msg("Size of reconstructed variable before suffix addition: ", adapt_data.reconstructed_variable.size());
        DecompressionPrefixAdaptDataEGU<N> suffix_adapt_data(rl_suffix_indication, suffix_lengths, suffix_encodings);

        suffix_adapt_data.reconstructed_variable = adapt_data.reconstructed_variable;
        suffix_adapt_data.count_adaptation_step_ = 1;
        suffix_adapt_data.domain = var_iter->domain;
        suffix_adapt_data.current_layout = var_iter->current_layout;

        cmc_debug_msg("Suffix decompression starts");

        /* Decode and apply the prefixes to the data */
           suffix_adapt_data.InitializeCompressionIteration(t8_forest_get_local_num_elements(mesh));

            mesh = t8_forest_new_adapt(mesh, DecodeSuffixEGU<N>, 0, 0, static_cast<void*>(&suffix_adapt_data));

            suffix_adapt_data.FinalizeCompressionIteration();

        cmc_debug_msg("Suffix decompression finished");

        cmc_debug_msg("Size of comrpession values: ", suffix_adapt_data.reconstructed_variable.size());


        //TODO: remove below
        #if 1
        cmc_debug_msg("Get data as float:");
        std::vector<float> decompressed_data_new = GetDataAsType<4>(suffix_adapt_data.reconstructed_variable);

        cmc_debug_msg("Remove this file IO");
        FILE* file2 = fopen("compr_initial_data.bin", "rb");
        std::vector<float> initial_data2(decompressed_data_new.size(), static_cast<float>(0.0));

        fread(initial_data2.data(), sizeof(float), decompressed_data_new.size(), file2);
        fclose(file2);

        float max_err2 = 0.0;
        for (int ii = 0; ii < decompressed_data_new.size(); ++ii)
        {
            //cmc_debug_msg("Decompressed Data: ", decompressed_data_new[ii], ", initial data: ", initial_data2[ii], " und Error: ", std::abs(decompressed_data_new[ii] - initial_data2[ii]));
            uint32_t inval{0};
            std::memcpy(&inval, &initial_data2[ii], 4);
            uint32_t deval{0};
            std::memcpy(&deval, &decompressed_data_new[ii], 4);
            cmc_debug_msg("Initial Value: ", std::bitset<8*4>(inval), " als value: ", initial_data2[ii]);
            cmc_debug_msg("Decompr Value: ",std::bitset<8*4>(deval), " als value: ", decompressed_data_new[ii]);
            cmc_debug_msg("Error Decompression: ", std::abs(decompressed_data_new[ii] - initial_data2[ii]));
            if (std::abs(decompressed_data_new[ii] - initial_data2[ii]) > max_err2)
            {
                max_err2 = std::abs(decompressed_data_new[ii] - initial_data2[ii]);
            }
        }


        cmc_debug_msg("Max absolute error is: ", max_err2);
        #endif
    }
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

void
PrefixCompressionData::WriteCompressedDataEGU(const std::string& file_name) const
{
    /* Create the netCDF file to which the cmompressed data will be written */
    cmc_assert(GetMpiSize() <= 1);

    std::vector<int> num_bytes_prefix_indications;
    num_bytes_prefix_indications.reserve(compression_variables_.size());

    std::vector<int> num_bytes_prefix_lengths;
    num_bytes_prefix_lengths.reserve(compression_variables_.size());

    std::vector<int> num_bytes_prefix_encodings;
    num_bytes_prefix_encodings.reserve(compression_variables_.size());

    std::vector<int> num_bytes_run_length_suffix_indicator;
    num_bytes_run_length_suffix_indicator.reserve(compression_variables_.size());

    std::vector<int> num_bytes_suffix_lenghts;
    num_bytes_suffix_lenghts.reserve(compression_variables_.size());

    std::vector<int> num_bytes_suffix_encodings;
    num_bytes_suffix_encodings.reserve(compression_variables_.size());

    std::vector<std::vector<uint8_t>> serialized_variables;
    serialized_variables.reserve(compression_variables_.size());

    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        auto [nb_prefix_indication_bytes,
              nb_prefix_lengths, nb_prefix_encodings,
              nb_run_length_suffix_indicator,
              nb_suffix_lengths, nb_suffix_encodings,
              serialized_variable] = var_iter->GatherSerializedCompressionDataEGU();
        
        num_bytes_prefix_indications.push_back(nb_prefix_indication_bytes);
        num_bytes_prefix_lengths.push_back(nb_prefix_lengths);
        num_bytes_prefix_encodings.push_back(nb_prefix_encodings);
        num_bytes_run_length_suffix_indicator.push_back(nb_run_length_suffix_indicator);
        num_bytes_suffix_lenghts.push_back(nb_suffix_lengths);
        num_bytes_suffix_encodings.push_back(nb_suffix_encodings);
        serialized_variables.push_back(std::move(serialized_variable));
    }

    /* Create an array holding the IDs of the variables and dimensions */
    int variable_ids[compression_variables_.size()];
    int variable_dim_ids[compression_variables_.size()];

    /* Create a new netCDF File */
    int ncid;
    /* We need to create a CDF-5 file in order to use unisgned data types (needed for the encoded mesh refinements) */
    int err = nc__create(file_name.c_str(), NC_CLOBBER|NC_CDF5, NC_SIZEHINT_DEFAULT, NULL, &ncid);
    NcCheckError(err);

    const std::string var_dim_name{"crb"}; // number of 'compressed elements'

    /* Write a dimension for each variable and it's refinement bits */
    int index = 0;
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter, ++index)
    {
        const std::string var_dimension_name = var_dim_name + std::to_string(index);

        //const size_t var_dim_length = (var_iter->size() != 0 ? var_iter->size() : 1);
        const size_t var_dim_length = serialized_variables[index].size();

        /* Define a dimension for the comrpessed elements */
        err = nc_def_dim(ncid, var_dimension_name.c_str(), var_dim_length, variable_dim_ids + index);
        NcCheckError(err);

        /* Define the variable holding the compressed data */
        err = nc_def_var(ncid, var_iter->GetName().c_str(), CmcTypeToNcType(CmcType::Uint8_t), 1, &(variable_dim_ids[index]), variable_ids + index);
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
            err = nc_put_att(ncid, variable_ids[index], "pcl", NC_INT, 1, &pre_compression_layout);
            //err = nc_put_att(ncid, variable_ids[index], "p", NC_INT, 1, &pre_compression_layout);
            NcCheckError(err);
        }

        #if 0
        if (perform_default_lossy_compression_)
        {
            /* In case a default lossy comrpession has been applied before, we save the num bytes for the initial mesh reconstruction */
            err = nc_put_att(ncid, variable_ids[index], "rfi", NC_INT, 1, &num_bytes_refinement_indications[index]);
            NcCheckError(err);
        }
        #endif

        /* Set the number of bytes for the prefix indication*/
        err = nc_put_att(ncid, variable_ids[index], "pfi", NC_INT, 1, &num_bytes_prefix_indications[index]);
        NcCheckError(err);
        
        /* Set the number of bytes for the prefix lengths */
        err = nc_put_att(ncid, variable_ids[index], "pfl", NC_INT, 1, &num_bytes_prefix_lengths[index]);
        NcCheckError(err);

        /* Set the number of bytes for the prefix encodings*/
        err = nc_put_att(ncid, variable_ids[index], "pfe", NC_INT, 1, &num_bytes_prefix_encodings[index]);
        NcCheckError(err);

        /* Set the number of bytes for the prefix encodings*/
        err = nc_put_att(ncid, variable_ids[index], "sfi", NC_INT, 1, &num_bytes_run_length_suffix_indicator[index]);
        NcCheckError(err);

        /* Set the number of bytes for the suffix lengths */
        err = nc_put_att(ncid, variable_ids[index], "sfl", NC_INT, 1, &num_bytes_suffix_lenghts[index]);
        NcCheckError(err);

        /* Set the number of bytes for the suffix encodings*/
        err = nc_put_att(ncid, variable_ids[index], "sfe", NC_INT, 1, &num_bytes_suffix_encodings[index]);
        NcCheckError(err);
    }

    /* All dimensions and variables have been defined */
    /* Therefore, we are leaving the 'define-mode' and switch to the data mode */
    err = nc_enddef(ncid);
    NcCheckError(err);

    int var_index = 0;
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter, ++var_index)
    {
        /* Write the variable's data */
        err = nc_put_var(ncid, variable_ids[var_index], static_cast<void*>(serialized_variables[var_index].data()));
        NcCheckError(err);
    }

    /* All data has been written. Therefore, the file may be closed */
    err = nc_close(ncid);
    NcCheckError(err);

    #if 0
    /* For testing purposes try to decompress directly */
    PrefixDecompressionData pd_data;

    pd_data.num_bytes_refinement_indications = -1; //Needs to be completely deleted
    pd_data.serialized_variable = serialized_variables[0];
    pd_data.num_bytes_prefix_indications = num_bytes_prefix_indications[0];
    pd_data.num_bytes_prefix_lengths = num_bytes_prefix_lengths[0];
    pd_data.num_bytes_prefix_encodings = num_bytes_prefix_encodings[0];
    pd_data.num_bytes_run_length_suffix_indicator = num_bytes_run_length_suffix_indicator[0];
    pd_data.num_bytes_suffix_lengths = num_bytes_suffix_lenghts[0];
    pd_data.num_bytes_suffix_encodings = num_bytes_suffix_encodings[0];
    pd_data.domain = compression_variables_.begin()->GetGlobalDomain();
    pd_data.current_layout = static_cast<DataLayout>(1); 
    pd_data.missing_value = static_cast<float>(-32767.0);
    //pd_data.missing_value = static_cast<float>(-9.0E+33); //MPTRAC missing value
    pd_data.type = CmcType::Float;
    pd_data.precompression_layout = static_cast<DataLayout>(1); 
    //pd_data.id = 10;//MPTRAC t
    pd_data.id = 5;
    std::vector<PrefixDecompressionData> pd_vec;
    pd_vec.push_back(std::move(pd_data));

    DecompressVariableEGU<4>(pd_vec);

    #endif
}



}
