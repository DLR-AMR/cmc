#include "lossy/cmc_full_amr_compression.hxx"
#include "t8code/cmc_t8_adapt_callbacks.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_vector_view.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_reader.hxx"
#include "netcdf/cmc_nc_writer.hxx"
#endif

namespace cmc
{

namespace full_amr
{

int
Compressor::GetMpiSize() const
{
    int comm_size{1};
    int err = MPI_Comm_size(comm_, &comm_size);
    MPICheckError(err);
    return comm_size;
}

void
Compressor::Setup()
{
    compression_data_->SplitVariables();

    compression_data_->CheckConsistencyOfInputVariables();

    compression_data_->BuildInitialMesh();

    compression_data_->DistributeDataOnInitialMesh();

    compression_data_->ApplyScalingAndOffset();

    compression_data_->SetupVariablesForCompression();
    //compression_data_->TransformInputToCompressionVariables();

    compression_data_->IndicateToKeepInitialMeshAndData();
}


void
Compressor::Compress()
{
    /* Most liekly the compression on the inital data is only an increase of a few percent in spared elements, but due to the search on the large forest it takes more time */
    /* Perform the default lossy compression */
    //compression_data_->CompressByAdaptiveCoarseningWithRegardToInitialData();
    /* Perform the default lossy compression */
    compression_data_->CompressByAdaptiveCoarsening(CompressionMode::OneForOne);
    cmc_debug_msg("\n\nAdaptive Coarsening has been applied\n\n");
    //compression_data_->WriteVTKFilePerVariable("ac_t2m_data");

    /* Transform the data to byte variables */
    compression_variables_ = compression_data_->GetByteVariablesForCompression();
    cmc_debug_msg("\n\nVariables have been transformed to byte variables\n\n");
    //std::exit(1);
    /* Get the coarsening/refienment bits for all variables */
    ac_indications_ = compression_data_->TransferIndicationBits();
    cmc_debug_msg("\n\nTransfer Indication bits has been performed \n\n");
    /* Afterwards, we do not need the compression_data anymore since we have transfered all necessary data */
    compression_data_.reset(nullptr);

    /* Store the initial compressed mesh */
    for (auto cv_iter = compression_variables_.begin(); cv_iter != compression_variables_.end(); ++cv_iter)
    {
        cv_iter->StoreInitialMesh();
    }
    cmc_debug_msg("\n\nInitial mesh storage has been indicated\n\n");
    /* Perform trail truncation until the error trhesholds are exhausted */
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        //var_iter->PerformTailTruncationRegardingUncompressedStates();
        var_iter->PerformTailTruncationOnInitialData();
        //var_iter->PerformTailTruncation();
    }
    cmc_debug_msg("\n\nTail truncation has been applied\n\n");

    /* Afterwards, we create prefixes in the tree hierachy */
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        /* We create the adapt data based on the compression settings, the forest and the variable to consider during the adaptation/coarsening */
        PrefixAdaptData adapt_data{*var_iter};

        /* Iterate until all prefixes have been extracted (up until the the root of the mesh) */
        while (adapt_data.IsCompressionProgressing())
        {
            /* Allocate memory for the prefix extraction and set up evertything needed for the coarsening process */
            adapt_data.InitializeCompressionIteration();

            /* Perform a coarsening iteration and find prefix and refinement bits */
            t8_forest_t adapted_forest = t8_forest_new_adapt(adapt_data.GetCurrentMesh(), ExtractCommonPrefixes, 0, 0, static_cast<void*>(&adapt_data));

            /* After the prefixes have been extracted and 'stored' on the coarser mesh in this iteration,
             * we update the mesh that they are defined on */
            adapt_data.SetCurrentMesh(adapted_forest);

            adapt_data.FinalizeCompressionIteration();
        }
    }
    cmc_debug_msg("\n\nCompression is completely finihsed \n\n");
}


void
Compressor::WriteCompressedData(const std::string& file_name, const int time_step) const
{
    NcWriter writer(file_name, NC_NETCDF4); //oder NC_CDF5
    writer.ReserveVariables(2 * compression_variables_.size());

    int ac_bits_index = 0;
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter, ++ac_bits_index)
    {
        /* Write the compressed byte variable */
        writer.AddVariable(var_iter->WriteCompressedData(time_step, SuffixEncoding::ArithmeticLengthEncoding));
        /* Write the refinement bits from the adaptive coarsening steps */
        //writer.AddVariable(CreateRefinementBitsVariable(*var_iter, time_step, ac_indications_[ac_bits_index].ac_indicator_bits));
    }

    writer.AddGlobalAttribute(NcAttribute(kCompressionSchemeAttrName, CmcUniversalType(static_cast<CompressionSchemeType>(CompressionScheme::PrefixExtraction))));
    writer.Write();
}


}
}
