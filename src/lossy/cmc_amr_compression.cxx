#include "lossy/cmc_amr_compression.hxx"
#include "t8code/cmc_t8_adapt_callbacks.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_span.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#include "netcdf/cmc_nc_io.hxx"
#include "netcdf/cmc_nc_reader.hxx"
#include "netcdf/cmc_nc_writer.hxx"
#endif

#include <vector>

namespace cmc
{

namespace amr
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
Compressor::Setup(const bool with_default_lossy_amr_compression)
{
    compression_data_->SplitVariables();

    compression_data_->CheckConsistencyOfInputVariables();

    compression_data_->BuildInitialMesh();

    compression_data_->DistributeDataOnInitialMesh();

    compression_data_->ApplyScalingAndOffset();

    compression_data_->SetupVariablesForCompression();

    if (!with_default_lossy_amr_compression)
    {
        /* In case a default lossy compression needs to be performed first, we have to postpone the transfomration to byte variables */
        compression_variables_ = compression_data_->GetByteVariablesForCompression();
    } else
    {
        perform_default_lossy_compression_ = true;
    }
}


void
Compressor::Compress()
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

    /* Perform trail truncation until the error thresholds are exhausted */
    //for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    //{
    //    var_iter->PerformTailTruncation();
    //}

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
}

void
Compressor::WriteCompressedData(const std::string& file_name, const int time_step) const
{
    std::vector<NcVariable> vars;
    vars.reserve(compression_variables_.size());
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        vars.push_back(var_iter->WriteCompressedData(time_step));
    }

    NcWriter writer(file_name, NC_NETCDF4); //oder NC_CDF5

    writer.AddVariable(vars.back());
    writer.AddGlobalAttribute(NcAttribute(kCompressionSchemeAttrName, CmcUniversalType(static_cast<CompressionSchemeType>(CompressionScheme::PrefixExtraction))));
    writer.Write();   
}

}

}
