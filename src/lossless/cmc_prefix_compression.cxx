#include "lossless/cmc_prefix_compression.hxx"
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

#include <vector>

namespace cmc
{

namespace prefix
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
Compressor::SetSplitVariable(const SplitVariable& split_var)
{
    split_variables_.push_back(split_var);
}

void
Compressor::SetSplitVariable(SplitVariable&& split_var)
{
    split_variables_.push_back(std::move(split_var));
}

void
Compressor::Setup()
{
    if (split_variables_.size() > 0)
    {
        CompressionSettings settings;
        settings.SetSplitVariables(split_variables_);

        compression_data_->SetCompressionSettings(std::move(settings));
        compression_data_->SplitVariables();
    }

    compression_data_->CheckConsistencyOfInputVariables();

    compression_data_->BuildInitialMesh();

    compression_data_->DistributeDataOnInitialMesh();

    compression_data_->ApplyScalingAndOffset();

    compression_data_->TransformInputToCompressionVariables();

    compression_data_->FilterDataAsDifferences();

    compression_variables_ = compression_data_->GetByteVariablesForCompression();

    compression_data_.reset(nullptr); 
}


void
Compressor::Compress()
{
    /* Afterwards, we create prefixes in the tree hierachy */
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        //var_iter->XORConsecutiveValues();
        /* (Try to) release the initial data which already has been transformed */
        var_iter->KeepInitialData(false);
        //var_iter->PrintCompressionValues();
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
    NcWriter writer(file_name, NC_NETCDF4);
    writer.ReserveVariables(compression_variables_.size());
    
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        writer.AddVariable(var_iter->WriteCompressedData(time_step, SuffixEncoding::Plain));
    }

    writer.AddGlobalAttribute(NcAttribute(kCompressionSchemeAttrName, CmcUniversalType(static_cast<CompressionSchemeType>(CompressionScheme::PrefixExtraction))));
    writer.Write();
}

}

}
