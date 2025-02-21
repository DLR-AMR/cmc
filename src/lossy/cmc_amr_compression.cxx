#include "lossy/cmc_amr_compression.hxx"
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
Compressor::Setup()
{
    compression_data_->SplitVariables();

    compression_data_->CheckConsistencyOfInputVariables();

    compression_data_->BuildInitialMesh();

    compression_data_->DistributeDataOnInitialMesh();

    compression_data_->ApplyScalingAndOffset();

    compression_data_->SetupVariablesForCompression();
    //compression_data_->TransformInputToCompressionVariables();
}


void
Compressor::Compress()
{

    /* Perform the default lossy compression */
    compression_data_->CompressByAdaptiveCoarsening(CompressionMode::OneForOne);

    //compression_data_->WriteVTKFilePerVariable("ac_t2m_data");

    //std::exit(1);
    /* Transform the data to byte variables */
    compression_variables_ = compression_data_->GetByteVariablesForCompression();
    //std::exit(1);
    /* Get the coarsening/refienment bits for all variables */
    ac_indications_ = compression_data_->TransferIndicationBits();
    
    /* Afterwards, we do not need the compression_data anymore since we have transfered all necessary data */
    compression_data_.reset(nullptr);

    /* Store the initial compressed mesh */
    for (auto cv_iter = compression_variables_.begin(); cv_iter != compression_variables_.end(); ++cv_iter)
    {
        cv_iter->StoreInitialMesh();
    }

    /* Perform trail truncation until the error trhesholds are exhausted */
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter)
    {
        var_iter->PerformTailTruncation();
    }

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

static
nc::Variable
CreateRefinementBitsVariable(const ByteVar& var, const int time_step, const std::vector<bit_map::BitMap>& levelwise_indication_bits)
{
    /* Generate a (potentially new) context information for the variable */
    const int global_context_info = (var.GetGlobalContextInformation() != kNoGlobalContext ? var.GetGlobalContextInformation() : 0);

    /* Create a netCDF variable to put out */
    std::string var_name = var.GetName() + "_" + std::to_string(global_context_info) + "_" + std::to_string(time_step) + "ref_bits";

    /* Count the bytes needed by the refinement bits */
    size_t num_bytes = 0;
    for (auto iter = levelwise_indication_bits.rbegin(); iter != levelwise_indication_bits.rend(); ++iter)
    {
        num_bytes += iter->size_bytes();
    }

    /* We add a size_t for each refinement level indicating the number of bits encoded */
    num_bytes += sizeof(size_t) * levelwise_indication_bits.size();

    nc::SpecificVariable<uint8_t> indication_bits {var_name, var.GetID()};
    indication_bits.Reserve(num_bytes);

    /* We store a size_t of the number of bits in a stream */
    std::vector<uint8_t> serialized_bit_count;
    serialized_bit_count.reserve(sizeof(size_t));

    /* We copy the data to the byte stream (in reverese in order to reconstruct it more easily during the decompression) */
    for (auto iter = levelwise_indication_bits.rbegin(); iter != levelwise_indication_bits.rend(); ++iter)
    {
        /* We serialize the bit count and store it temporarily */
        const size_t num_bits = iter->size();
        serialized_bit_count.clear();
        PushBackValueToByteStream(serialized_bit_count, num_bits);

        /* Push back the bit count */
        indication_bits.PushBack(serialized_bit_count);
        
        /* Afterwards, we are storing the actual refinement bits for this level */
        indication_bits.PushBack(iter->GetByteData());
    }

    /* Assign some attributes to it */
    std::vector<nc::Attribute> attributes;
    attributes.emplace_back("id", var.GetID());
    attributes.emplace_back("time_step", time_step);
    attributes.emplace_back("add_refinements", static_cast<int>(levelwise_indication_bits.size()));

    return nc::Variable(std::move(indication_bits), std::move(attributes));
}

void
Compressor::WriteCompressedData(const std::string& file_name, const int time_step) const
{
    nc::Writer writer(file_name, NC_NETCDF4); //oder NC_CDF5
    writer.ReserveVariables(2 * compression_variables_.size());

    int ac_bits_index = 0;
    for (auto var_iter = compression_variables_.begin(); var_iter != compression_variables_.end(); ++var_iter, ++ac_bits_index)
    {
        /* Write the compressed byte variable */
        writer.AddVariable(var_iter->WriteCompressedData(time_step, SuffixEncoding::LengthEncoding));
        /* Write the refinement bits from the adaptive coarsening steps */
        //writer.AddVariable(CreateRefinementBitsVariable(*var_iter, time_step, ac_indications_[ac_bits_index].ac_indicator_bits));
    }

    writer.AddGlobalAttribute(nc::Attribute(nc::kCompressionSchemeAttrName, CmcUniversalType(static_cast<nc::CompressionSchemeType>(nc::CompressionScheme::PrefixExtraction))));
    writer.Write();
}

}

}
