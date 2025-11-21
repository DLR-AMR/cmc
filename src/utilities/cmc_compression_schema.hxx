#ifndef CMC_COMPRESSION_SCHEMA_HXX
#define CMC_COMPRESSION_SCHEMA_HXX

namespace cmc
{

enum CompressionSchema
{
    AdaptiveCoarsening,
    TrimmedMultiResResiduals,
    PrefixExtraction,
    PrefixExtractionPlainSuffixes,
    MultiResExtraction,
    EmbeddedPrefixExtraction,
    EmbeddedPrefixExtractionTrimmedSuffixes,
    EmbeddedMultiResExtractionLinReconstruction,
    LossyEmbeddedMultiResExtractionNearestNeighborReconstruction,
    EmbeddedPrefixExtractionPlainSuffixes,
    EmbeddedQuantizedPrefixExtraction,
    EmbeddedMultiResExtraction,
    EmbeddedTrimmedMultiResExtraction,
    EmbeddedMultiResExtractionTrimmedResiduals,
    _TestEmbeddedPCP4Extraction,
    _TestPCP4Extraction,
    PatchPrefixExtractionPlainSuffixes
};

}


#endif /* !CMC_COMPRESSION_SCHEMA_HXX */
