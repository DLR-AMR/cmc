#ifndef CMC_T8_MPI_HXX
#define CMC_T8_MPI_HXX
/**
 * @file cmc_t8_morton.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <iterator>

namespace cmc
{


class DataOffsets
{
    using iterator = std::vector<MortonIndex>::iterator;
    using const_iterator = std::vector<MortonIndex>::const_iterator;
public:

    DataOffsets() = delete;
    DataOffsets(const size_t num_elements)
    : offsets_(num_elements){};

    DataOffsets(const DataOffsets& other) = default;
    DataOffsets& operator=(const DataOffsets& other) = default;
    DataOffsets(DataOffsets&& other) = default;
    DataOffsets& operator=(DataOffsets&& other) = default;

    template<typename U> auto operator[](U index) const -> std::enable_if_t<std::is_integral_v<U>, MortonIndex> {return offsets_[index];};
    template<typename U> auto operator[](U index) -> std::enable_if_t<std::is_integral_v<U>, MortonIndex&> {return offsets_[index];};

    iterator Begin() { return offsets_.begin(); };
    iterator End() { return offsets_.end(); };
    const_iterator Begin() const { return offsets_.begin(); };
    const_iterator End() const { return offsets_.end(); };
    const_iterator CBegin() const { return offsets_.cbegin(); };
    const_iterator CEnd() const { return offsets_.cend(); };

    size_t size() const {return offsets_.size();};

    MortonIndex* data()
    {
        return offsets_.data();
    };
private:
    std::vector<MortonIndex> offsets_;
};

DataOffsets
GatherGlobalDataOffsets(const AmrMesh mesh, const MPI_Comm comm);

}

#endif /* !CMC_T8_MPI_HXX */
