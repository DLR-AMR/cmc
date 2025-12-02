#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"

namespace cmc
{

[[noreturn]] void cmc_exit(const char* _err_msg, const char* _location)
{
    std::cout << "CMC_EXIT is invoked..." << std::endl << _err_msg << std::endl << "Error Occurence:"  << _location << std::endl;
    std::exit(EXIT_FAILURE);
}

[[noreturn]] void cmc_abort(const char* _err_msg, const char* _location)
{
    std::cout << "CMC_ABORT is invoked..." << std::endl << _err_msg << std::endl << "Error Occurence:"  << _location << std::endl;
    std::abort();
}


int
GetDimensionalityOfDataLayout(const DataLayout layout)
{
    if (layout > DataLayout::LayoutUndefined &&
        layout < DataLayout::_InternEnd1DLayouts)
    {
        /* If one-diemnsional */
        return 1;
    } else if (layout > DataLayout::_InternEnd1DLayouts &&
               layout < DataLayout::_InternEnd2DLayouts)
    {
        /* If two-dimensional */
        return 2;
    } else if (layout > DataLayout::_InternEnd2DLayouts &&
               layout < DataLayout::_InternEnd3DLayouts)
    {
        /* If three-dimensional */
        return 3;
    } else if (layout > DataLayout::_InternEnd3DLayouts &&
               layout < DataLayout::_InternEnd4DLayouts)
    {
        /* If four-dimensional */
        return 4;
    } else
    {
        /* Return an error if a unsupported DataLayout has been supplied */
        cmc_err_msg("The dimensionality of the supplied layout could not be retrieved.");
        return -1;
    }
}


std::vector<Dimension>
GetDimensionVectorFromLayout(const DataLayout layout)
{
    switch (layout)
    {
        /* 1D layouts */
        case DataLayout::Lon_:
            return std::vector<Dimension>{Dimension::Lon};
        break;
        case DataLayout::Lat_:
            return std::vector<Dimension>{Dimension::Lat};
        break;
        case DataLayout::Lev_:
            return std::vector<Dimension>{Dimension::Lev};
        break;
        case DataLayout::Time_:
            return std::vector<Dimension>{Dimension::Time};
        break;
        /* 2D layouts */
        case DataLayout::Lat_Lon:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lon};
        break;
        case DataLayout::Lon_Lat:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lat};
        break;
        case DataLayout::Lat_Lev:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lev};
        break;
        case DataLayout::Lev_Lat:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lat};
        break;
        case DataLayout::Lon_Lev:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lev};
        break;
        case DataLayout::Lev_Lon:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lon};
        break;
        /* 3D layouts */
        case DataLayout::Lat_Lon_Lev:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lon, Dimension::Lev};
        break;
        case DataLayout::Lat_Lev_Lon:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lev, Dimension::Lon};
        break;
        case DataLayout::Lev_Lat_Lon:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lat, Dimension::Lon};
        break;
        case DataLayout::Lev_Lon_Lat:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lon, Dimension::Lat};
        break;
        case DataLayout::Lon_Lev_Lat:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lev, Dimension::Lat};
        break;
        case DataLayout::Lon_Lat_Lev:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lat, Dimension::Lev};
        break;
        /* 4D layouts */
        case DataLayout::Time_Lev_Lat_Lon:
            return std::vector<Dimension>{Dimension::Time, Dimension::Lev, Dimension::Lat, Dimension::Lon};
        break;
        case DataLayout::Time_Lev_Lon_Lat:
            return std::vector<Dimension>{Dimension::Time, Dimension::Lev, Dimension::Lon, Dimension::Lat};
        break;
        case DataLayout::Time_Lat_Lev_Lon:
            return std::vector<Dimension>{Dimension::Time, Dimension::Lat, Dimension::Lev, Dimension::Lon};
        break;
        case DataLayout::Time_Lat_Lon_Lev:
            return std::vector<Dimension>{Dimension::Time, Dimension::Lat, Dimension::Lon, Dimension::Lev};
        break;
        case DataLayout::Time_Lon_Lev_Lat:
            return std::vector<Dimension>{Dimension::Time, Dimension::Lon, Dimension::Lev, Dimension::Lat};
        break;
        case DataLayout::Time_Lon_Lat_Lev:
            return std::vector<Dimension>{Dimension::Time, Dimension::Lon, Dimension::Lat, Dimension::Lev};
        break;
        case DataLayout::Lev_Time_Lat_Lon:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Time, Dimension::Lat, Dimension::Lon};
        break;
        case DataLayout::Lev_Time_Lon_Lat:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Time, Dimension::Lon, Dimension::Lat};
        break;
        case DataLayout::Lev_Lat_Time_Lon:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lat, Dimension::Time, Dimension::Lon};
        break;
        case DataLayout::Lev_Lat_Lon_Time:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lat, Dimension::Lon, Dimension::Time};
        break;
        case DataLayout::Lev_Lon_Time_Lat:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lon, Dimension::Time, Dimension::Lat};
        break;
        case DataLayout::Lev_Lon_Lat_Time:
            return std::vector<Dimension>{Dimension::Lev, Dimension::Lon, Dimension::Lat, Dimension::Time};
        break;
        case DataLayout::Lat_Time_Lev_Lon:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Time, Dimension::Lev, Dimension::Lon};
        break;
        case DataLayout::Lat_Time_Lon_Lev:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Time, Dimension::Lon, Dimension::Lev};
        break;
        case DataLayout::Lat_Lev_Time_Lon:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lev, Dimension::Time, Dimension::Lon};
        break;
        case DataLayout::Lat_Lev_Lon_Time:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lev, Dimension::Lon, Dimension::Time};
        break;
        case DataLayout::Lat_Lon_Time_Lev:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lon, Dimension::Time, Dimension::Lev};
        break;
        case DataLayout::Lat_Lon_Lev_Time:
            return std::vector<Dimension>{Dimension::Lat, Dimension::Lon, Dimension::Lev, Dimension::Time};
        break;
        case DataLayout::Lon_Time_Lev_Lat:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Time, Dimension::Lev, Dimension::Lat};
        break;
        case DataLayout::Lon_Time_Lat_Lev:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Time, Dimension::Lat, Dimension::Lev};
        break;
        case DataLayout::Lon_Lev_Time_Lat:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lev, Dimension::Time, Dimension::Lat};
        break;
        case DataLayout::Lon_Lev_Lat_Time:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lev, Dimension::Lat, Dimension::Time};
        break;
        case DataLayout::Lon_Lat_Time_Lev:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lat, Dimension::Time, Dimension::Lev};
        break;
        case DataLayout::Lon_Lat_Lev_Time:
            return std::vector<Dimension>{Dimension::Lon, Dimension::Lat, Dimension::Lev, Dimension::Time};
        break;
        default:
            cmc_err_msg("The dimension vector could not be retrieved from the DataLayout since it is not recognized.");
            return std::vector<Dimension>();
    }
}

DataLayout
GetDataLayoutAfterRemoval(const DataLayout initial_layout, const Dimension removed_dimension)
{
    /* Currently, it is only possible to remove a dimension from a three dimensional layout
     * since two is the minimum dimenison of data to be compressed */
    cmc_assert(initial_layout > DataLayout::_InternEnd2DLayouts && initial_layout < DataLayout::_InternEnd3DLayouts);

    switch (initial_layout)
    {
        case DataLayout::Lat_Lon_Lev:
            switch (removed_dimension)
            {
                case Dimension::Lev:
                    return DataLayout::Lat_Lon;
                break;
                case Dimension::Lat:
                    return DataLayout::Lon_Lev;
                break;
                case Dimension::Lon:
                    return DataLayout::Lat_Lev;
                break;
                default:
                    cmc_err_msg("An undefined dimension has been supplied.");
                    return DataLayout::LayoutUndefined;
            }
        break;
        case DataLayout::Lat_Lev_Lon:
            switch (removed_dimension)
            {
                case Dimension::Lev:
                    return DataLayout::Lat_Lon;
                break;
                case Dimension::Lat:
                    return DataLayout::Lev_Lon;
                break;
                case Dimension::Lon:
                    return DataLayout::Lat_Lev;
                break;
                default:
                    cmc_err_msg("An undefined dimension has been supplied.");
                    return DataLayout::LayoutUndefined;
            }
        break;
        case DataLayout::Lev_Lat_Lon:
            switch (removed_dimension)
            {
                case Dimension::Lev:
                    return DataLayout::Lat_Lon;
                break;
                case Dimension::Lat:
                    return DataLayout::Lev_Lon;
                break;
                case Dimension::Lon:
                    return DataLayout::Lev_Lat;
                break;
                default:
                    cmc_err_msg("An undefined dimension has been supplied.");
                    return DataLayout::LayoutUndefined;
            }
        break;
        case DataLayout::Lev_Lon_Lat:
            switch (removed_dimension)
            {
                case Dimension::Lev:
                    return DataLayout::Lon_Lat;
                break;
                case Dimension::Lat:
                    return DataLayout::Lev_Lon;
                break;
                case Dimension::Lon:
                    return DataLayout::Lev_Lat;
                break;
                default:
                    cmc_err_msg("An undefined dimension has been supplied.");
                    return DataLayout::LayoutUndefined;
            }
        break;
        case DataLayout::Lon_Lev_Lat:
            switch (removed_dimension)
            {
                case Dimension::Lev:
                    return DataLayout::Lon_Lat;
                break;
                case Dimension::Lat:
                    return DataLayout::Lon_Lev;
                break;
                case Dimension::Lon:
                    return DataLayout::Lev_Lat;
                break;
                default:
                    cmc_err_msg("An undefined dimension has been supplied.");
                    return DataLayout::LayoutUndefined;
            }
        break;
        case DataLayout::Lon_Lat_Lev:
            switch (removed_dimension)
            {
                case Dimension::Lev:
                    return DataLayout::Lon_Lat;
                break;
                case Dimension::Lat:
                    return DataLayout::Lon_Lev;
                break;
                case Dimension::Lon:
                    return DataLayout::Lat_Lev;
                break;
                default:
                    cmc_err_msg("An undefined dimension has been supplied.");
                    return DataLayout::LayoutUndefined;
            }
        break;
        default:
            cmc_err_msg("The supplied layout is not supported for removal of a dimension.");
            return DataLayout::LayoutUndefined;
    }
}

std::string
GetDimensionName(const Dimension dimension)
{
    switch (dimension)
    {
        case Dimension::Lon:
            return std::string{"lon"};
        break;
        case Dimension::Lat:
            return std::string{"lat"};
        break;
        case Dimension::Lev:
            return std::string{"lev"};
        break;
        case Dimension::Time:
            return std::string{"time"};
        break;
        default:
            cmc_err_msg("A not supported dimension has been supplied.");
            return std::string{"dim_err"};
    }
}

DataLayout
GetDefaultDataLayout(const int dimensionality)
{
    switch (dimensionality)
    {
        case 1:
            return DataLayout::Lon_;
        break;
        case 2:
            return DataLayout::Lat_Lon;
        break;
        case 3:
            return DataLayout::Lev_Lat_Lon;
        break;
        case 4:
            return DataLayout::Time_Lev_Lat_Lon;
        break;
        default:
            std::cout << "[cmc] ERROR: The specified dimensionality is not supported (only 1D, 2D and 3D)." << std::endl;
            std::exit(EXIT_FAILURE);
            return DataLayout::LayoutUndefined;
        break;
    }
}

}
