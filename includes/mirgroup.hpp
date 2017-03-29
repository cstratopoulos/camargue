/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Utility structures for running safe Gomory cut separation.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_MIRGROUP_H
#define CMR_MIRGROUP_H

#include "config.hpp"

#if CMR_HAVE_SAFEGMI

#ifndef DO_SAFE_MIR_DBL
#define DO_SAFE_MIR_DBL 1
#endif

#ifndef SAFE_MIR_DEBUG_LEVEL
#define SAFE_MIR_DEBUG_LEVEL DBG_LEVEL_HIGH
#endif

#include <safemir/src/cutmaster_slvr.hpp>
#include <safemir/src/ds_cuts.hpp>

#include "util.hpp"

#include <memory>
#include <vector>

namespace CMR {
namespace Sep {

/**@name Safe MIR info deleters.
 * A collection of structs for unique_ptr deleters to manage safemir data
 * structures, as per //http://stackoverflow.com/a/17906532/6516346
 */
///@{

/// Deleter for constraint matrix system.
struct SystemDeleter {
    void operator()(CUTSsystem_t<double> *P) {
        CUTSfreeSystem<double>(&P);
    }
};

//// Deleter for variable info.
struct VinfoDeleter {
    void operator()(CUTSvarInfo_t<double> *P) {
        CUTSfreeVarInfo<double>(&P);
    }
};

///@}


/// Memory-managed access to classes needed during safe GMI separation.
/// See LP::Relaxation::init_mir_data for what is effectively the constructor
/// for this struct.
struct MIRgroup {
    MIRgroup()
    {
        settings.do_mir = true;
        settings.use_log = false;

        settings.max_rank = 1;
        settings.maxcutsperround = 500;

        settings.mir_cutsperfather = 10;
        settings.mir_cutsperround = 1000;
        settings.mir_maxscalings = 15;
        settings.mir_tmin = 1;
        settings.mir_tmax = 3;

        settings.save_cuts = 0;
    }

    SLVRcutterSettings_t settings;

    std::unique_ptr<CUTSsystem_t<double>, SystemDeleter> constraint_matrix;

    std::unique_ptr<CUTSsystem_t<double>, SystemDeleter> tableau_rows;

    std::unique_ptr<CUTSvarInfo_t<double>, VinfoDeleter> var_info;

    util::c_array_ptr<double> full_x;

    std::vector<char> vartype;
    std::vector<double> vranking;

};

}
}

#endif //CMR_HAVE_SAFEGMI
#endif //CMR_MIRGROUP_H
