#ifndef CMR_MIRGROUP_H
#define CMR_MIRGROUP_H

#include "config.hpp"

#if CMR_HAVE_SAFEGMI

#include "util.hpp"

#include <memory>
#include <vector>

#ifndef DO_SAFE_MIR_DBL
#define DO_SAFE_MIR_DBL 1
#endif

#ifndef SAFE_MIR_DEBUG_LEVEL
#define SAFE_MIR_DEBUG_LEVEL DBG_LEVEL_HIGH
#endif

#include <safemir/src/cutmaster_slvr.hpp>
#include <safemir/src/ds_cuts.hpp>

namespace CMR {
namespace Sep {

//http://stackoverflow.com/a/17906532/6516346
struct SystemDeleter {
    void operator()(CUTSsystem_t<double> *P) {
        CUTSfreeSystem<double>(&P);
    }
};

struct VinfoDeleter {
    void operator()(CUTSvarInfo_t<double> *P) {
        CUTSfreeVarInfo<double>(&P);
    }
};

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
