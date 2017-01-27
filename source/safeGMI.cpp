#include "safeGMI.hpp"

#if CMR_HAVE_SAFEGMI

#ifndef DO_SAFE_MIR_DBL
#define DO_SAFE_MIR_DBL 1
#endif

#ifndef SAFE_MIR_DEBUG_LEVEL
#define SAFE_MIR_DEBUG_LEVEL DBG_LEVEL_HIGH
#endif

#include <safemir/src/util_cuts.hpp>
#include <safemir/src/cutmaster_slvr.cpp>
#include <safemir/src/ds_cuts.cpp>
#include <safemir/src/gen_mir.cpp>

#include "util.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using std::unique_ptr;

using std::vector;

using std::cout;
using std::cerr;

using std::runtime_error;
using std::exception;

namespace CMR {
namespace Sep {

SafeGomory::SafeGomory(LP::Relaxation &rel, const vector<double> &_tour_edges,
                       const vector<double> &lp_edges) try
    : lp_relax(rel), gmi_q(15), tour_edges(_tour_edges), frac_x(lp_edges)
{
    lp_relax.init_mir_data(mir_data);
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("SafeGomory constructor failed.");
}

bool SafeGomory::find_cuts()
{
    runtime_error err("Problem in SafeGomory::find_cuts.");

    int num_added = 0;
    int total_num_added = 0;

    using sp_rowlist = CUTSsprowlist_t<double>;
    using rowlist_elem = CUTSrowListElem_t<double>;
    using sp_row = CUTSsprow_t<double>;
    
    struct RowlistDeleter {
        void operator()(sp_rowlist **p) {
            CUTSfreeRowList(p);
        }
    };

    sp_rowlist **gen_cuts;
    if (CUTSnewRowList(gen_cuts)) {
        cerr << "CUTSnewRowList failed.\n";
        throw err;
    }

    if (SLVRcutter_iter(0,
                        &mir_data.settings,
                        *mir_data.constraint_matrix,
                        *mir_data.tableau_rows,
                        &total_num_added, &num_added,
                        mir_data.full_x.get(),
                        *mir_data.var_info,
                        false, NULL,
                        lp_relax.num_cols(),
                        *gen_cuts,
                        &mir_data.vranking[0])) {
        CUTSfreeRowList(gen_cuts);
        cerr << "SLVRcutter_iter failed.\n";
        throw err;
    }

    unique_ptr<sp_rowlist *, RowlistDeleter> generated_cuts(gen_cuts);

    if ((*generated_cuts)->size == 0) {
        cout << "\tNo safe Gomory cuts generated.\n";
        return false;
    }

    using SparseRow = Sep::SparseRow;
    vector<SparseRow> primal_found;

    for(rowlist_elem *it = (*generated_cuts)->first; it; it = it->next) {
        sp_row *row = it->row;
        double rhs = row->rhs;
        int nz = row->nz;

        if (row->sense != 'G') {
            cerr << "Somehow got non >= cut.\n";
            throw err;
        }

        double lp_viol = 0.0;
        double tour_act = 0.0;

        if (CUTScomputeViolation(row, &frac_x[0], &lp_viol) ||
            CUTScomputeActivity(row, &tour_edges[0], &tour_act)) {
            cerr << "CUTScomputeActivity or CUTScomputeViolation failed.\n";
            throw err;
        }

        if (tour_act == rhs && lp_viol >= Epsilon::Cut) {
            SparseRow primal_row;
            cout << "\tFound GMI cut with viol " << lp_viol << "\n";
            try {
                primal_row.rmatind.resize(nz);
                primal_row.rmatval.resize(nz);
                primal_row.sense = 'G';
                primal_row.rhs = rhs;
                primal_row.lp_viol = lp_viol;

                for (int i = 0; i < nz; ++i) {
                    primal_row.rmatind[i] = row->rowind[i];
                    primal_row.rmatval[i] = row->rowval[i];
                }

                primal_found.push_back(primal_row);
            } CMR_CATCH_PRINT_THROW("allocating/pushing primal row", err);
        }
    }

    if (primal_found.empty()) {
        cout << "\tFound safe Gomory cuts but none were tight.\n";
        return false;
    }

    // sort by lexicographically preferring violation, then sparsity
    int numcols = lp_relax.num_cols();
    std::sort(primal_found.begin(), primal_found.end(),
              [numcols](const SparseRow &a, const SparseRow &b)
              {
                  return
                  std::make_tuple(a.lp_viol, numcols - a.rmatind.size())
                  >
                  std::make_tuple(b.lp_viol, numcols - b.rmatind.size());
              });

    primal_found.resize(gmi_q.q_capacity);

    try {
        for (SparseRow &a : primal_found)
            gmi_q.push_back(std::move(a));
    } CMR_CATCH_PRINT_THROW("putting found cuts in cut q", err);

    return true;    
}

}
}

#endif //CMR_HAVE_SAFEGMI
