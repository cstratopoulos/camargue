/* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */
/** @file
 * @brief Separating with cut metamorphoses.
 *
\* * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * */

#ifndef CMR_META_SEP_H
#define CMR_META_SEP_H

#include "active_tour.hpp"
#include "hypergraph.hpp"
#include "datagroups.hpp"
#include "cc_lpcuts.hpp"

#include <iostream>
#include <unordered_map>
#include <utility>
#include <vector>

namespace CMR {
namespace Sep {

class MetaCuts {
public:
    MetaCuts(const ExternalCuts &EC_,
             const std::vector<Graph::Edge> &core_edges_,
             const LP::ActiveTour &active_tour_,
             Data::SupportGroup &s_dat);

    /// Categories of implemented cut metamorphoses.
    enum class Type {
        Decker,
        Handling,
        Teething,
    };

    /// Choose what kind of cuts will be found by find_cuts.
    void set_type(Type t) { meta_type = t; }

    bool find_cuts();

    LPcutList &metacuts_q() { return meta_q; }

    using HGitr = std::vector<HyperGraph>::const_iterator;

    bool verbose = true;
    bool filter_primal = true;

private:
    Type meta_type = Type::Decker;

    /// Price all comb-like cuts in LP; return true iff some are of interest.
    bool price_combs();

    /// Used by attempt_sep to determine if \p num_paths is high enough.
    bool above_threshold(int num_paths);

    bool attempt_sep(); //!< Should separation be attempted at all.

    void price_cliques(); //!< Price all the Cliques in the LP.

    bool pure_comb(CCtsp_lpcut_in &c); //!< Is \p c a pure comb.

    std::unordered_map<Clique, double> clique_vals;
    std::vector<HGitr> interest_combs;

    const ExternalCuts &EC;

    const std::vector<Graph::Edge> &core_edges;
    const LP::ActiveTour &active_tour;
    Data::SupportGroup &supp_data;

    TourGraph TG;
    std::vector<int> perm_elist;

    LPcutList meta_q;
};

std::ostream &operator<<(std::ostream &os, MetaCuts::Type t);

}
}

#endif
