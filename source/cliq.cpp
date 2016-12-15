#include "cliq.hpp"
#include "err_util.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>

using std::cout;
using std::cerr;

using std::unordered_map;
using std::vector;

using std::runtime_error;
using std::logic_error;
using std::exception;

using lpclique = CCtsp_lpclique;
using lpcut_in = CCtsp_lpcut_in;

namespace CMR {
namespace Sep {

Clique::Clique(const lpclique &cc_clq,
               const vector<int> &saved_tour, const vector<int> &saved_perm,
               const vector<int> &current_tour)
try
{
    for (int i = 0; i < cc_clq.segcount; ++i) {
        segment seg(cc_clq.nodes[i].lo, cc_clq.nodes[i].hi);
        bool range_agrees = true;

        for (int k = seg.start; k <= seg.end; ++k) {
            if (saved_tour[k] != current_tour[k]) {
                range_agrees = false;
                break;
            }
        }

        if (range_agrees) {
            seglist.push_back(seg);
            continue;
        }

        vector<int> seg_nodes;
        seg_nodes.reserve(seg.size());
        
        for (int k = seg.start; k <= seg.end; ++k)
            seg_nodes.push_back(saved_perm[current_tour[k]]);

        std::sort(seg_nodes.begin(), seg_nodes.end());

        int k = 0;

        while (k < seg_nodes.size()) {
            int low = seg_nodes[k];
            
            while ((k < (seg_nodes.size() - 1)) &&
                   (seg_nodes[k + 1] == (seg_nodes[k] + 1)))
                ++k;
            
            seglist.push_back(segment(low, seg_nodes[k++]));
        }
    }

    std::sort(seglist.begin(), seglist.end(), std::greater<segment>());
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Clique CCtsp_lpclique constructor failed.");
}

Clique::Clique(int start, int end,
               const vector<int> &saved_tour, const vector<int> &saved_perm,
               const vector<int> &current_tour) try
{
    bool range_agrees = true;
    segment seg(start, end);
    
    for (int k = start; k <= end; ++k)
        if (saved_tour[k] != current_tour[k]) {
            range_agrees = false;
            break;
        }

    if (range_agrees) {
        seglist.push_back(seg);
        return;
    }

    vector<int> seg_nodes;
    seg_nodes.reserve(seg.size());

    for (int k = seg.start; k <= seg.end; ++k)
        seg_nodes.push_back(saved_perm[current_tour[k]]);

    std::sort(seg_nodes.begin(), seg_nodes.end());

    int k = 0;

    while (k < seg_nodes.size()) {
        int low = seg_nodes[k];
            
        while ((k < (seg_nodes.size() - 1)) &&
               (seg_nodes[k + 1] == (seg_nodes[k] + 1)))
            ++k;
            
        seglist.push_back(segment(low, seg_nodes[k++]));
    }    
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Clique seg constructor failed.");
}

Clique::Clique(std::vector<int> &nodes, const std::vector<int> &perm)
try
{
    if (nodes.empty())
        throw logic_error("Tried to construct empty clique.");

    std::sort(nodes.begin(), nodes.end(),
              [&perm](int n1, int n2) -> bool {
                  return perm[n1] < perm[n2];
              });    

    int i = 0;

    while (i < nodes.size()) {
        int low = perm[nodes[i]];

        while ((i < (nodes.size() - 1)) &&
               (perm[nodes[i + 1]] == (perm[nodes[i]] + 1)))
            ++i;

        seglist.push_back(segment(low, perm[nodes[i++]]));        
    }

    std::sort(seglist.begin(), seglist.end(), std::greater<segment>());
} catch (const exception &e) {
    cerr << e.what() << "\n";
    throw runtime_error("Clique nodelist constructor failed.");
}

vector<int> Clique::node_list(const vector<int> &saved_tour) const
{
    vector<int> result;

    for (const segment &seg : seglist)
        for (int k = seg.start; k <= seg.end; ++k) 
            result.push_back(saved_tour[k]);

    return result;
}

CliqueBank::CliqueBank(const vector<int> &tour, const vector<int> &perm)
try : saved_tour(tour), saved_perm(perm) {} catch (const exception &e) {
    throw runtime_error("CliqueBank constructor failed.");
}

Clique::Ptr CliqueBank::add_clique(const Clique &clq)
{
    if (bank.count(clq) == 0)
        bank[clq] = std::make_shared<Clique>(clq);

    return bank[clq];
}

Clique::Ptr CliqueBank::add_clique(const CCtsp_lpclique &cc_clq,
                                    const vector<int> &tour)
{
    return add_clique(Clique(cc_clq, saved_tour, saved_perm, tour));
}

Clique::Ptr CliqueBank::add_clique(int start, int end, const vector<int> &tour)
{
    return add_clique(Clique(start, end, saved_tour, saved_perm, tour));
}

Clique::Ptr CliqueBank::add_clique(vector<int> &nodes)
{
    return add_clique(Clique(nodes, saved_perm));
}

void CliqueBank::del_clique(Clique::Ptr &clq_ptr)
{
    if (!clq_ptr)
        return;

    Clique &clq = *clq_ptr;

    CliqueBank::Itr find_it = bank.find(clq);

    if (find_it == bank.end())
        return;

    clq_ptr.reset();

    if (find_it->second.use_count() <= 1)
        bank.erase(find_it);
}

}
}
