#include "edgehash.hpp"

#include <stdexcept>

extern "C" {
#include <concorde/INCLUDE/util.h>
}

using std::exception;
using std::runtime_error;

using std::vector;

namespace CMR {
namespace util {

struct EdgeHash::eh_impl {
    eh_impl(int size);
    ~eh_impl();

    CCutil_edgehash *ptr() { return &eh; }

    CCutil_edgehash eh;
    static constexpr double factor{1.5};
};

EdgeHash::eh_impl::eh_impl(int size)
{
    if (CCutil_edgehash_init(&eh, factor * size))
        throw runtime_error("CCutil_edgehash_init failed");
}

EdgeHash::eh_impl::~eh_impl() { CCutil_edgehash_free(&eh); }

EdgeHash::EdgeHash(int size) try
    : eh_pimpl(util::make_unique<eh_impl>(size))
{} catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    throw std::runtime_error("EdgeHash constructor failed.");
}

EdgeHash::~EdgeHash() {}

/// Add the edge with endpoints \p end1 and \p end2 to the EdgeHash, with
/// lookup value set to \p val.
void EdgeHash::add(int end1, int end2, int val)
{
    if (CCutil_edgehash_add(eh_pimpl->ptr(), end1, end2, val))
        throw runtime_error("CCutil_edgehash_add failed.");
}

///If the edge with endpoints \p end1 \p end2 is already in the EdgeHash, set
/// its value to \p val.
void EdgeHash::set(int end1, int end2, int val)
{
    if (CCutil_edgehash_set(eh_pimpl->ptr(), end1, end2, val))
        throw runtime_error("CCutil_edgehash_set failed.");
}

/// Erase the edge with endpoints \p end1 \p end2 from the EdgeHash.
/// @remark we ignore errors from Concorde, hence do not consider it an error
/// to erase no element from an empty hash, or to attempt to erase an element
/// which was not there in the first place.
void EdgeHash::erase(int end1, int end2)
{
    CCutil_edgehash_del(eh_pimpl->ptr(), end1, end2);
}

vector<Graph::Edge> EdgeHash::get_all()
{
    int ecount = 0;
    int *elist;
    int *elen;

    if (CCutil_edgehash_getall(eh_pimpl->ptr(), &ecount, &elist, &elen))
        throw runtime_error("CCutil_edgehash_getall failed.");

    util::c_array_ptr<int> elist_owner(elist);
    util::c_array_ptr<int> elen_owner(elen);

    vector<Graph::Edge> result(ecount);

    for (int i = 0; i < ecount; ++i)
        result[i] = Graph::Edge(elist[2 * i], elist[(2 * i) +1], elen[i]);

    return result;
}

void EdgeHash::clear() { CCutil_edgehash_delall(eh_pimpl->ptr()); }

int EdgeHash::get_val(int end1, int end2)
{
    int val = 0;

    if (CCutil_edgehash_find(eh_pimpl->ptr(), end1, end2, &val) == -1)
        return -1;

    return val;
}


}
}
