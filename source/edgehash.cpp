#include "edgehash.hpp"

#include <stdexcept>

using std::exception;
using std::runtime_error;
using std::logic_error;

using std::vector;

namespace CMR {
namespace util {

EdgeHash::EdgeHash(int size) try 
{
    if (CCutil_edgehash_init(&eh, factor * size))
        throw std::runtime_error("CCutil_edgehash_init failed.");
    
} catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    throw std::runtime_error("EdgeHash constructor failed.");
}

EdgeHash::~EdgeHash() { CCutil_edgehash_free(&eh); }

void EdgeHash::add(int end1, int end2, int val)
{
    if (CCutil_edgehash_add(&eh, end1, end2, val))
        throw runtime_error("CCutil_edgehash_add failed.");
}

void EdgeHash::set(int end1, int end2, int val)
{
    if (CCutil_edgehash_set(&eh, end1, end2, val))
        throw runtime_error("CCutil_edgehash_set failed.");
}

void EdgeHash::erase(int end1, int end2)
{
    if (CCutil_edgehash_del(&eh, end1, end2))
        throw runtime_error("CCutil_edgehash_del failed.");
}

vector<Edge> EdgeHash::get_all()
{
    int ecount = 0;
    int *elist;
    int *elen;

    if (CCutil_edgehash_getall(&eh, &ecount, &elist, &elen))
        throw runtime_error("CCutil_edgehash_getall failed.");

    util::c_array_ptr elist_owner(elist);
    util::c_array_ptr elen_owner(elen);

    vector<Edge> result(ecount);

    for (int i = 0; i < ecount; ++i)
        result[i] = Edge(elist[2 * i], elist[(2 * i) +1], elen[i]);

    return result;
}

void EdgeHash::clear() { CCutil_edgehash_delall(&eh); }

int EdgeHash::get_val(int end1, int end2)
{
    int val = 0;

    if (CCutil_edgehash_find(&eh, end1, end2, &val) == -1)
        return -1;

    return val;
}


}
}
