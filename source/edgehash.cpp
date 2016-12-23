#include "edgehash.hpp"

namespace CMR {
namespace util {

EdgeHash::EdgeHash(int size) try :
    hash_handle(new CCutil_edgehash, CCutil_edgehash_free)
{
    if (CCutil_edgehash_init(hash_handle.get(), factor * size)) {
        hash_handle.reset(nullptr);
        throw std::runtime_error("CCutil_edgehash_init failed.");
    }
} catch (const std::exception &e) {
    std::cerr << e.what() << "\n";
    throw std::runtime_error("EdgeHash constructor failed.");
}

}
}
