#ifndef CMR_FIXED64_H
#define CMR_FIXED64_H

extern "C" {
#include <concorde/INCLUDE/bigguy.h>
}

namespace CMR {
namespace util {

class Fixed64 {
public:
    Fixed64() = default;
    Fixed64(int i) : bg(CCbigguy_itobigguy(i)) {}
    Fixed64(double d) : bg(CCbigguy_dtobigguy(d)) {};

    double to_d() const { return CCbigguy_bigguytod(bg); }

    Fixed64 &operator+=(Fixed64 f) { CCbigguy_add(&bg, f.bg); return *this; }
    Fixed64 &operator-=(Fixed64 f) { CCbigguy_sub(&bg, f.bg); return *this; }

    void add_mult(Fixed64 f, int m) { CCbigguy_addmult(&bg, f.bg, m); }

    Fixed64 ceil() const
        { Fixed64 result; result.bg = CCbigguy_ceil(bg); return result; }
    

    bool operator<(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) == -1; }

    bool operator==(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) == 0; }

    bool operator>(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) == 1; }


private:
    CCbigguy bg;
};

inline Fixed64 operator+(Fixed64 a, Fixed64 b) { return a += b; }
inline Fixed64 operator-(Fixed64 a, Fixed64 b) { return a -= b; }

}
}

#endif
