#ifndef CMR_FIXED64_H
#define CMR_FIXED64_H

extern "C" {
#include <concorde/INCLUDE/bigguy.h>
}

#include <iostream>

namespace CMR {
namespace util {

/// A 64-bit fixed precision number type "implemented" as a wrapper to CCbigguy.
class Fixed64 {
public:
    Fixed64() = default; //!< Default constructor.
    
    Fixed64(int i) :
        bg(CCbigguy_itobigguy(i)) {} //!< Construct from integer.
    
    Fixed64(double d) :
        bg(CCbigguy_dtobigguy(d)) {}; //!< Construct from double. 

    
    double to_d() const
        { return CCbigguy_bigguytod(bg); } //!< Convert to double.

    Fixed64 &operator+=(Fixed64 f)
        { CCbigguy_add(&bg, f.bg); return *this; } //!< Plus increment.
    
    Fixed64 &operator-=(Fixed64 f)
        { CCbigguy_sub(&bg, f.bg); return *this; } //!< Minus decrement.

    friend void add_mult(Fixed64 &f, const Fixed64 &g, int m);
    
    /// Integer ceiling.
    Fixed64 ceil() const
        { Fixed64 result; result.bg = CCbigguy_ceil(bg); return result; }
    

    bool operator<(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) == -1; } //!< Less than operator.

    bool operator==(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) == 0; } //!< Equality operator.

    bool operator!=(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) != 0; }

    bool operator>(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) == 1; } //!< Greater than operator.

    bool operator<=(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) <= 0; }

    bool operator>=(const Fixed64 &f) const
        { return CCbigguy_cmp(bg, f.bg) >= 0; }


private:
    CCbigguy bg; //!< The underlying Concorde structure.
};

inline Fixed64 operator+(Fixed64 a, Fixed64 b) { return a += b; }
inline Fixed64 operator-(Fixed64 a, Fixed64 b) { return a -= b; }

inline void add_mult(Fixed64 &f, const Fixed64 &g, int m)
{
    CCbigguy_addmult(&f.bg, g.bg, m);
}

inline void add_mult(double &d, const double &g, int m)
{
    d += m * g;
}

inline std::ostream &operator<<(std::ostream &os, const Fixed64 &f)
{
    os << (f.to_d());
    return os;
}

}
}

#endif
