#ifndef ZIGGURAT_NORMAL_HPP
#define ZIGGURAT_NORMAL_HPP

#include <cmath>
#include <cstddef>
#include <vector>

#include <rqs/generate_random.hpp>

#include "tail_dist.hpp"




// An arbitrary precision Ziggurat implementation.
//
// This is a relatively straightforward adaptation of the original algorithm
// (Marsaglia and Tsang, 2000) which is generic over the RNG and the floating
// point type, with a user-specified bit width W for the underlying random
// integer.
//
// Bit width W should not be greater than that of the longest unsigned integer
// available on the platform and must be able to hold more than the 7 bits
// required by the Ziggurat algorithm to generate the table index.
//
// The RNG passed as an argument to operator() should ideally have min=0 and
// max=2^N-1 (arbitrary N), but min=1 and/or max=2^N-2 are tolerated.
//
// If the raw random numbers generate by the RNG have less than W bits, 2 or
// more numbers are generated, as appropriate.
//
// Apart from the obvious, the main differences with the original algorithm are:
//  - a bug in the original algorithm was fixed; this algorithm applied abs()
//    to a signed integer which could occasionally have the smallest possible
//    negative value, the result of which is undefined since the absolute
//    value of the smallest negative integer cannot be represented by a positive
//    value and in fact abs() will then typically return a _negative_ integer in
//    such case, causing a bug in the acceptance test; surprisingly the fix made
//    the code slightly faster than the original algorithm on gcc both and
//    clang,
//  - the correlation flaw of the original algorithm is fixed (Doornik, 2005),
//  - both the sign and the table index are determined from the upper RNG bits
//    because many RNGs have worse quality lower bits (e.g. xorshift family
//    RNGs).
//
template<typename RealType, int W>
class ZigguratNormalDistribution
{
private:
    using IntType = typename rqs::integer_traits<W>::int_fast_t;
    using UIntType = typename rqs::integer_traits<W>::uint_fast_t;

public:
    ZigguratNormalDistribution() : tail_dist_(3.442619855899)
    {
        const std::size_t n = 128;
        const RealType xt = 3.442619855899;
        const RealType weight = 9.91256303526217e-3;
        const IntType scale = static_cast<IntType>(UIntType(1) << (W-8));

        w_.resize(n);
        k_.resize(n);
        f_.resize(n);

        // Build the tables.
        f_[0] = 1;
        f_[n-1] = std::exp(-RealType(0.5)*xt*xt);
        w_[0] = weight/(f_[n-1]*scale);
        w_[n-1] = xt/scale;
        k_[0] = static_cast<IntType>(xt/w_[0]);
        k_[1] = 0;
        
        RealType xmax_prev = xt;
        for(std::size_t i=n-2; i!=0; --i)
        {
            RealType xmax = std::sqrt(-RealType(2)*std::log(weight/xmax_prev + f_[i+1]));
            k_[i+1] = static_cast<IntType>((xmax/xmax_prev)*scale);
            w_[i] = xmax/scale;
            f_[i] = std::exp(-RealType(0.5)*xmax*xmax);
            xmax_prev = xmax;
        }
    }
    

    template<class G>
    RealType operator()(G& g)
    {
        using std::abs;
        using std::exp;

        while (true) {
            auto r = rqs::generate_random_integer<UIntType, W>(g);
            // Integer u is made of bits 0:(W-8); unlike the original ziggurat,
            // we avoid relying on casting to obtain a potentially negative
            // signed integer from an unsigned integer since the C standard does
            // not guarantee the portability of overflowing casts. Instead, we
            // cast an unsigned integer that is small enough to be representable
            // as a signed integer and then apply a complement to obtain a
            // potentially negative integer.
            // Since the signed integer can hold at least W-1 digits while the
            // actual value has at most W-8 digits, using abs() on a negative
            // value is guarranteed to give a representable positive value
            // (this fixes a bug in the original ziggurat where abs() could be
            // potentially applied to the smallest signed integer, with
            // undefined consequences).
            // With both gcc and clang this appears to be even slightly faster
            // than the original ziggurat code.
            constexpr UIntType m_mask = (UIntType(1) << (W - 7)) - 1;
            constexpr IntType complement =
                static_cast<IntType>((UIntType(1) << (W - 8)) - 1);
            IntType u = complement - IntType(r & m_mask);
            // The table index is made of bits (W-7):(W-1).
            std::size_t i = r >> (W - 7);
            
            if (abs(u)<k_[i])
                return u*w_[i];

            if (i==0)
                return u>0 ? tail_dist_(g) : -tail_dist_(g);

            RealType x = u*w_[i];
            RealType v = rqs::generate_random_real<RealType, W>(g);
            if (f_[i] + v*(f_[i-1] - f_[i]) <= exp(RealType(-0.5)*x*x))
                return x;
        }
    }

private:
    std::vector<RealType> w_;
    std::vector<IntType> k_;
    std::vector<RealType> f_;
    NormalTailDistribution<RealType, W> tail_dist_;
};

#endif // ZIGGURAT_NORMAL_HPP

