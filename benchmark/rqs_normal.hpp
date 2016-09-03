#ifndef RQS_NORMAL_HPP
#define RQS_NORMAL_HPP

#include <array>
#include <cstddef>
#include <cmath>
#include <limits>
#include <vector>

#include <rqs/util.hpp>
#include <rqs/generate_random.hpp>

#include "tail_dist.hpp"


// Non-normalized pdf of the normal distribution.
template<typename RealType>
struct Pdf {
    RealType operator()(RealType x) {
        return std::exp(-RealType(0.5)*x*x);
    }
};


// Derivative of the non-normalized pdf.
template<typename RealType>
struct DPdf {
    RealType operator()(RealType x) {
        return -x*std::exp(-RealType(0.5)*x*x);
    }
};


// RQS-based central normal distribution.
template<typename RealType, int W, int N>
class RqsNormalDistribution
{
private:
    using UIntType = typename rqs::integer_traits<W>::uint_least_t;

    struct Datum
    {
        RealType x;
        UIntType scaled_fratio;
        RealType scaled_finf;
        RealType scaled_dx;
    };


public:
    RqsNormalDistribution() :
        tail_dist_(1)
    {

        // Set the table size.
        const std::size_t n = std::size_t(1) << N;
        data_.resize(n);
        
        // For high precision (W large), the position of the tail can be chosen
        // rather freely; for low W values, though, the tail position should be
        // chosen such that the area of the tail relatively to the whole area
        // sampled (upper rectangles + tail) is a multiple of 1/2^(W-N-1),
        // otherwise the tail sampling probability will not be accurate.
        // The magical values that are closest to the empirical optimum of the
        // tail position (around 3.25) are tabulated for the common cases
        // N=7 and N=8.
        RealType xtail;
        if (N==7 && W>=11) {
            RealType magic_xtail[] =
                { 1.532095304, 1.859950459, 2.150455371, 2.413614185,
                  2.655703474, 2.880953316, 3.092363645, 3.292145211 };
            xtail = magic_xtail[std::min(W - 11, 7)];
        }
        else if (N==8 && W>=12) {
            RealType magic_xtail[] =
                { 1.533103263, 1.861331463, 2.152146391, 2.415553089,
                  2.657829951, 2.883210552, 3.094702254, 3.294526271 };
            xtail = magic_xtail[std::min(W - 12, 7)];
        }
        else {
            // let's hope W is large...
            xtail = 3.25;
        }
        tail_dist_ = NormalTailDistribution<RealType, W>(xtail);

        // The tail is sampled with Marsaglia's algorithm.
        const RealType sqrt_pi_over_two = 1.2533141373155001;
        RealType tail_area = sqrt_pi_over_two*
            std::erfc(xtail/sqrt(RealType(2.0)));

        // Compute the quantiles.
        const RealType rel_tol = std::numeric_limits<RealType>::epsilon()*1e4;

        auto x_guess = rqs::trapezoidal_rule_equipartition(
            Pdf<RealType>(), RealType(0), xtail, n);

        auto p = rqs::newton_monotonic_function_partition(
            Pdf<RealType>(), DPdf<RealType>(),
            x_guess.begin(), x_guess.end(),
            rel_tol);
        
        // Compute the tail switch, i.e. an integer threshold such that when
        // drawing a random integer r, the probability:
        //  P(r>=switch)
        // expresses the probability to sample the tail.
        RealType upper_quadrature_area = 0.0;
        for (std::size_t i=0; i!=n; ++i) {
            upper_quadrature_area += (p.x[i+1] - p.x[i])*p.fsup[i];
        }
        tail_switch_ = static_cast<UIntType>(
            std::round(RealType(UIntType(1) << (W - N - 1)) *
            (upper_quadrature_area/(tail_area + upper_quadrature_area))));
        
        // Compute the tables.
        for (std::size_t i=0; i!=n; ++i) {
            data_[i].x = p.x[i];
            data_[i].scaled_fratio = static_cast<UIntType>(
                (p.finf[i]/p.fsup[i])*tail_switch_);
            data_[i].scaled_finf = p.fsup[i]/tail_switch_;
            data_[i].scaled_dx = (p.x[i+1] - p.x[i])/data_[i].scaled_fratio;
        }
    }
    

    template<class G>
    RealType operator()(G& g)
    {
        using std::exp;

        while (true)
        {
            // Generate a table index, a sign bit and a positive value from a
            // single random number.
            auto r = rqs::generate_random_integer<UIntType, W>(g);
            // Mantissa is made of bits 0:(W-N-2).
            constexpr UIntType m_mask = (UIntType(1) << (W - N - 1)) - 1;
            UIntType u = r & m_mask;
            // The table index is made of bits (W-N-1):(W-2).
            constexpr std::size_t i_mask = (std::size_t(1) << N) - 1;
            auto i = std::size_t(r >> (W - N - 1)) & i_mask;
            // Sign is bit (W-1).
            int s = r >> (W - 1) ? 1 : -1;
            
            const Datum& d = data_[i];
            // Note that the following test will also fail if 'u' is greater
            // than the tail switch value since all fratio values are lower
            // than the switch value.
            if (u<=d.scaled_fratio)
                return s*d.x + s*d.scaled_dx*u;
            
            // Should the tail be sampled?
            if (u>=tail_switch_)
                return s*tail_dist_(g);

            // Otherwise it is a wedge, test y<f(x) for rejection sampling.
            RealType v = rqs::generate_random_real<RealType, W>(g); // v in [0,1)
            RealType x = d.x + v*(d.scaled_dx*d.scaled_fratio);
            if ((u*d.scaled_finf) < exp(RealType(-0.5)*x*x))
                return s*x;
        }
    }


private:
     std::vector<Datum> data_;
     UIntType tail_switch_;
     NormalTailDistribution<RealType, W> tail_dist_;
};

#endif // RQS_NORMAL_HPP
