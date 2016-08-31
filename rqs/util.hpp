#ifndef RQS_UTIL_HPP
#define RQS_UTIL_HPP

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <limits>
#include <utility>
#include <vector>


namespace rqs { // Rectangular quantiles sampling namespace

namespace detail {

// Beware: for efficiency the diagonal and RHS terms are modified in place.
template<typename T>
void solve_tridiagonal_system(
    const std::vector<T>& a,
    std::vector<T>& b,
    const std::vector<T>& c,
    std::vector<T>& rhs,
    std::vector<T>& sol )
{
    std::size_t m = a.size();
    
    // eliminate the sub-diagonal
    for (std::size_t i=1; i!=m; ++i) {
        T pivot = a[i]/b[i-1];
        b[i] -= pivot*c[i-1];
        rhs[i] -= pivot*rhs[i-1];
    }
    
    // solve the remaining upper bidiagonal system
    sol[m-1] = rhs[m-1]/b[m-1];
    for (std::size_t i=m-2; ; --i) {
        sol[i] = (rhs[i] - c[i]*sol[i+1])/b[i];
        if (i==0) break;
    }
}        

} // namespace detail



/// Multiple roots finder based on the bisection method.
///
/// The intervals (x0, x1) within which each root should be searched must be
/// provided as std::pair or std::tuple within an STL-compliant sequence.
/// The values of f(x0) and f(x1) at the bounds of each interval should have
/// opposite signs, failing which the returned vector will not contain any root
/// for such interval (the size of the returned vector is the number of roots
/// actually found).
/// Each returned root x is within +/-tol of the actual solution where tol is
/// the specified tolerance.
///
template<typename T, class Func, class InputIt>
std::vector<T> bisect_roots(
    Func f,
    InputIt interval_first,
    InputIt interval_last,
    T tol)
{
    std::vector<T> roots;
    for (InputIt interval=interval_first; interval!=interval_last; ++interval) {
        T x0 = std::get<0>(*interval);
        T x1 = std::get<1>(*interval);
        T xm = 0.5*(x0+x1);
        auto y0 = f(x0);
        auto y1 = f(x1);
        if (y0*y1<0.0) {
            while (std::abs(xm-x0)>tol) {
                auto ym = f(xm);
                if (y0*ym<0.0) {
                    x1 = xm;
                    y1 = ym;
                }
                else {
                    x0 = xm;
                    y0 = ym;
                }
                xm = 0.5*(x0+x1);
            }
            roots.push_back(xm);
        }
    }
    return roots;
}


/// Single root finder based on the bisection method.
///
/// See the documentation for the multiple finder root version. The returned
/// vector contains the solution if found or is empty otherwise.
///
template<typename T, class Func>
std::vector<T> bisect_roots(
    Func f,
    std::pair<T, T> interval,
    T tol )
{
    std::vector<std::pair<T, T>> intervals{interval};
    return bisect_roots<T>(f, intervals.begin(), intervals.end(), tol);
}


/// Compute an approximate qualtile partition using the trapezoidal rule.
///
/// The trapezoidal rule is applied to the given function over a regular grid
/// with 'nb_points' grid points (including outer and inner nodes) to divide
/// interval (x0, x1) into 'nb_partitions' partitions of approximately equal
/// areas.
/// The returned vector is the set of x coordinates.
template<typename T, class Func>
std::vector<T> trapezoidal_rule_approximate_partition(
    Func f,
    T x0,
    T x1,
    std::size_t nb_partitions,
    std::size_t nb_points )
{
    // convenient aliases
    auto& n = nb_points;
    auto& m = nb_partitions;
    
    // compute the curve
    T dx = (x1-x0)/(n-1);
    std::vector<T> x(n);
    std::vector<T> y(n);
    for (std::size_t i=0; i!=(n-1); ++i) {
        x[i] = x0 + i*dx;
        y[i] = f(x[i]);
    }
    x[n-1] = x1;
    y[n-1] = f(x1);
    
    // total area (scaled by 1/dx)
    T s = 0.5*(y[0] + y[n-1]);
    for (std::size_t i=1; i!=(n-2); ++i) {
        s += y[i];
    }
    
    // choose abscissas that split the curve area into equal partitions
    std::vector<T> xp(m+1);
    xp[0] = x0;
    xp[m] = x1;   
    {
        T al = 0.0;
        T ar = 0.5*(y[0] + y[1]);
        std::size_t i=0;
        for (std::size_t j=1; j!=m; ++j) {
            T a = s*(static_cast<T>(j)/static_cast<T>(m));
            while (a>ar) {
                ++i;
                al = ar;
                ar += 0.5*(y[i] + y[i+1]);
            }
            xp[j] = x[i] + (x[i+1]-x[i])*((a-al)/(ar-al));
        }
    }
    
    return xp;
}


/// Compute an approximate qualtile partition using the trapezoidal rule.
///
/// This overload sets the number of grid points to the requested number of
/// partitions.
///
template<typename T, class Func>
std::vector<T> trapezoidal_rule_approximate_partition(
    Func f,
    T x0,
    T x1,
    std::size_t nb_partitions )
{
    return trapezoidal_rule_approximate_partition<T, Func>(
        f, x0, x1, nb_partitions, nb_partitions );
}




/// Upper and lower rectangular quadrature of a function (see Riemann sums).
///
/// A rectangular quadrature similar to that used to compute Riemann sums.
/// The 'yinf' and 'ysup' member variables represent respectively lower
/// and upper quadratures, i.e. rectangles which height is the infimum and
/// supremum of a function over the support of the rectangle.
///
template<typename T>
struct quadrature
{
    std::vector<T> x;
    std::vector<T> yinf;
    std::vector<T> ysup;
};


/// Compute a quadrature with even upper rectangles areas using Newton's method.
///
/// A multivariate Newton method is used to determine the abscissas 'x' such
/// that the areas of all upper quadrature rectangles are equal. For faster
/// convergence it requires a reasonable initial estimate of 'x'.
/// The function derivative and an ordered sequence of the inner function extrema
/// (boundary points excluded) must as well be provided.
/// The tolerance is the maximum relative dispersion of upper rectangle areas,
/// computed as the difference between the largest and smallest area divided by
/// the average area.
/// The maximum number of iterations for the Newtow method may be optionally
/// specified.
///
template<typename T, class Func, class DFunc, class InputIt1, class InputIt2>
quadrature<T> newton_quantile_quadrature(
    Func f,
    DFunc df,
    InputIt1 x_initial_first,
    InputIt1 x_initial_last,
    InputIt2 x_extremum_first,
    InputIt2 x_extremum_last,
    T tol,
    T relax = 1.0,
    std::size_t max_iter = 128 )
{
    // quadrature object and convenient aliases
    quadrature<T> q;
    auto& x    = q.x;
    auto& yinf = q.yinf;
    auto& ysup = q.ysup;
    
    // initialization of the quadrature and extrema vectors
    x.assign(x_initial_first, x_initial_last);
    auto n = x.size() - 1;
    yinf.resize(n);
    ysup.resize(n);
    
    std::vector<std::pair<T,T>> extrema;
    while (x_extremum_first!=x_extremum_last) {
        if ((*x_extremum_first-x.front())*(*x_extremum_first-x.back())<=0.0) {
            extrema.push_back(
                std::pair<T,T>(*x_extremum_first, f(*x_extremum_first)) );
        }
        ++x_extremum_first;
    }
    
    // define the main vectors and pre-compute edge values
    std::vector<T> y(n+1);
    std::vector<T> dx(n-1);
    std::vector<T> dy_dx(n+1);
    std::vector<T> dysup_dxl(n), dysup_dxr(n);
    std::vector<T> minus_s(n-1);
    std::vector<T> ds_dxc(n-1), ds_dxl(n-1), ds_dxr(n-1);
    
    y.front() = f(x.front());
    y.back()  = f(x.back());
    dy_dx.front() = 0.0;
    dy_dx.back()  = 0.0;
    
    std::size_t iter = 0;
    while(true)
    {
        // compute the values at inner points
        for (std::size_t i=1; i!=n; ++i) {
            y[i] = f(x[i]);
            dy_dx[i] = df(x[i]);
        }
        
        // determine the supremum ysup of y in the range (x[i], x[i+1]),
        // the partial derivatives of ysup with respect to x[i] and x[i+1],
        // the minimum and maximum partition areas and the total area
        auto extremum = extrema.begin();
        T max_area = 0.0;
        T min_area = std::numeric_limits<T>::max();
        T sum_area = 0.0;
        for (std::size_t i=0; i!=n; ++i) {
            if (y[i]>y[i+1]) {
                ysup[i] = y[i];
                dysup_dxl[i] = dy_dx[i];
                dysup_dxr[i] = 0.0;
            }
            else {
                ysup[i] = y[i+1];
                dysup_dxl[i] = 0.0;
                dysup_dxr[i] = dy_dx[i+1];
            }
            
            // check if there are extrema within the (x[i], x[i+1]) range
            while ( extremum!=extrema.end() &&
                    (extremum->first-x[i])*(extremum->first-x[i+1])<=0.0 )
            {
                if (extremum->second > ysup[i]) {
                    ysup[i] = extremum->second;
                    dysup_dxl[i] = 0.0;
                    dysup_dxr[i] = 0.0;
                }
                ++extremum;
            }
            
            T area = ysup[i]*std::abs(x[i+1]-x[i]);
            max_area = std::max(area, max_area);
            min_area = std::min(area, min_area);
            sum_area += area;
        }
        
        // check convergence
        if ((max_area-min_area)<tol*(sum_area/n)) {
            // determine the infimum yinf of y in the range (x[i], x[i+1])
            extremum = extrema.begin();
            for (std::size_t i=0; i!=n; ++i) {
                if (y[i]>y[i+1]) {
                    yinf[i] = y[i+1];
                }
                else {
                    yinf[i] = y[i];
                }
                
                // check if there are extrema within the (x[i], x[i+1]) range
                while ( extremum!=extrema.end() &&
                    (extremum->first-x[i])*(extremum->first-x[i+1])<=0.0 )
                {
                    if (extremum->second < yinf[i]) {
                        yinf[i] = extremum->second;
                    }
                    ++extremum;
                }                
            }
            break;
        }
        
        if (++iter>max_iter) {
            q.x.clear();
            q.yinf.clear();
            q.ysup.clear();
            
            break;
        }
        
        // area difference between neigboring rectangles and partial
        // derivatives of s with respect to x[i], x[i+1] and x[i+2]
        for (std::size_t i=0; i!=(n-1); ++i) { 
            minus_s[i] = ysup[i]*(x[i+1]-x[i]) - ysup[i+1]*(x[i+2]-x[i+1]);
            
            ds_dxl[i] = ysup[i] - (x[i+1]-x[i])*dysup_dxl[i];
            ds_dxc[i] = (x[i+2]-x[i+1])*dysup_dxl[i+1]
                      - (x[i+1]-x[i])*dysup_dxr[i]
                      - (ysup[i] + ysup[i+1]);
            ds_dxr[i] = ysup[i+1] + (x[i+2]-x[i+1])*dysup_dxr[i+1];
        }
        
        // solve the tri-diagonal system S + (dS/dX)*dX = 0 with:
        //         | ds0/dx1 ds0/dx2    0     ...                    0     |
        //         | ds1/dx1 ds1/dx2 ds1/dx3    0     ...            0     |
        // dS/dX = |    0    ds2/dx2 ds2/dx3 ds2/dx4    0     ...    0     |
        //         |                       ...                             |
        //         |    0     ...     0    ds(n-2)/dx(n-2) ds(n-2)/dx(n-2) |
        //
        //
        // and:
        //      | dx1     |         | minus_s0     |
        // dX = | ...     |    -S = | ...    |
        //      | dx(n-1) |         | minus_s(n-2) |
        
        detail::solve_tridiagonal_system<T>(
            ds_dxl, ds_dxc, ds_dxr, minus_s, dx );
        
        // for the sake of convergence stability, updated positions are
        // constrained within the bounds set by former neighbors positions
        {
            T x0 = x[0];
            for (std::size_t i=1; i!=n; ++i) {
                std::pair<T,T> x_range = std::minmax(x0, x[i+1]);
                x0 = x[i];
                x[i] = std::max(x[i] + relax*dx[i-1], x_range.first);
                x[i] = std::min(x[i], x_range.second);
            }
        }
    }
    
    // voila
    return q;
}


/// Tail of a Weibull distribution generated by inversion sampling.
///
/// Generates the tail of a Weibull distribution such that:
///
///  f(x;a,b,c) = s*((x-c)/b)*exp[-((x-c)/b)^a] for x/b > x0/b
///  f(x;a,b,c) = 0 otherwise
///
/// where 'a' is strictly positive. The scale parameter 'b' may be positive for
/// a tail extending to +inf and negative for a tail extending to -inf.
/// The (positive) normalization constant 's' need not be specified.
///
template<typename T>
class weibull_tail_distribution
{
public:
    using result_type = T;
    using param_type = weibull_tail_distribution<T>;
    
    
    weibull_tail_distribution(T a=1.0, T b=1.0, T c=0.0, T x0=0.0) :
        inv_a_(1.0/a), b_(b), c_(c), x0_(x0), alpha_(std::pow((x0-c)/b, a))
    {}
    
    
    template<class RngType>
    result_type operator()(RngType& g)
    {
        constexpr auto digits = std::numeric_limits<T>::digits;
        
        T r = std::generate_canonical<T,digits>(g);
        return c_ + b_*std::pow(alpha_ - std::log(1.0-r), inv_a_);
    }
    
    
    void reset()
    {}
    
    
    param_type param() const
    {
        return *this;
    }
    
    
    void param(const param_type& params)
    {
        *this = params;
    }
    
    
    result_type min() const
    {
        T minus_inf = std::numeric_limits<T>::is_iec559 ?
                      -std::numeric_limits<T>::infinity()
                    : std::numeric_limits<T>::min();
        return b_<0 ? minus_inf : x0_;
    }
    
    
    result_type max() const
    {
        T plus_inf = std::numeric_limits<T>::is_iec559 ?
                     std::numeric_limits<T>::infinity()
                   : std::numeric_limits<T>::min();
        return b_>0 ? plus_inf : x0_;
    }
    
    
    result_type a() const
    {
        return 1.0/inv_a_;
    }
    
    
    result_type b() const
    {
        return b_;
    }
    
    
    result_type c() const
    {
        return c_;
    }
    
    
private:
    T inv_a_;
    T b_;
    T c_;
    T x0_;
    T alpha_;
};




// weibull_pdf<T>
// Weibull probability density function with optional weighting.
//
// The Weibull pdf is defined as:
//
//  f(x;a,b,c) = w*a/|b|*((x-c)/b)*exp[-((x-c)/b)^a] for x/b > x0/b
//  f(x;a,b,c) = 0 otherwise
//
// where 'a' is strictly positive, 'b' may be negative and 'w' is an optional
// weighting factor.

template<typename T>
class weibull_pdf
{
public:
    using result_type = T;
    
    
    weibull_pdf(T a=1.0, T b=1.0, T c=0.0, T w=1.0) :
        a_(a), inv_b_(1.0/b), c_(c), s_(w*std::abs(a/b))
    {}
    
    
    result_type operator()(T x) const
    {
        T y = (x - c_)*inv_b_;
        if (y<0) return T(0);
        T z = std::pow(y, a_-1.0);
        
        return s_*z*std::exp(-y*z);
    }
    
    result_type total_weight() const
    {
        return s_/std::abs(a_*inv_b_);
    }
    
    result_type tail_weight(T x0) const
    {
        T z0 = std::pow((x0 - c_)*inv_b_, a_);
        
        return s_*std::exp(-z0)/(a_*inv_b_);
    }
    
    
private:
    T a_;
    T inv_b_;
    T c_;
    T s_;
};



} // namespace rqs


#endif // RQS_UTIL_HPP

