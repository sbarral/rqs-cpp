# RQS: rectangular quantiles sampling

C++ implementation of a blazingly fast algorithm to sample arbitrary continuous
probability distributions solely from their probability density function (PDF)
and using only one random number per sample.


## Why?

The RQS algorithm has a number of advantages over other methods, including:

* **speed**: it is comparable in speed to the 
[ziggurat](https://en.wikipedia.org/wiki/Ziggurat_algorithm) algorithm
[[1]](#references) and thus typically many times faster than
[inversion sampling](https://en.wikipedia.org/wiki/Inverse_transform_sampling),

* **quality**: the entropy of the generated distribution is
[noticeably higher](#quality)
than that of the ziggurat and quite close to that of inversion sampling,

* **convenience**: there is no need to know the cumulative distribution function
(CDF) and even less its inverse; the lookup tables can be automatically
generated with a multivariate newton method provided within the library,
requiring only the PDF and its derivative (but even the PDF derivative could be
omitted with the use of e.g. a multivariate secant method),

* **versatility**: unlike the ziggurat, it can be applied to a wide range of
distributions including non-monotonic or asymmetric distributions; also,
the *x* position beyond which the tail is sampled with a fallback algorithm
can be arbitrarily adjusted without consideration for the statistical weight
of the tail,

* **thread-safety**: unlike methods such as the Box-Muller transform
or the polar method, it is inherently stateless and can thus be safely and
efficiently shared between threads (assuming of course thread-local RNGs).


## How does it work?

This algorithm was born from necessity after realizing that there did not seem
to exist convenient and fast methods to generate arbitrary continuous
distributions for which the inverse of the CDF has no analytical expression
or cannot otherwise be computed efficiently.
The ziggurat algorithm came close, but it cannot cope with non-monotonic
distributions and is quite inconvenient as a ready-to-use algorithm due to its
complex set-up and its constraints on the tail area (the tail is included in the
bottom layer which is required to have the same total area as other layers).

The underlying idea of RQS is simple: instead of sampling horizontal 
rectangular quantiles like in the ziggurat algorithm, vertical rectangles are
used instead.
This idea is nothing exotic in itself and I was therefore not surprised to
find that it had already been investigated [[2]](#references). What was
apparently overlooked, however, is that just like the ziggurat algorithm it
can be very efficiently implemented using only one random number per sample.

For simplicity, let us consider a non-symmetric distribution with compact
support.
The algorithm starts by computing majorizing and minorizing functions for the
PDF using vertical rectangles. This is quite similar to the procedure used to
compute upper and lower
[Riemann sums](https://en.wikipedia.org/wiki/Riemann_sum),
but instead of splitting the *x* interval evenly, it is divided so that
all rectangles of the majorizing functions have equal areas. This information
is tabulated as part of pre-processing, typically in a 128 or 256-elements
table.

Sampling is based on an optimized rejection algorithm that uses the
above-described majorizing function. Ignoring for now some optimizations of the
actual implementation, random variates are generated as follows:

1. a *W*-bit random integer is generated,

2. the upper rectangular quantile to be sample is determined by extracting *N*
bits (typically *N*=7 or *N*=8) from the random integer,

3. the remaining bits (*W*-*N*) are used to generate a real number *y* between 0
and *H* where *H* is the height of the upper rectangle,

4. if *y* is larger than the height *h* of the lower rectangle then jump to 6
(rejection), otherwise proceed to 5 (acceptance),

5. \[*rectangle sampling: here is the trick*\]: we are left with *y*, a nice
number uniformly distributed between in \[0, *h*) with almost the same entropy
as the untested *y* —after all we have only lost *log2(H/h)* bits— so why throw
it to the garbage?
Instead, we transform linearly *y* from \[0, *h*) to the interval
\[*xᵢ* and *xᵢ₊₁*) and return the number *x* thus produced.

6. \[*wedge sampling: this happens rarely*\] generate another random number
*x* within the interval \[*xᵢ* and *xᵢ₊₁* ) of the current rectangle and check
whether *y* is within the PDF envelope; if yes, return *x*, otherwise jump
back to 1.

A legitimate question would be whether we have not lost the *N* bits of
precision in generating the quantile index.
It is useful to note, however, that the generation of *x* actually requires
less precision precisely because once the quantile has been determined, we
have in effect already narrowed down the interval to \[*xᵢ* and *xᵢ₊₁*), and
since there are 2*ᴺ* such intervals, we need this much less precision.
In fact, the method can be shown to be nothing else than a fast inversion
sampling method for the majorizing function with the addition of rejection
sampling of the wedges, so the precision loss is in theory relatively moderate
and mostly imputable to wedge rejection.

For tailed distributions, the algorithm admits an externally supplied (fallback)
tail distribution, just like the ziggurat algorithm, but does not constrain
the *x* threshold beyond which the fallback needs to be used. It can also use
instead a majorizing PDF and apply rejection sampling in the tail, in which case
a truncated Weibull distribution often makes a good majorizing function with a
quite fast inverse transform generation.
For this reason, the truncated Weibull distribution is provided as part of this
library.

To include a tail, step 4 above is modified to generate a real number *y*
between 0 and *H/\[1-P(tail)\]* where *P(tail)* is the tail sampling probability
(i.e. the area of the tail divided by the total area of the majorizing function,
tail included).
If *y* is greater than *H*, then the tail is sampled.

Note that for the sake of efficiency, in the actual implementation *y* is an
integer and the various rejection thresholds are appropriately scaled so as
to only use integer arithmetic until *x* is to be computed.


## Quality

In the spirit of the C++ standard library, this library is RNG-agnostic so
potential quality issued related to the RNG itself will not be investigated.
Also, since this is a rejection method, goodness of fit is assumed without
demonstration.

The main question becomes therefore: how much of the quality of the original
RNG does it actually retain? Or alternatively, how "uniformly" (in a loose
meaning) would the numbers be distributed if we used an ideal RNG?
From this perspective, inversion sampling can be considered the reference
method: baring rounding issues and provided that the floating point mantissa
has more digits that the random number, one can always map a number generated
by inversion sampling back to the original random number using the CDF.

One way to crudely assess the quality is by sampling the distribution with a
high quality RNG, using many more samples than the range 2*ᵂ* of the RNG. If we
set for now our ambitions and *W* low, at say 16 bits, this is fairly easy to
do.

Without further ado, this is the result of the collection of 10⁹ normal variates
in 1001 bins, generated with the ziggurat and the RQS algorithms with a 16-bit RNG
(actually, a truncated mt19937).

![Normal variates from ziggurat and RQS](https://sbarral.github.io/rqs/fig/fig_16bit_comparison.png)

Yes, quite a difference.

It is important to note that these are the converged numerical PDFs: all
irregularities are inherent to the algorithms used and are in no way related
to the RNG: adding more samples or using another (good) RNG would not change
these plots in any visible manner. A testimony to convergence can incidentally
be found in the symmetry of the scatter plots.

But can we quantify this, especially at higher precision *W*? To some extent,
yes. A test that can point to poor distribution of the variates at small scales
without requiring too much computing power is the Knuth collision test
[[3]](#references) where *n* balls are randomly thrown into 2*ᵈ* urns.
Even with *n*<2*ᵈ*, the number of collisions (i.e. balls that get into an
already filled urn) can be sufficient to form the basis for a statistical
test that assesses uniformity defects.

Following Doornik [[4]](#references) we first generate normal variates and then
use the normal CDF to transform them back to a random number that *should* be
uniformly distributed in \[0, 1).
The Knuth collision test is then performed for several values of *d*, keeping
*n* smaller than 2*ᵈ* by a (somewhat arbitrary) factor of 256.
The RNG is a mt19937 with *W*=32 bits and all floating point computations
are performed in double precision.
The ziggurat and RQS algorithm both use 128-entry tables.
The inversion sampling test used for comparison is actually an ideal case: it
is simulated by simply transforming the RNG integer directly to a uniform number
in \[0, 1), which is akin to performing inversion sampling and then transforming
back the variate with the normal CDF using infinite floating point precision.

So here are the results, with each test repeated 10 times for each value of *d*.
See the [source](benchmark/collision.cpp) for details (note that
the Boost special functions library is required).

![p-values for d between 26 and 33](https://sbarral.github.io/rqs/fig/fig_pvalue.png)

Without going into details, if the variates are well distributed then the
10 p-values (scatter plot) should be uniformly distributed between 0 and 1 and
their mean value (solid line) should be close to 0.5. When all p-values go
low (the traditional threshold is 5%), this is a sure sign that the test has
failed.

So what does this plot say?

Expectedly all algorithms fail at *d*=33 since with *W*=32 bits of information
it is obviously impossible to evenly distribute balls into 2*³³* urns.
Simulated inversion sampling succeeds at *d*=32, but this was also expected
since this ideal inversion sampling test is de-facto a measure of the quality
of mt19937 (which is known to be good).

The quality of the RQS algorithm is very satisfactory with a loss of only 2
effective bits compared to ideal inversion sampling. The ziggurat algorithm is
clearly worse with a 5-bit effective loss compared to ideal inversion sampling,
thus confirming the hint given by our former visual comparison of distributions
generated with *W*=16.


## Speed

The ziggurat was of course used as the reference algorithm for speed since it
is widely considered the fastest method to generate normal variates of
reasonable quality.
In order to ensure a fair comparison, its
[implementation](benchmark/ziggurat_normal.hpp)
was heavily optimized and appears to be even a tiny bit faster than the
original (assuming the same RNG of course), with the advantage that it is
RNG-agnostic and generic over both the floating point type and *W*.
An (apparently unknown) bug and a known issue [[4]](#references) were also
fixed along the way; see comments in the source.

The RQS implementation reflects most of the feature of the ziggurat
implementation but it is additionally generic over the number of bits *N* of the
table index, unlike the ziggurat which is hardcoded with *N*=7.
In order to compare apples to apples, all tests are performed with *N*=7 for the
RQS algorithm as well.

The speed benchmark was run on a i5-3230M (yes I know — I accept donations).
It sums 10000 normal variates with either `std::mt19937` (for *W*=32) or
`std::mt19937_64` (for *W*=64) as RNG.
Each of the reported timing is an average over 1000 tests.
Please feel free to contribute your timings using the
[benchmarking code](benchmark/timing.cpp). Note that this requires the
[nonius](https://nonius.io) microbenchmarking library which in turn needs a
couple of the Boost libraries.

|                                   | gcc 4.9.2 (-O2)  | clang 3.5.0 (-O2) |
| --------------------------------- | ---------------- | ----------------- |
| Ziggurat, W=32                    |  68.7 µs         |  108.0 µs         |
| RQS, W=32                         |  74.8 µs   (+9%) |  119.7 µs  (+11%) |
| Ziggurat, W=64                    |  79.5 µs         |  102.9 µs         |
| RQS, W=64                         |  81.0 µs   (+2%) |  113.4 µs  (+10%) |
| libstdc++ (polar transform), W=64 | 417.9 µs (+426%) |  939.0 µs (+812%) |


Despite a few oddities like clang taking longer to compute with a 32-bit RNG
than with a 64 bit RNG, the picture is very consistent: the difference between
the ziggurat and the RQS is barely noticeable with only a 2% to 11% penalty for
the RQS, compared to an abysmal +426% and +812% for the native libstdc++
implementation of the normal distribution which uses the polar transform.


## Todo

The library is not quite complete yet as I am finishing up the design of a
generic, ready-to-use RQS templated class with specializations for
symmetric and/or tailed distributions, with and without tail rejection sampling.

In the meantime you may have a look at the
[sample implementation](benchmark/rqs_normal.hpp) used in the benchmark above.
Some doc and examples will be coming soon too, stay tuned... And hopefully a
Rust implementation someday.


## References

\[1\] G. Marsaglia and W. W. Tsang,
["The ziggurat method for generating random variables"]
(https://www.jstatsoft.org/article/view/v005i08),
*Journal of Statistical Software 5*, 1-7 (2000).

\[2\] R. Zhang and L. M. Leemis,
["Rectangles algorithm for generating normal variates"]
(http://dx.doi.org/10.1002/nav.21474),
*Naval Research Logistics 59(1)*, 52-57 (2012).

\[3\] D. E. Knuth, "The art of computer programming",
*Vol. 2*, 3rd ed., Addison-Wesley (1997).

\[4\] J. A. Doornik
["An improved ziggurat method to generate normal random samples"]
(http://www.doornik.com/research/ziggurat.pdf),
Technical note, University of Oxford (2005).

## License

This software is licensed under the [Apache License, Version 2.0](LICENSE-APACHE), the
[MIT license](LICENSE-MIT) or the [Boost license](LICENSE-BOOST) at your option.

Copyright (c) 2016 Serge Barral.

