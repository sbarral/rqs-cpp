#include <iostream>
#include <random>
#include <vector>

#define NONIUS_RUNNER
#include <nonius_single.h++>

#include "rqs_normal.hpp"
#include "ziggurat_normal.hpp"


NONIUS_BENCHMARK("standard library implementation (64-bit)", [](nonius::chronometer meter) {
    using GeneratorType = std::mt19937;

    std::size_t nb_iter = 10000;
    GeneratorType rng;
    std::normal_distribution<double> normal_dist;

    meter.measure([&] {
        double s = 0.0;
        for (std::size_t i=0; i!=nb_iter; ++i) {
            s += normal_dist(rng);
        }
        return s;
    });
})

NONIUS_BENCHMARK("ziggurat (32-bit)", [](nonius::chronometer meter) {
    using GeneratorType = std::mt19937;

    std::size_t nb_iter = 10000;
    GeneratorType rng;
    ZigguratNormalDistribution<double, 32> normal_dist;

    meter.measure([&] {
        double s = 0.0;
        for (std::size_t i=0; i!=nb_iter; ++i) {
            s += normal_dist(rng);
        }
        return s;
    });
})

NONIUS_BENCHMARK("ziggurat (64-bit)", [](nonius::chronometer meter) {
    using GeneratorType = std::mt19937_64;

    std::size_t nb_iter = 10000;
    GeneratorType rng;
    ZigguratNormalDistribution<double, 64> normal_dist;

    meter.measure([&] {
        double s = 0.0;
        for (std::size_t i=0; i!=nb_iter; ++i) {
            s += normal_dist(rng);
        }
        return s;
    });
})
NONIUS_BENCHMARK("RQS (32-bit)", [](nonius::chronometer meter) {
    using GeneratorType = std::mt19937;

    std::size_t nb_iter = 10000;
    GeneratorType rng;
    RqsNormalDistribution<double, 32, 7> normal_dist;

    meter.measure([&] {
        double s = 0.0;
        for (std::size_t i=0; i!=nb_iter; ++i) {
            s += normal_dist(rng);
        }
        return s;
    });
})

NONIUS_BENCHMARK("RQS (64-bit)", [](nonius::chronometer meter) {
    using GeneratorType = std::mt19937_64;

    std::size_t nb_iter = 10000;
    GeneratorType rng;
    RqsNormalDistribution<double, 64, 7> normal_dist;

    meter.measure([&] {
        double s = 0.0;
        for (std::size_t i=0; i!=nb_iter; ++i) {
            s += normal_dist(rng);
        }
        return s;
    });
})

