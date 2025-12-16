#include "ts/time_series.hpp"
#include <chrono>
#include <cmath>

int main() {
    ts::TimeSeries s;

    using clock = std::chrono::system_clock;
    auto t0 = clock::now();

    for (int i = 0; i < 45; ++i) {
        s.push_back(t0 + std::chrono::hours(24 * i), 100.0 + std::sin(i / 3.0));
    }

    auto out = s
        .lag(1)
        .mutate([](double x){ return x * 1.01; })
        .group_by_month().mean();

    out.print();
    return 0;
}
