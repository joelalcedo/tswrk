#include "ts/time_series.hpp"
#include <chrono>

int main() {
    ts::TimeSeries s;

    using clock = std::chrono::system_clock;
    auto t0 = clock::now();

    for (int i = 0; i < 15; ++i) {
        s.push_back(t0 + std::chrono::hours(24 * i), 2.04 + 0.2 * i);
    }

    ts::head(s, 10);
    ts::tail(s, 5);
    // ts::print(s);

    return 0;
}
