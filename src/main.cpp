#include "ts/time_series.hpp"

#include <chrono>
#include <iostream>

int main() {
    ts::TimeSeries s;

    using clock = std::chrono::system_clock;
    auto t0 = clock::now();

    s.push_back(t0, 100.0);
    s.push_back(t0 + std::chrono::hours(24), 101.5);

    std::cout << "size=" << s.size() << "\n";
    std::cout << "first value=" << s.value_at(0) << "\n";
    std::cout << "second value=" << s.value_at(1) << "\n";
    return 0;
}