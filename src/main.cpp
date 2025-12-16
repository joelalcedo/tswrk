#include "ts/frame.hpp"
#include <chrono>
#include <cmath>

int main() {
    ts::Frame df;

    using clock = std::chrono::system_clock;
    auto t0 = clock::now();

    for (int i = 0; i < 45; ++i) df.index.push_back(t0 + std::chrono::hours(24 * i));

    std::vector<double> px, vol;
    for (int i = 0; i < 45; ++i) {
        px.push_back(100.0 + std::sin(i / 3.0));
        vol.push_back((i % 7) + 1.0);
    }

    df.add_col("px", px);
    df.add_col("vol", vol);

    auto out = df
      .arrange_time()
      .filter([](const ts::Frame::Row& r){ return r("vol") >= 3.0; })
      .mutate("px2", [](const ts::Frame::Row& r){ return r("px") * r("px"); })
      .lag("px", "px_l1", 1)
      .group_by_time(ts::Frame::freq::month)
      .mean("px", "px_monthly_mean");

    out.tail(10);
}
