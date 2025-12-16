#pragma once

#include <chrono>
#include <stdexcept>
#include <utility>
#include <vector>

class TimeSeries {
public:
    using time_point = std::chrono::system_clock::time_point;
    using point = std::pair<time_point, double>;

    TimeSeries() = default;

    void push_back(time_point t, double value) {
        data_.emplace_back(t, value);
    }

    std::size_t size() const noexcept { return data_.size(); }
    bool empty() const noexcept { return data_.empty(); }

    const point& at(std::size_t i) const { return data_.at(i); }
    const point& operator[](std::size_t i) const { return data_[i]; }


    time_point time_at(std::size_t i) const { return data_.at(i).first; }
    double value_at(std::size_t i) const { return data_.at(i).second; }

    auto begin() const noexcept { return data_.begin(); }
    auto end() const noexcept { return data_.end(); }

private:
    std::vector<point> data_;
};

