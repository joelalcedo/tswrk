#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace ts {

class TimeSeries {
public:
    using time_point = std::chrono::system_clock::time_point;
    using point = std::pair<time_point, double>;

    void push_back(time_point t, double value) { data_.emplace_back(t, value); }

    std::size_t size() const noexcept { return data_.size(); }
    bool empty() const noexcept { return data_.empty(); }

    const point& at(std::size_t i) const { return data_.at(i); }

    time_point time_at(std::size_t i) const { return data_.at(i).first; }
    double value_at(std::size_t i) const { return data_.at(i).second; }

    // Prints first n rows (like head(x, n))
    void head(std::size_t n = 10, std::ostream& os = std::cout) const {
        const std::size_t end = std::min(n, size());
        print_range(0, end, os);
    }

    // Prints last n rows (like tail(x, n))
    void tail(std::size_t n = 10, std::ostream& os = std::cout) const {
        const std::size_t sz = size();
        const std::size_t start = (sz > n) ? (sz - n) : 0;
        print_range(start, sz, os);
    }

    // Prints all rows
    void print(std::ostream& os = std::cout) const {
        print_range(0, size(), os);
    }

private:
    std::vector<point> data_;

    static std::string format_date(time_point tp) {
        using namespace std::chrono;
        const auto d = floor<days>(tp);
        const year_month_day ymd{sys_days{d}};

        std::ostringstream ss;
        ss << std::setw(4) << std::setfill('0') << int(ymd.year()) << "-"
           << std::setw(2) << std::setfill('0') << unsigned(ymd.month()) << "-"
           << std::setw(2) << std::setfill('0') << unsigned(ymd.day());
        return ss.str();
    }

    void print_range(std::size_t start, std::size_t end, std::ostream& os) const {
        // preserve stream formatting
        const auto old_flags = os.flags();
        const auto old_prec  = os.precision();

        os << "date       | value\n";
        os << "-----------|----------\n";

        os << std::fixed << std::setprecision(2);
        for (std::size_t i = start; i < end; ++i) {
            os << format_date(data_[i].first) << " | " << data_[i].second << "\n";
        }

        os.flags(old_flags);
        os.precision(old_prec);
    }
};

inline void head(const TimeSeries& s, std::size_t n = 10, std::ostream& os = std::cout) { s.head(n, os); }
inline void tail(const TimeSeries& s, std::size_t n = 10, std::ostream& os = std::cout) { s.tail(n, os); }
inline void print(const TimeSeries& s, std::ostream& os = std::cout) { s.print(os); }

}
