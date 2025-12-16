#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <type_traits>
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

    // ---------- verbs (return new series) ----------

    TimeSeries lag(std::size_t k = 1,
                   double fill = std::numeric_limits<double>::quiet_NaN()) const {
        TimeSeries out;
        out.data_.reserve(size());
        for (std::size_t i = 0; i < size(); ++i) {
            const double v = (i >= k) ? data_[i - k].second : fill;
            out.data_.emplace_back(data_[i].first, v);
        }
        return out;
    }

    TimeSeries lead(std::size_t k = 1,
                    double fill = std::numeric_limits<double>::quiet_NaN()) const {
        TimeSeries out;
        out.data_.reserve(size());
        for (std::size_t i = 0; i < size(); ++i) {
            const double v = (i + k < size()) ? data_[i + k].second : fill;
            out.data_.emplace_back(data_[i].first, v);
        }
        return out;
    }

    // mutate: transform values (supports fn(double) OR fn(time_point,double) OR fn(point))
    template <class Fn>
    TimeSeries mutate(Fn fn) const {
        TimeSeries out;
        out.data_.reserve(size());

        for (std::size_t i = 0; i < size(); ++i) {
            const auto& p = data_[i];
            double nv = 0.0;

            if constexpr (std::is_invocable_r_v<double, Fn, double>) {
                nv = fn(p.second);
            } else if constexpr (std::is_invocable_r_v<double, Fn, time_point, double>) {
                nv = fn(p.first, p.second);
            } else if constexpr (std::is_invocable_r_v<double, Fn, const point&>) {
                nv = fn(p);
            } else {
                static_assert(sizeof(Fn) == 0,
                              "mutate(fn): fn must return double and accept (double) or (time_point,double) or (point)");
            }

            out.data_.emplace_back(p.first, nv);
        }
        return out;
    }

    // ---------- group_by (time series style) ----------

    class MonthlyGroup {
    public:
        explicit MonthlyGroup(const TimeSeries& s) : s_(s) {}

        TimeSeries mean() const { return summarize(Mode::Mean); }
        TimeSeries sum() const { return summarize(Mode::Sum); }
        TimeSeries last() const { return summarize(Mode::Last); }

    private:
        enum class Mode { Mean, Sum, Last };
        const TimeSeries& s_;

        static std::chrono::year_month month_key(time_point tp) {
            using namespace std::chrono;
            const auto d = floor<days>(tp);
            const year_month_day ymd{sys_days{d}};
            return ymd.year() / ymd.month();
        }

        TimeSeries summarize(Mode mode) const {
            TimeSeries out;
            if (s_.empty()) return out;

            auto cur = month_key(s_.data_[0].first);
            double acc = 0.0;
            std::size_t n = 0;
            double lastv = s_.data_[0].second;
            time_point lastt = s_.data_[0].first;

            auto flush = [&] {
                double v = 0.0;
                if (mode == Mode::Mean) v = (n ? acc / static_cast<double>(n) : std::numeric_limits<double>::quiet_NaN());
                if (mode == Mode::Sum)  v = acc;
                if (mode == Mode::Last) v = lastv;
                out.push_back(lastt, v); // stamp group with last timestamp in that month
            };

            for (const auto& [t, v] : s_.data_) {
                const auto k = month_key(t);
                if (k != cur) {
                    flush();
                    cur = k;
                    acc = 0.0;
                    n = 0;
                }
                acc += v;
                ++n;
                lastv = v;
                lastt = t;
            }
            flush();
            return out;
        }
    };

    MonthlyGroup group_by_month() const { return MonthlyGroup{*this}; }

    // ---------- printing (chainable) ----------

    const TimeSeries& head(std::size_t n = 10, std::ostream& os = std::cout) const {
        const std::size_t end = std::min(n, size());
        print_range(0, end, os);
        return *this;
    }

    const TimeSeries& tail(std::size_t n = 10, std::ostream& os = std::cout) const {
        const std::size_t sz = size();
        const std::size_t start = (sz > n) ? (sz - n) : 0;
        print_range(start, sz, os);
        return *this;
    }

    const TimeSeries& print(std::ostream& os = std::cout) const {
        print_range(0, size(), os);
        return *this;
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

// optional free functions (R-like)
inline void head(const TimeSeries& s, std::size_t n = 10, std::ostream& os = std::cout) { s.head(n, os); }
inline void tail(const TimeSeries& s, std::size_t n = 10, std::ostream& os = std::cout) { s.tail(n, os); }
inline void print(const TimeSeries& s, std::ostream& os = std::cout) { s.print(os); }

} 
