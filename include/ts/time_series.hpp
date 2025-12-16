#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace ts {

class TimeSeries {
public:
    using time_point = std::chrono::system_clock::time_point;
    using point      = std::pair<time_point, double>;

    enum class freq { day, week, month, year };
    enum class keep { first, last };

    // ---- basic storage ----
    void push_back(time_point t, double value) { data_.emplace_back(t, value); }
    std::size_t size() const noexcept { return data_.size(); }
    bool empty() const noexcept { return data_.empty(); }

    const point& at(std::size_t i) const { return data_.at(i); }
    time_point time_at(std::size_t i) const { return data_.at(i).first; }
    double value_at(std::size_t i) const { return data_.at(i).second; }

    // ---- printing ----
    const TimeSeries& head(std::size_t n = 10, std::ostream& os = std::cout) const {
        print_range(0, std::min(n, size()), os);
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

    // ---- dplyr-ish verbs (chainable) ----

    // filter(predicate): predicate can be (double) or (time_point,double) or (point)
    template <class Pred>
    TimeSeries filter(Pred pred) const {
        TimeSeries out;
        out.data_.reserve(size());

        for (const auto& p : data_) {
            bool keep_row = false;

            if constexpr (std::is_invocable_r_v<bool, Pred, double>) {
                keep_row = pred(p.second);
            } else if constexpr (std::is_invocable_r_v<bool, Pred, time_point, double>) {
                keep_row = pred(p.first, p.second);
            } else if constexpr (std::is_invocable_r_v<bool, Pred, const point&>) {
                keep_row = pred(p);
            } else {
                static_assert(sizeof(Pred) == 0,
                              "filter(pred): pred must return bool and accept (double) or (time_point,double) or (point)");
            }

            if (keep_row) out.data_.push_back(p);
        }

        return out;
    }

    // arrange by time or value
    TimeSeries arrange_time(bool desc = false) const {
        TimeSeries out = *this;
        auto cmp = [&](const point& a, const point& b) {
            return desc ? (a.first > b.first) : (a.first < b.first);
        };
        std::stable_sort(out.data_.begin(), out.data_.end(), cmp);
        return out;
    }

    TimeSeries arrange_value(bool desc = false) const {
        TimeSeries out = *this;
        auto cmp = [&](const point& a, const point& b) {
            return desc ? (a.second > b.second) : (a.second < b.second);
        };
        std::stable_sort(out.data_.begin(), out.data_.end(), cmp);
        return out;
    }

    // slice by [start, end) indices
    TimeSeries slice(std::size_t start, std::size_t end) const {
        if (start > end) throw std::runtime_error("slice: start > end");
        start = std::min(start, size());
        end   = std::min(end, size());

        TimeSeries out;
        out.data_.reserve(end - start);
        for (std::size_t i = start; i < end; ++i) out.data_.push_back(data_[i]);
        return out;
    }

    TimeSeries slice_head(std::size_t n = 5) const { return slice(0, std::min(n, size())); }
    TimeSeries slice_tail(std::size_t n = 5) const {
        const std::size_t sz = size();
        const std::size_t start = (sz > n) ? (sz - n) : 0;
        return slice(start, sz);
    }

    // distinct by timestamp (keep first/last)
    TimeSeries distinct_time(keep which = keep::first) const {
        if (empty()) return {};

        // easiest correct behavior: sort by time first
        TimeSeries sorted = arrange_time(false);

        TimeSeries out;
        out.data_.reserve(sorted.size());

        if (which == keep::first) {
            out.data_.push_back(sorted.data_.front());
            for (std::size_t i = 1; i < sorted.size(); ++i) {
                if (sorted.data_[i].first != sorted.data_[i - 1].first) out.data_.push_back(sorted.data_[i]);
            }
        } else { // keep::last
            // walk, but overwrite last occurrence
            out.data_.push_back(sorted.data_.front());
            for (std::size_t i = 1; i < sorted.size(); ++i) {
                if (sorted.data_[i].first != sorted.data_[i - 1].first) {
                    out.data_.push_back(sorted.data_[i]);
                } else {
                    out.data_.back() = sorted.data_[i];
                }
            }
        }

        return out;
    }

    // mutate(fn): fn returns new double; can accept (double) or (time_point,double) or (point)
    template <class Fn>
    TimeSeries mutate(Fn fn) const {
        TimeSeries out;
        out.data_.reserve(size());

        for (const auto& p : data_) {
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

    // transmute for a single-column series is effectively the same as mutate()
    template <class Fn>
    TimeSeries transmute(Fn fn) const { return mutate(std::move(fn)); }

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

    // ---- group_by + summarise ----
    class Grouped {
    public:
        enum class agg { mean, sum, first, last, min, max, count };

        TimeSeries summarise(agg a) const {
            TimeSeries out;
            out.data_.reserve(groups_.size());

            for (const auto& g : groups_) {
                const std::size_t start = g.start;
                const std::size_t end   = g.end;

                if (start >= end) continue;

                double v = std::numeric_limits<double>::quiet_NaN();

                if (a == agg::count) {
                    v = static_cast<double>(end - start);
                } else if (a == agg::first) {
                    v = s_.data_[start].second;
                } else if (a == agg::last) {
                    v = s_.data_[end - 1].second;
                } else if (a == agg::sum || a == agg::mean) {
                    double sum = 0.0;
                    for (std::size_t i = start; i < end; ++i) sum += s_.data_[i].second;
                    v = (a == agg::mean) ? (sum / static_cast<double>(end - start)) : sum;
                } else if (a == agg::min) {
                    v = s_.data_[start].second;
                    for (std::size_t i = start + 1; i < end; ++i) v = std::min(v, s_.data_[i].second);
                } else if (a == agg::max) {
                    v = s_.data_[start].second;
                    for (std::size_t i = start + 1; i < end; ++i) v = std::max(v, s_.data_[i].second);
                }

                out.data_.emplace_back(g.stamp, v); // stamp = last timestamp in group
            }

            return out;
        }

        TimeSeries mean()  const { return summarise(agg::mean); }
        TimeSeries sum()   const { return summarise(agg::sum); }
        TimeSeries first() const { return summarise(agg::first); }
        TimeSeries last()  const { return summarise(agg::last); }
        TimeSeries min()   const { return summarise(agg::min); }
        TimeSeries max()   const { return summarise(agg::max); }
        TimeSeries count() const { return summarise(agg::count); }

    private:
        friend class TimeSeries;

        struct group_info {
            std::size_t start{};
            std::size_t end{};     // one past
            time_point  stamp{};   // last time in the group
        };

        const TimeSeries& s_;
        freq f_;
        std::vector<group_info> groups_;

        static bool is_sorted_time(const std::vector<point>& x) {
            for (std::size_t i = 1; i < x.size(); ++i)
                if (x[i].first < x[i - 1].first) return false;
            return true;
        }

        static std::int64_t key(time_point tp, freq f) {
            using namespace std::chrono;
            const auto d = floor<days>(tp);
            const auto days_since = d.time_since_epoch().count(); // integer days since epoch

            if (f == freq::day)  return days_since;
            if (f == freq::week) return days_since / 7;

            const year_month_day ymd{sys_days{d}};
            const int y = int(ymd.year());
            const unsigned m = unsigned(ymd.month());

            if (f == freq::month) return static_cast<std::int64_t>(y) * 12 + static_cast<std::int64_t>(m) - 1;
            // year
            return static_cast<std::int64_t>(y);
        }

        Grouped(const TimeSeries& s, freq f) : s_(s), f_(f) {
            if (!is_sorted_time(s_.data_))
                throw std::runtime_error("group_by(): TimeSeries must be sorted by time (call arrange_time() first).");

            if (s_.empty()) return;

            std::size_t start = 0;
            std::int64_t cur = key(s_.data_[0].first, f_);

            for (std::size_t i = 1; i < s_.size(); ++i) {
                const auto k = key(s_.data_[i].first, f_);
                if (k != cur) {
                    groups_.push_back(group_info{start, i, s_.data_[i - 1].first});
                    start = i;
                    cur = k;
                }
            }
            groups_.push_back(group_info{start, s_.size(), s_.data_.back().first});
        }
    };

    Grouped group_by(freq f) const { return Grouped(*this, f); }

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

} 
