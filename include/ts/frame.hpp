#pragma once

#include <algorithm>
#include <chrono>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <limits>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>
#include <memory>


namespace ts {

class Frame {
public:
    using time_point = std::chrono::system_clock::time_point;

    enum class freq { day, week, month, year };
    enum class agg  { mean, sum, min, max, first, last, count };

    // -------- basic structure --------
    std::vector<time_point> index;

    std::size_t nrows() const noexcept { return index.size(); }

    void add_col(std::string name, std::vector<double> values) {
        if (!index.empty() && values.size() != index.size())
            throw std::runtime_error("add_col: column length must match index length");
        if (cols_.count(name) == 0) order_.push_back(name);
        cols_[std::move(name)] = std::move(values);
    }

    const std::vector<double>& col(const std::string& name) const {
        auto it = cols_.find(name);
        if (it == cols_.end()) throw std::runtime_error("Missing column: " + name);
        return it->second;
    }

    std::vector<double>& col_mut(const std::string& name) {
        auto it = cols_.find(name);
        if (it == cols_.end()) throw std::runtime_error("Missing column: " + name);
        return it->second;
    }

    bool has_col(const std::string& name) const { return cols_.count(name) != 0; }
    const std::vector<std::string>& col_names() const noexcept { return order_; }

    void require_aligned() const {
        for (const auto& name : order_) {
            const auto& v = cols_.at(name);
            if (v.size() != index.size())
                throw std::runtime_error("Column '" + name + "' not aligned with index.");
        }
    }

    // Row view for lambdas in filter/mutate
    class Row {
    public:
        Row(const Frame* f, std::size_t i) : f_(f), i_(i) {}
        double operator()(const std::string& c) const { return f_->col(c).at(i_); }
        const time_point& t() const { return f_->index.at(i_); }
        std::size_t i() const noexcept { return i_; }
    private:
        const Frame* f_;
        std::size_t i_;
    };

    // -------- printing --------
    const Frame& head(std::size_t n = 10, std::ostream& os = std::cout) const {
        print_range(0, std::min(n, nrows()), os);
        return *this;
    }
    const Frame& tail(std::size_t n = 10, std::ostream& os = std::cout) const {
        const auto sz = nrows();
        const auto start = (sz > n) ? (sz - n) : 0;
        print_range(start, sz, os);
        return *this;
    }
    const Frame& print(std::ostream& os = std::cout) const {
        print_range(0, nrows(), os);
        return *this;
    }

    // -------- dplyr-like verbs (chainable) --------

    // select(cols...)
    Frame select(const std::vector<std::string>& keep) const {
        require_aligned();
        Frame out;
        out.index = index;
        for (const auto& name : keep) {
            out.add_col(name, col(name));
        }
        return out;
    }

    // rename("old","new")
    Frame rename(const std::string& old_name, const std::string& new_name) const {
        require_aligned();
        if (!has_col(old_name)) throw std::runtime_error("rename: missing column: " + old_name);
        if (has_col(new_name)) throw std::runtime_error("rename: new name already exists: " + new_name);

        Frame out = *this;
        out.cols_[new_name] = out.cols_.at(old_name);
        out.cols_.erase(old_name);

        for (auto& n : out.order_) {
            if (n == old_name) { n = new_name; break; }
        }
        return out;
    }

    // filter(row_pred) where pred(Row) -> bool
    template <class Pred>
    Frame filter(Pred pred) const {
        require_aligned();
        Frame out;
        // pre-create columns
        out.index.reserve(nrows());
        for (const auto& name : order_) out.add_col(name, {});

        // reserve columns
        for (const auto& name : order_) out.col_mut(name).reserve(nrows());

        for (std::size_t i = 0; i < nrows(); ++i) {
            Row r{this, i};
            if (!pred(r)) continue;
            out.index.push_back(index[i]);
            for (const auto& name : order_) out.col_mut(name).push_back(col(name)[i]);
        }
        return out;
    }

    // mutate("new_col", fn(Row)->double)
    template <class Fn>
    Frame mutate(std::string new_col, Fn fn) const {
        require_aligned();
        Frame out = *this;
        std::vector<double> v;
        v.reserve(nrows());
        for (std::size_t i = 0; i < nrows(); ++i) v.push_back(fn(Row{this, i}));
        out.add_or_replace_(std::move(new_col), std::move(v));
        return out;
    }

    // across: apply fn(double)->double to each named col, writing col+suffix
    template <class Fn>
    Frame mutate_across(const std::vector<std::string>& cols,
                        std::string suffix,
                        Fn fn) const {
        require_aligned();
        Frame out = *this;
        for (const auto& c : cols) {
            const auto& x = col(c);
            std::vector<double> y;
            y.reserve(nrows());
            for (double xi : x) y.push_back(fn(xi));
            out.add_or_replace_(c + suffix, std::move(y));
        }
        return out;
    }

    // arrange by time or by a numeric column
    Frame arrange_time(bool desc = false) const {
        require_aligned();
        auto idx = row_index_();
        auto cmp = [&](std::size_t a, std::size_t b) {
            return desc ? (index[a] > index[b]) : (index[a] < index[b]);
        };
        std::stable_sort(idx.begin(), idx.end(), cmp);
        return permute_(idx);
    }

    Frame arrange(const std::string& by_col, bool desc = false) const {
        require_aligned();
        const auto& x = col(by_col);
        auto idx = row_index_();
        auto cmp = [&](std::size_t a, std::size_t b) {
            return desc ? (x[a] > x[b]) : (x[a] < x[b]);
        };
        std::stable_sort(idx.begin(), idx.end(), cmp);
        return permute_(idx);
    }

    // lag/lead for a given column
    Frame lag(const std::string& src, std::string dst, std::size_t k = 1,
              double fill = std::numeric_limits<double>::quiet_NaN()) const {
        require_aligned();
        const auto& x = col(src);
        std::vector<double> y(nrows(), fill);
        for (std::size_t i = k; i < nrows(); ++i) y[i] = x[i - k];
        Frame out = *this;
        out.add_or_replace_(std::move(dst), std::move(y));
        return out;
    }

    Frame lead(const std::string& src, std::string dst, std::size_t k = 1,
               double fill = std::numeric_limits<double>::quiet_NaN()) const {
        require_aligned();
        const auto& x = col(src);
        std::vector<double> y(nrows(), fill);
        for (std::size_t i = 0; i + k < nrows(); ++i) y[i] = x[i + k];
        Frame out = *this;
        out.add_or_replace_(std::move(dst), std::move(y));
        return out;
    }

    // -------- group_by / summarise --------

    class Grouped {
        public:
            Frame summarise(const std::string& value_col, agg a, std::string out_name = "") const {
                base_->require_aligned();
                const auto& v = base_->col(value_col);
                Frame out;

                out.index.reserve(groups_.size());
                for (const auto& k : keys_) out.add_col(k, {});
                const std::string sname = out_name.empty() ? (value_col + "__" + agg_name_(a)) : out_name;
                out.add_col(sname, {});

                for (const auto& g : groups_) {
                    const auto start = g.start;
                    const auto end   = g.end;
                    if (start >= end) continue;

                    out.index.push_back(base_->index[order_[end - 1]]);

                    const auto first_row = order_[start];
                    for (const auto& k : keys_) out.col_mut(k).push_back(base_->col(k)[first_row]);

                    out.col_mut(sname).push_back(agg_apply_(v, order_, start, end, a));
                }

                return out;
            }

            Frame mean (const std::string& value_col, std::string out_name = "") const { return summarise(value_col, agg::mean,  std::move(out_name)); }
            Frame sum  (const std::string& value_col, std::string out_name = "") const { return summarise(value_col, agg::sum,   std::move(out_name)); }
            Frame min  (const std::string& value_col, std::string out_name = "") const { return summarise(value_col, agg::min,   std::move(out_name)); }
            Frame max  (const std::string& value_col, std::string out_name = "") const { return summarise(value_col, agg::max,   std::move(out_name)); }
            Frame first(const std::string& value_col, std::string out_name = "") const { return summarise(value_col, agg::first, std::move(out_name)); }
            Frame last (const std::string& value_col, std::string out_name = "") const { return summarise(value_col, agg::last,  std::move(out_name)); }

            Frame count(std::string out_name = "n") const {
                Frame out;
                out.index.reserve(groups_.size());
                for (const auto& k : keys_) out.add_col(k, {});
                out.add_col(out_name, {});

                for (const auto& g : groups_) {
                    out.index.push_back(base_->index[order_[g.end - 1]]);
                    const auto first_row = order_[g.start];
                    for (const auto& k : keys_) out.col_mut(k).push_back(base_->col(k)[first_row]);
                    out.col_mut(out_name).push_back(static_cast<double>(g.end - g.start));
                }
                return out;
            }

        private:
            friend class Frame;

            struct group { std::size_t start{}, end{}; }; // [start,end)

            std::shared_ptr<Frame> keep_alive_; // if non-null, owns the base frame
            const Frame* base_;                 // always valid as long as keep_alive_ or caller frame lives
            std::vector<std::string> keys_;
            std::vector<std::size_t> order_;
            std::vector<group> groups_;

            static std::string agg_name_(agg a) {
                switch (a) {
                    case agg::mean:  return "mean";
                    case agg::sum:   return "sum";
                    case agg::min:   return "min";
                    case agg::max:   return "max";
                    case agg::first: return "first";
                    case agg::last:  return "last";
                    case agg::count: return "count";
                }
                return "agg";
            }

            static double agg_apply_(const std::vector<double>& v,
                                    const std::vector<std::size_t>& order,
                                    std::size_t start, std::size_t end,
                                    agg a) {
                auto at = [&](std::size_t j) -> double { return v[order[j]]; };

                if (a == agg::count) return static_cast<double>(end - start);

                const double first = at(start);
                const double last  = at(end - 1);

                if (a == agg::first) return first;
                if (a == agg::last)  return last;

                if (a == agg::min) {
                    double m = first;
                    for (std::size_t j = start + 1; j < end; ++j) m = std::min(m, at(j));
                    return m;
                }
                if (a == agg::max) {
                    double m = first;
                    for (std::size_t j = start + 1; j < end; ++j) m = std::max(m, at(j));
                    return m;
                }

                double s = 0.0;
                for (std::size_t j = start; j < end; ++j) s += at(j);
                return (a == agg::mean) ? (s / static_cast<double>(end - start)) : s;
            }

            void build_() {
                order_ = base_->row_index_();

                std::vector<const std::vector<double>*> keycols;
                keycols.reserve(keys_.size());
                for (const auto& k : keys_) keycols.push_back(&base_->col(k));

                auto cmp = [&](std::size_t a, std::size_t b) {
                    for (std::size_t j = 0; j < keycols.size(); ++j) {
                        const double va = (*keycols[j])[a];
                        const double vb = (*keycols[j])[b];
                        if (va < vb) return true;
                        if (va > vb) return false;
                    }
                    return base_->index[a] < base_->index[b];
                };

                std::stable_sort(order_.begin(), order_.end(), cmp);

                auto same_keys = [&](std::size_t a, std::size_t b) {
                    for (std::size_t j = 0; j < keycols.size(); ++j)
                        if ((*keycols[j])[a] != (*keycols[j])[b]) return false;
                    return true;
                };

                if (order_.empty()) return;
                std::size_t start = 0;
                for (std::size_t i = 1; i < order_.size(); ++i) {
                    if (!same_keys(order_[i - 1], order_[i])) {
                        groups_.push_back(group{start, i});
                        start = i;
                    }
                }
                groups_.push_back(group{start, order_.size()});
            }

            Grouped(const Frame& base, std::vector<std::string> keys)
                : keep_alive_(nullptr), base_(&base), keys_(std::move(keys)) {
                build_();
            }

            Grouped(std::shared_ptr<Frame> owned, std::vector<std::string> keys)
                : keep_alive_(std::move(owned)), base_(keep_alive_.get()), keys_(std::move(keys)) {
                build_();
            }
        };


    Grouped group_by(const std::vector<std::string>& keys) const {
        require_aligned();
        for (const auto& k : keys) (void)col(k); // validate exists
        return Grouped(*this, keys);
    }

    // group_by_time: adds derived key columns (__ts_year, __ts_month, __ts_week, __ts_day) then groups
    Grouped group_by_time(freq f) const {
        require_aligned();

        auto tmp = std::make_shared<Frame>(*this);

        std::vector<double> year(nrows()), month(nrows()), day(nrows()), week(nrows());
        for (std::size_t i = 0; i < nrows(); ++i) {
            using namespace std::chrono;
            const auto d = floor<days>(tmp->index[i]);
            const year_month_day ymd{sys_days{d}};
            const int y = int(ymd.year());
            const unsigned m = unsigned(ymd.month());
            const unsigned dd = unsigned(ymd.day());

            year[i]  = static_cast<double>(y);
            month[i] = static_cast<double>(m);
            day[i]   = static_cast<double>(dd);
            week[i]  = static_cast<double>(d.time_since_epoch().count() / 7);
        }

        std::vector<std::string> keys;
        if (f == freq::year) {
            tmp->add_or_replace_("__ts_year", std::move(year));
            keys = {"__ts_year"};
        } else if (f == freq::month) {
            tmp->add_or_replace_("__ts_year", std::move(year));
            tmp->add_or_replace_("__ts_month", std::move(month));
            keys = {"__ts_year", "__ts_month"};
        } else if (f == freq::week) {
            tmp->add_or_replace_("__ts_year", std::move(year));
            tmp->add_or_replace_("__ts_week", std::move(week));
            keys = {"__ts_year", "__ts_week"};
        } else {
            tmp->add_or_replace_("__ts_year", std::move(year));
            tmp->add_or_replace_("__ts_month", std::move(month));
            tmp->add_or_replace_("__ts_day", std::move(day));
            keys = {"__ts_year", "__ts_month", "__ts_day"};
        }

        return Grouped(tmp, keys);
    }

    // -------- join (time index left_join) --------
    Frame left_join_time(const Frame& rhs, std::string suffix = "_rhs") const {
        require_aligned();
        rhs.require_aligned();

        // map rhs time -> row index (requires unique rhs timestamps)
        std::unordered_map<std::int64_t, std::size_t> rhs_map;
        rhs_map.reserve(rhs.nrows());

        for (std::size_t i = 0; i < rhs.nrows(); ++i) {
            const auto k = time_key_(rhs.index[i]);
            if (!rhs_map.emplace(k, i).second)
                throw std::runtime_error("left_join_time: rhs has duplicate timestamps");
        }

        Frame out = *this;

        // add rhs columns (excluding any name collisions; suffix them)
        for (const auto& name : rhs.order_) {
            std::string out_name = name;
            if (out.has_col(out_name)) out_name += suffix;

            std::vector<double> v(nrows(), std::numeric_limits<double>::quiet_NaN());
            const auto& rv = rhs.col(name);

            for (std::size_t i = 0; i < nrows(); ++i) {
                const auto k = time_key_(index[i]);
                auto it = rhs_map.find(k);
                if (it != rhs_map.end()) v[i] = rv[it->second];
            }
            out.add_or_replace_(std::move(out_name), std::move(v));
        }

        return out;
    }

private:
    std::unordered_map<std::string, std::vector<double>> cols_;
    std::vector<std::string> order_;

    void add_or_replace_(std::string name, std::vector<double> values) {
        if (!index.empty() && values.size() != index.size())
            throw std::runtime_error("add_or_replace: column length must match index length");
        if (cols_.count(name) == 0) order_.push_back(name);
        cols_[std::move(name)] = std::move(values);
    }

    static std::string format_date_(time_point tp) {
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
        require_aligned();

        os << "date       ";
        for (const auto& name : order_) os << "| " << name << " ";
        os << "\n";

        os << "-----------";
        for (std::size_t i = 0; i < order_.size(); ++i) os << "+---------";
        os << "\n";

        const auto old_flags = os.flags();
        const auto old_prec  = os.precision();

        os << std::fixed << std::setprecision(4);
        for (std::size_t i = start; i < end; ++i) {
            os << format_date_(index[i]) << " ";
            for (const auto& name : order_) {
                os << "| " << col(name)[i] << " ";
            }
            os << "\n";
        }

        os.flags(old_flags);
        os.precision(old_prec);
    }

    std::vector<std::size_t> row_index_() const {
        std::vector<std::size_t> idx(nrows());
        std::iota(idx.begin(), idx.end(), 0);
        return idx;
    }

    Frame permute_(const std::vector<std::size_t>& idx) const {
        Frame out;
        out.index.reserve(nrows());
        for (auto i : idx) out.index.push_back(index[i]);

        for (const auto& name : order_) {
            const auto& x = col(name);
            std::vector<double> y;
            y.reserve(nrows());
            for (auto i : idx) y.push_back(x[i]);
            out.add_col(name, std::move(y));
        }
        return out;
    }

    static std::int64_t time_key_(time_point tp) {
        // seconds since epoch as key (good enough for daily/hourly data; can change later)
        const auto sec = std::chrono::time_point_cast<std::chrono::seconds>(tp);
        return sec.time_since_epoch().count();
    }
};

}
