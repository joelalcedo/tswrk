
# tswrk

`tswrk` is an experimental C++20 time‑series wrangling toolkit.  It takes
inspiration from the R tidyverse: **dplyr**, which organises data
manipulation around a small set of consistent verbs, and **magrittr**,
which popularised the pipe operator for left‑to‑right composition of
operations.  In the dplyr philosophy, each common data manipulation task
is expressed as a verb (e.g. `filter()`, `arrange()`, `mutate()`)https://dplyr.tidyverse.org/articles/dplyr.html#:~:text=rows,as_tibble.
tswrk aims to provide similar verbs and a pipeline syntax in modern C++.

> **Status:** early, evolving prototype.  APIs may change.

## Building

tswrk requires a C++20 compiler (tested on AppleClang 17.0) and CMake ≥3.22.

```
cmake -S . -B build
cmake --build build
./build/tswrk
```

This uses the `CMakeLists.txt` included in the repository to build the
executable from `src/main.cpp` and bring in headers from `include/`.

## Quickstart examples

### Multi‑column frame

tswrk’s core type is a **Frame**, representing a time index and a set
of numeric columns.  You add columns with `add_col(name, values)` and
then chain verbs to transform the frame.  A typical pipeline might look
like this:

```cpp
#include "ts/frame.hpp"
#include <chrono>
#include <cmath>

int main() {
  ts::Frame df;

  // Build a daily index of 45 days starting from now
  using clock = std::chrono::system_clock;
  auto t0 = clock::now();
  for (int i = 0; i < 45; ++i)
    df.index.push_back(t0 + std::chrono::hours(24 * i));

  // Add some columns
  std::vector<double> px, vol;
  for (int i = 0; i < 45; ++i) {
    px.push_back(100.0 + std::sin(i / 3.0));
    vol.push_back((i % 7) + 1.0);
  }
  df.add_col("px", px);
  df.add_col("vol", vol);

  // Chain dplyr‑style verbs
  auto out = df
    .arrange_time()                           // sort by time
    .filter([](const ts::Frame::Row& r){ return r("vol") >= 3.0; })
    .mutate("px2", [](const ts::Frame::Row& r){ return r("px") * r("px"); })
    .lag("px", "px_l1", 1)                 // lag the px column
    .group_by_time(ts::Frame::freq::month)    // group by month
    .mean("px", "px_monthly_mean");        // summarise

  out.tail(10); // print the last 10 rows
}
```

### Single‑series time series

For univariate series there is a lightweight `TimeSeries` class.  It
stores `(time, value)` pairs and supports similar verbs:

```cpp
#include "ts/time_series.hpp"
#include <chrono>
#include <cmath>

int main() {
  ts::TimeSeries s;
  using clock = std::chrono::system_clock;
  auto t0 = clock::now();
  for (int i = 0; i < 45; ++i)
    s.push_back(t0 + std::chrono::hours(24 * i), 100.0 + std::sin(i / 3.0));

  auto monthly = s
    .arrange_time()
    .filter([](double v){ return v > 0; })
    .mutate([](double v){ return v * 1.01; })
    .lag(1)
    .group_by(ts::TimeSeries::freq::month)
    .mean();

  monthly.tail(10);
}
```

## Current API (what exists)

Below is a snapshot of implemented verbs.  Items checked have been
implemented; unchecked boxes remain TODO.

### Frame (multi‑column table)

| Functionality | Status |
|---|---|
| `add_col(name, values)`, `col(name)`, `col_mut(name)` | ✅ |
| `select({cols…})` | ✅ |
| `rename(old, new)` | ✅ |
| `filter(pred(Row)→bool)` | ✅ |
| `mutate(new_col, fn(Row)→double)` | ✅ |
| `mutate_across(cols, suffix, fn)` | ✅ |
| `arrange_time(desc=false)` | ✅ |
| `arrange(by_col, desc=false)` | ✅ |
| `lag(src, dst, k, fill)` and `lead(src, dst, k, fill)` | ✅ |
| `group_by({key_cols…})` and `group_by_time(freq)` | ✅ |
| `summarise(value_col, agg, out_name)` with helpers `.mean()`, `.sum()`, `.min()`, `.max()`, `.first()`, `.last()`, `.count()` | ✅ |
| `left_join_time(rhs, suffix="_rhs")` | ✅ |
| `head(n)`, `tail(n)`, `print()` | ✅ |

### TimeSeries (univariate)

| Functionality | Status |
|---|---|
| store `(time, value)` pairs | ✅ |
| `push_back(time, value)` | ✅ |
| `arrange_time(desc=false)` | ✅ |
| `filter(pred)` | ✅ |
| `mutate(fn)` | ✅ |
| `lag(k, fill)` and `lead(k, fill)` | ✅ |
| `group_by(freq)` and summarise (`mean`, `sum`, `min`, `max`, `first`, `last`, `count`) | ✅ |
| `head(n)`, `tail(n)`, `print()` | ✅ |

## Roadmap for dplyr/magrittr parity

dplyr organises verbs into operations on rows, columns and groupshttps://dplyr.tidyverse.org/articles/dplyr.html#:~:text=Single%20table%20verbs.  The items below
represent the gap between tswrk and full dplyr/magrittr.  We check off
what’s done and outline future modules.

### Core verbs

- [x] `filter()` – keep rows satisfying a predicate.
- [x] `arrange()`/`arrange_time()` – reorder rows by column or time.
- [x] `mutate()` – add or transform columns.
- [x] `summarise()` – collapse groups to summary statistics (mean, sum, etc.).
- [x] `select()`/`rename()` – include/exclude or rename columns.
- [ ] `distinct()` – deduplicate rows on key columns.
- [ ] `slice()`, `slice_head()`, `slice_tail()`, `slice_sample()` – row slicing.
- [ ] `transmute()` – keep only newly created columns.
- [ ] `relocate()` – reorder columns.
- [ ] `across()` with multiple functions and column specifications.

### Joins and two‑table verbs

- [x] `left_join_time(rhs)` – join on the time index.
- [ ] `inner_join_time`, `right_join_time`, `full_join_time`.
- [ ] joins on arbitrary key columns (not just time).
- [ ] inequality/rolling joins (asof joins) for nearest match.
- [ ] overlap joins (interval/range matches).

### Grouping extensions

- [x] grouping by explicit key columns (`group_by()`).
- [x] grouping by time buckets (`group_by_time()`).
- [ ] grouped `mutate()`/`transmute()` semantics (compute within group).
- [ ] window functions (rolling mean, rank, row_number, etc.).
- [ ] `ungroup()`.

### Pipeline syntax

magrittr introduced `%>%` to avoid nested function calls and make
pipelines read left‑to‑righthttps://dplyr.tidyverse.org/articles/dplyr.html#:~:text=All%20of%20the%20dplyr%20functions,the%20pipe%20operator%20as%20%E2%80%9Cthen%E2%80%9D.  tswrk currently
uses method chaining.  Planned:

- [ ] `operator|` or a `pipe()` adaptor for forward composition.
- [ ] tee/side effect pipes for logging intermediate results.

### Type and NA handling

- [ ] support for integer, boolean and string columns (currently only double).
- [ ] missing values semantics (NaN vs optional values).
- [ ] schema / variant‑based column storage.

### Data I/O and reshaping modules

- [ ] `tsio`: reading/writing CSV, Parquet or Arrow for frames and series.
- [ ] `tstidy`: pivoting/reshaping (analogous to `pivot_longer`/`pivot_wider`).
- [ ] `tslubridate`: date parsing, time zones, calendar arithmetic.
- [ ] `tsstringr`: string helpers for column names and parsing.

### Plotting: towards a ggplot2 analogue

In the future, tswrk may grow a plotting subsystem inspired by
**ggplot2**.  This could follow a grammar of graphics:

- [ ] `tsplot`: an API that takes a data frame, aesthetic mappings
  (`aes(x=…, y=…)`) and layers (`geom_line`, `geom_point`, `geom_bar`,
  `geom_histogram`) to build charts.
- [ ] scales: continuous, discrete and date/time scales.
- [ ] facets for small multiples.
- [ ] theming: a `tsthemes` module with built‑in themes (minimal,
  dark, economist, WSJ, fivethirtyeight, etc.).

### Quality and developer experience

- [ ] unit tests (e.g. Catch2/GoogleTest) and continuous integration.
- [ ] benchmarks for filter/mutate/aggregate performance.
- [ ] code formatting (`clang-format`) and static analysis (`clang-tidy`).
- [ ] example programs under `examples/` illustrating common pipelines.
- [ ] documentation generation (e.g. Doxygen) and a rendered docs site.

## Known issues

- Grouped pipelines currently store references to temporaries; this can
  cause use‑after‑free errors if intermediate objects are destroyed.  A
  `shared_ptr` strategy is being developed.
- Visual Studio Code may show “identifier ‘days’ is undefined” for
  chrono types; ensure the IntelliSense C++ language standard is set
  to C++20.  The CMake build will compile successfully even when
  IntelliSense warns.
