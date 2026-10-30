#include "gtest/gtest.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <numbers>
#include <random>
#include <set>
#include <vector>

#include "osr/extract/area_routing.h"
#include "osr/extract/point_in_polygon.h"

using namespace osr;

namespace {

constexpr auto kX0 = 20'000'000;  // 2.0 E
constexpr auto kY0 = 480'000'000;  // 48.0 N

struct pt {
  std::int64_t x_, y_;
};

using ring_spec = std::vector<pt>;

// Vertex ids are assigned consecutively over all rings, starting at 0, and
// continue with the entries that are not on the outline.
area make_area(std::vector<ring_spec> const& rings,
               std::set<std::uint64_t> const& entries,
               ring_spec const& interior = {}) {
  auto a = area{};
  auto id = std::uint64_t{0U};
  for (auto const& r : rings) {
    auto& ring = a.rings_.emplace_back();
    for (auto const& p : r) {
      ring.push_back({.id_ = id,
                      .pos_ = point{static_cast<std::int32_t>(kX0 + p.x_),
                                    static_cast<std::int32_t>(kY0 + p.y_)},
                      .entry_ = entries.contains(id)});
      ++id;
    }
  }
  for (auto const& p : interior) {
    a.interior_.push_back(
        {.id_ = id++,
         .pos_ = point{static_cast<std::int32_t>(kX0 + p.x_),
                       static_cast<std::int32_t>(kY0 + p.y_)},
         .entry_ = true});
  }
  return a;
}

std::vector<std::vector<std::uint64_t>> way_ids(area const& a,
                                                area_routing_result const& r) {
  auto ways = std::vector<std::vector<std::uint64_t>>{};
  for (auto const& w : r.ways_) {
    auto& ids = ways.emplace_back();
    for (auto const& v : w) {
      ids.push_back(get_vertex(a, v).id_);
    }
  }
  return ways;
}

struct metric {
  explicit metric(area const& a)
      : x_scale_{kUnit * std::cos(a.rings_[0][0].pos_.as_latlng().lat_ *
                                  std::numbers::pi / 180.0)} {}

  double operator()(point const a, point const b) const {
    return std::hypot((a.lat_ - b.lat_) * x_scale_, (a.lng_ - b.lng_) * kUnit);
  }

  static constexpr auto kUnit = 6371000.0 * std::numbers::pi / 180.0 * 1e-7;
  double x_scale_;
};

using matrix = std::vector<std::vector<double>>;

constexpr auto kInf = std::numeric_limits<double>::max() / 4.0;

void floyd_warshall(matrix& m) {
  auto const n = m.size();
  for (auto k = 0U; k != n; ++k) {
    for (auto i = 0U; i != n; ++i) {
      for (auto j = 0U; j != n; ++j) {
        m[i][j] = std::min(m[i][j], m[i][k] + m[k][j]);
      }
    }
  }
}

struct flat {
  explicit flat(area const& a) {
    for (auto const& r : a.rings_) {
      ring_begin_.push_back(v_.size());
      for (auto const& v : r) {
        v_.push_back(v);
      }
      if (ring_begin_.size() ==
          std::min(a.rings_.size(), a.n_routable_rings_)) {
        n_routable_ = v_.size();
      }
    }
    ring_begin_.push_back(v_.size());
    n_outline_ = v_.size();
    n_routable_ = std::min(n_routable_, n_outline_);
    for (auto const& v : a.interior_) {
      v_.push_back(v);
    }
    for (auto const& l : a.lines_) {
      line_begin_.push_back(v_.size());
      for (auto const& v : l) {
        v_.push_back(v);
      }
    }
    n_rings_ = a.rings_.size();
  }

  std::size_t idx(area_routing_result::vertex_ref const v) const {
    return v.ring_ == area_routing_result::vertex_ref::kInterior
               ? n_outline_ + v.idx_
               : (v.ring_ < n_rings_
                      ? ring_begin_[v.ring_] + v.idx_
                      : line_begin_[v.ring_ - n_rings_] + v.idx_);
  }

  // Only for outline vertices.
  std::size_t next(std::size_t const i) const {
    auto const r = std::distance(
        begin(ring_begin_),
        std::upper_bound(begin(ring_begin_), end(ring_begin_), i)) - 1;
    return i + 1U == ring_begin_[static_cast<std::size_t>(r) + 1U]
               ? ring_begin_[static_cast<std::size_t>(r)]
               : i + 1U;
  }

  // Outline vertices, followed by the entries that are not on the outline
  // (which have to be strictly inside the area) and the line vertices.
  // Index = vertex id (see make_area).
  std::vector<area::vertex> v_;
  std::vector<std::size_t> ring_begin_, line_begin_;
  std::size_t n_outline_, n_rings_;
  std::size_t n_routable_{std::numeric_limits<std::size_t>::max()};
};

int sgn(std::int64_t const x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }

int orient(pt const a, pt const b, pt const c) {
  return sgn((b.x_ - a.x_) * (c.y_ - a.y_) - (b.y_ - a.y_) * (c.x_ - a.x_));
}

// Reference: all shortest distances inside the area between all vertices,
// from the unfiltered visibility graph. A segment is walkable if it neither
// crosses the outline nor passes through an outline vertex, and its midpoint
// is inside the area (even-odd rule on doubled coordinates). Entries that
// are not on the outline are just more points to connect, links are not used.
matrix reference_distances(area const& a) {
  auto const f = flat{a};
  auto const m = metric{a};
  auto const n_all = f.v_.size();
  auto const n = f.n_outline_;
  auto const at = [&](std::size_t const i) {
    return pt{2 * static_cast<std::int64_t>(f.v_[i].pos_.lat_),
              2 * static_cast<std::int64_t>(f.v_[i].pos_.lng_)};
  };

  auto const inside = [&](pt const p) {
    auto in = false;
    for (auto i = 0U; i != n; ++i) {
      auto const a1 = at(i);
      auto const a2 = at(f.next(i));
      if ((a1.y_ > p.y_) != (a2.y_ > p.y_)) {
        auto const lo = a1.y_ < a2.y_ ? a1 : a2;
        auto const hi = a1.y_ < a2.y_ ? a2 : a1;
        if (orient(lo, hi, p) > 0) {
          in = !in;
        }
      }
    }
    return in;
  };

  auto const walkable = [&](std::size_t const i, std::size_t const j) {
    auto const p = at(i);
    auto const q = at(j);
    for (auto k = 0U; k != n; ++k) {
      auto const c = at(k);
      if (orient(p, q, c) == 0 &&
          std::min(p.x_, q.x_) <= c.x_ && c.x_ <= std::max(p.x_, q.x_) &&
          std::min(p.y_, q.y_) <= c.y_ && c.y_ <= std::max(p.y_, q.y_) &&
          !(c.x_ == p.x_ && c.y_ == p.y_) && !(c.x_ == q.x_ && c.y_ == q.y_)) {
        return false;
      }
      auto const e = at(f.next(k));
      if (orient(p, q, c) * orient(p, q, e) < 0 &&
          orient(c, e, p) * orient(c, e, q) < 0) {
        return false;
      }
    }
    return inside(pt{(p.x_ + q.x_) / 2, (p.y_ + q.y_) / 2});
  };

  auto d = matrix(n_all, std::vector<double>(n_all, kInf));
  for (auto i = 0U; i != n_all; ++i) {
    d[i][i] = 0.0;
    if (i < n) {
      auto const nx = f.next(i);
      d[i][nx] = d[nx][i] = m(f.v_[i].pos_, f.v_[nx].pos_);
    }
    for (auto j = i + 1U; j != n_all; ++j) {
      if (walkable(i, j)) {
        d[i][j] = d[j][i] = m(f.v_[i].pos_, f.v_[j].pos_);
      }
    }
  }
  floyd_warshall(d);
  return d;
}

// Distances using only the generated ways (plus the outline if routable,
// plus the links). Way ends are identified by end_ids_: vertex ids, or keys
// of one side of a line vertex (extra nodes behind the vertices).
matrix result_distances(area const& a,
                        area_routing_result const& r,
                        bool const with_outline) {
  auto const f = flat{a};
  auto const m = metric{a};
  auto keys = std::vector<std::uint64_t>{};
  for (auto const& [from, to] : r.end_ids_) {
    for (auto const k : {from, to}) {
      if (k >= f.v_.size()) {
        keys.push_back(k);
      }
    }
  }
  std::ranges::sort(keys);
  keys.erase(std::unique(begin(keys), end(keys)), end(keys));
  auto const node = [&](std::uint64_t const k) {
    return k < f.v_.size()
               ? k
               : f.v_.size() + static_cast<std::size_t>(
                                   std::ranges::lower_bound(keys, k) -
                                   begin(keys));
  };
  auto const n = f.v_.size() + keys.size();
  auto d = matrix(n, std::vector<double>(n, kInf));
  for (auto i = 0U; i != n; ++i) {
    d[i][i] = 0.0;
    if (with_outline && i < f.n_routable_) {
      auto const nx = f.next(i);
      if (!f.v_[i].blocked_ && !f.v_[nx].blocked_) {
        d[i][nx] = d[nx][i] = m(f.v_[i].pos_, f.v_[nx].pos_);
      }
    }
  }
  for (auto const& l : a.links_) {
    // Like the ways: only between entries.
    if (l.from_ < n && l.to_ < n && f.v_[l.from_].entry_ &&
        f.v_[l.to_].entry_) {
      d[l.from_][l.to_] = std::min(d[l.from_][l.to_], l.length_);
      d[l.to_][l.from_] = d[l.from_][l.to_];
    }
  }
  auto n_ends = std::vector<unsigned>(n);
  EXPECT_EQ(r.ways_.size(), r.end_ids_.size());
  for (auto k = 0U; k != r.ways_.size(); ++k) {
    auto const& w = r.ways_[k];
    auto idx = std::vector<std::size_t>{};
    for (auto const& v : w) {
      idx.push_back(f.idx(v));
    }
    // Only the ends are connected to the graph.
    auto len = 0.0;
    for (auto i = 0U; i + 1U < idx.size(); ++i) {
      len += m(f.v_[idx[i]].pos_, f.v_[idx[i + 1U]].pos_);
    }
    auto const from = node(r.end_ids_[k].first);
    auto const to = node(r.end_ids_[k].second);
    d[from][to] = std::min(d[from][to], len);
    d[to][from] = d[from][to];
    ++n_ends[from];
    ++n_ends[to];
  }
  // Ways end at entries or at junctions of more than two ways.
  for (auto i = 0U; i != n; ++i) {
    EXPECT_TRUE(n_ends[i] == 0U ||
                (i < f.v_.size() && f.v_[i].entry_ && !f.v_[i].blocked_) ||
                n_ends[i] > 2U)
        << i;
  }
  floyd_warshall(d);
  return d;
}

// All entries are connected within the allowed stretch of the optimum.
void expect_optimal(area const& a, area_routing_options const& opt) {
  auto const r = compute_area_ways(a, opt);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto const f = flat{a};
  auto const ref = reference_distances(a);
  auto const res = result_distances(a, r, opt.ring_routable_);
  for (auto i = 0U; i != f.v_.size(); ++i) {
    for (auto j = 0U; j != f.v_.size(); ++j) {
      if (f.v_[i].entry_ && f.v_[j].entry_ && !f.v_[i].blocked_ &&
          !f.v_[j].blocked_) {
        auto const tolerance = 2e-3 + ref[i][j] * 2e-5;
        EXPECT_GE(res[i][j], ref[i][j] - tolerance)
            << "from " << i << " to " << j;
        EXPECT_LE(res[i][j], ref[i][j] * opt.max_stretch_ + tolerance)
            << "from " << i << " to " << j;
      }
    }
  }
}

// Star-shaped polygon around the origin with optional square holes and
// entries at random places strictly inside.
area random_area(std::mt19937& rng, unsigned const n_outer,
                 unsigned const n_holes, double const entry_share,
                 unsigned const n_interior = 0U) {
  auto real = std::uniform_real_distribution<double>{0.0, 1.0};
  auto rings = std::vector<ring_spec>{};

  auto& outer = rings.emplace_back();
  for (auto i = 0U; i != n_outer; ++i) {
    auto const angle =
        2.0 * std::numbers::pi * (i + 0.8 * real(rng)) / n_outer;
    auto const radius = 7000.0 + 6000.0 * real(rng);
    outer.push_back({std::llround(radius * std::cos(angle)),
                     std::llround(radius * std::sin(angle))});
  }

  // Holes on a 4x4 grid well inside the minimal radius (needs >= 12 outer
  // vertices to keep the outline away from the grid).
  auto cells = std::vector<int>{};
  for (auto i = 0; i != 16; ++i) {
    cells.push_back(i);
  }
  std::shuffle(begin(cells), end(cells), rng);
  for (auto h = 0U; h != n_holes; ++h) {
    auto const cx = -3000 + 2000 * (cells[h] % 4);
    auto const cy = -3000 + 2000 * (cells[h] / 4);
    auto const half = 200 + static_cast<int>(600 * real(rng));
    auto const skew = static_cast<int>(300 * real(rng));
    rings.push_back({{cx - half, cy - half + skew},
                     {cx + half, cy - half},
                     {cx + half - skew, cy + half},
                     {cx - half, cy + half}});
  }

  auto entries = std::set<std::uint64_t>{};
  auto n_vertices = std::uint64_t{0U};
  for (auto const& r : rings) {
    n_vertices += r.size();
  }
  for (auto i = std::uint64_t{0U}; i != n_vertices; ++i) {
    if (real(rng) < entry_share) {
      entries.insert(i);
    }
  }
  entries.insert(0U);
  entries.insert(n_outer / 2U);
  auto a = make_area(rings, entries);

  // Small grid: collinear points and points on lines of sight are likely.
  for (auto i = 0U; a.interior_.size() != n_interior && i != 10'000U; ++i) {
    auto const pos = point{
        static_cast<std::int32_t>(kX0 - 13000 + 1000 * (rng() % 27U)),
        static_cast<std::int32_t>(kY0 - 13000 + 1000 * (rng() % 27U))};
    if (locate_in_rings(a.rings_, pos, &area::vertex::pos_) ==
            point_location::kInside &&
        std::ranges::none_of(a.interior_, [&](area::vertex const& v) {
          return v.pos_.lat_ == pos.lat_ && v.pos_.lng_ == pos.lng_;
        })) {
      a.interior_.push_back({.id_ = n_vertices + a.interior_.size(),
                             .pos_ = pos,
                             .entry_ = true});
    }
  }
  return a;
}

}  // namespace

TEST(area_routing, square_diagonal) {
  auto const a = make_area(
      {{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 2U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{0U, 2U}}),
            way_ids(a, r));
}

TEST(area_routing, winding_order_does_not_matter) {
  auto const a = make_area(
      {{{0, 1000}, {1000, 1000}, {1000, 0}, {0, 0}}}, {0U, 2U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{0U, 2U}}),
            way_ids(a, r));
}

TEST(area_routing, neighbours_use_outline) {
  auto const a = make_area(
      {{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 1U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_TRUE(r.ways_.empty());

  auto const with_outline =
      compute_area_ways(a, {.ring_routable_ = false});
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{0U, 1U}}),
            way_ids(a, with_outline));
}

TEST(area_routing, collinear_entries_use_outline) {
  auto const a = make_area(
      {{{0, 0}, {500, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 2U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_TRUE(r.ways_.empty());
}

TEST(area_routing, blocked_vertex_cuts_outline) {
  auto a = make_area(
      {{{0, 0}, {500, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 2U});
  a.rings_[0][1].blocked_ = true;
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{0U, 1U, 2U}}),
            way_ids(a, r));
}

TEST(area_routing, blocked_vertex_is_no_entry) {
  // Nothing passes through 1: 0 and 2 get their own way.
  auto a = make_area(
      {{{0, 0}, {500, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 1U, 2U});
  a.rings_[0][1].blocked_ = true;
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{0U, 1U, 2U}}),
            way_ids(a, r));
}

TEST(area_routing, visibility_graph) {
  // Square: 4 sides + 2 diagonals.
  auto const square = make_area(
      {{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 1U, 2U, 3U});
  auto const r = compute_area_ways(square, {.all_visibility_edges_ = true});
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ(6U, r.ways_.size());
  EXPECT_TRUE(r.junctions_.empty());

  // L-shape: the reflex corner 3 is a junction, the convex corner 4 is only
  // geometry of the outline from 3 to 5.
  auto const l = make_area({{{0, 0},
                             {1000, 0},
                             {1000, 500},
                             {500, 500},
                             {500, 1000},
                             {0, 1000}}},
                           {0U, 1U, 2U, 5U});
  auto const rl = compute_area_ways(l, {.all_visibility_edges_ = true});
  ASSERT_EQ(area_routing_status::kOk, rl.status_);
  EXPECT_EQ((std::vector<std::uint64_t>{3U}), rl.junctions_);
  expect_optimal(l, {.max_stretch_ = 1.0, .all_visibility_edges_ = true});
}

TEST(area_routing, less_than_two_entries) {
  auto const a = make_area(
      {{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {1U});
  EXPECT_EQ(area_routing_status::kNoEntries, compute_area_ways(a).status_);
}

TEST(area_routing, bend_at_reflex_corner) {
  // L shape, reflex corner is vertex 3.
  auto const a = make_area(
      {{{0, 0}, {2000, 0}, {2000, 1000}, {1000, 1000}, {1000, 2000},
        {0, 2000}}},
      {1U, 4U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{1U, 3U, 4U}}),
            way_ids(a, r));
}

TEST(area_routing, path_through_entry_is_split) {
  // Same L shape, the reflex corner is an entry as well.
  auto const a = make_area(
      {{{0, 0}, {2000, 0}, {2000, 1000}, {1000, 1000}, {1000, 2000},
        {0, 2000}}},
      {1U, 3U, 4U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto ways = way_ids(a, r);
  std::ranges::sort(ways);
  // 3-4 is part of the outline.
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{1U, 3U}}), ways);
}

TEST(area_routing, around_hole) {
  // Entries left and right of a hole that blocks the direct line.
  auto const a = make_area({{{0, 0}, {4000, 0}, {4000, 1500}, {4000, 3000},
                             {0, 3000}, {0, 1500}},
                            {{1500, 1000}, {2500, 1000}, {2500, 2100},
                             {1500, 2100}}},
                           {2U, 5U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  // Passing below the hole is shorter: via vertices 6 and 7.
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{2U, 7U, 6U, 5U}}),
            way_ids(a, r));
  expect_optimal(a, {.max_stretch_ = 1.0});
}

TEST(area_routing, stretch) {
  // Going from 0 to 3 via 2 is 1% longer than the diagonal.
  auto const a = make_area(
      {{{0, 0}, {1000, 0}, {1000, 950}, {1000, 1000}, {0, 1000}}},
      {0U, 2U, 3U});

  auto const with_stretch = compute_area_ways(a, {.max_stretch_ = 1.02});
  ASSERT_EQ(area_routing_status::kOk, with_stretch.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{0U, 2U}}),
            way_ids(a, with_stretch));

  auto const exact = compute_area_ways(a, {.max_stretch_ = 1.0});
  ASSERT_EQ(area_routing_status::kOk, exact.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{0U, 2U}, {0U, 3U}}),
            way_ids(a, exact));
}

TEST(area_routing, ways_per_entry_limit) {
  auto const a = make_area(
      {{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 2U});
  auto const r = compute_area_ways(a, {.max_ways_per_node_ = 0U});
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_TRUE(r.ways_.empty());
  EXPECT_EQ(1U, r.n_capped_);
}

TEST(area_routing, too_many_candidates) {
  auto const a = make_area(
      {{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}}, {0U, 1U, 2U});
  EXPECT_EQ(area_routing_status::kTooLarge,
            compute_area_ways(a, {.max_candidates_ = 2U}).status_);
}

TEST(area_routing, random_areas_match_reference) {
  auto rng = std::mt19937{4711U};
  for (auto i = 0U; i != 300U; ++i) {
    auto const n_outer = 4U + rng() % 40U;
    auto const n_holes = n_outer < 12U ? 0U : rng() % 6U;
    auto const a = random_area(rng, n_outer, n_holes, 0.15);
    SCOPED_TRACE(i);
    for (auto const stretch : {1.0, 1.02, 1.2}) {
      for (auto const ring_routable : {true, false}) {
        expect_optimal(a, {.max_stretch_ = stretch,
                           .max_ways_per_node_ = 1000U,
                           .ring_routable_ = ring_routable});
      }
    }
  }
}

TEST(area_routing, random_areas_with_blocked_vertices_match_reference) {
  auto rng = std::mt19937{815U};
  for (auto i = 0U; i != 300U; ++i) {
    auto const n_outer = 4U + rng() % 40U;
    auto const n_holes = n_outer < 12U ? 0U : rng() % 6U;
    auto a = random_area(rng, n_outer, n_holes, 0.3);
    for (auto& ring : a.rings_) {
      for (auto& v : ring) {
        v.blocked_ = rng() % 5U == 0U;
      }
    }
    if (compute_area_ways(a).status_ == area_routing_status::kNoEntries) {
      continue;
    }
    SCOPED_TRACE(i);
    for (auto const stretch : {1.0, 1.02, 1.2}) {
      expect_optimal(a,
                     {.max_stretch_ = stretch, .max_ways_per_node_ = 1000U});
    }
  }
}

TEST(area_routing, interior_entry_sees_in_all_directions) {
  // Entry 4 in the center: straight to all corners.
  auto const a = make_area({{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}},
                           {0U, 1U, 2U, 3U}, {{500, 500}});
  auto const r = compute_area_ways(a, {.max_stretch_ = 1.0});
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ(5U, r.n_entries_);
  EXPECT_EQ(1U, r.n_interior_entries_);
  auto ways = way_ids(a, r);
  std::ranges::sort(ways);
  // The diagonals pass through the center: no ways 0-2 and 1-3.
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{
                {0U, 4U}, {1U, 4U}, {2U, 4U}, {3U, 4U}}),
            ways);
  expect_optimal(a, {.max_stretch_ = 1.0});
}

TEST(area_routing, interior_entries_only) {
  // No entry on the outline: the two entries inside are connected around
  // the reflex corner 3 of the L shape.
  auto const a = make_area(
      {{{0, 0}, {2000, 0}, {2000, 1000}, {1000, 1000}, {1000, 2000},
        {0, 2000}}},
      {}, {{1500, 500}, {500, 1500}});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ(2U, r.n_entries_);
  EXPECT_EQ(2U, r.n_interior_entries_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{6U, 3U, 7U}}),
            way_ids(a, r));
  expect_optimal(a, {.max_stretch_ = 1.0});
}

TEST(area_routing, interior_entry_does_not_block) {
  // Entry 4 is exactly on the diagonal between the entries 0 and 2, entry 5
  // is next to it. Neither blocks the line of sight; the diagonal is split
  // at entry 4.
  auto const a = make_area({{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}},
                           {0U, 2U}, {{300, 300}, {700, 600}});
  auto const r = compute_area_ways(a, {.max_stretch_ = 1.0});
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ(6U, r.n_visibility_edges_);
  auto ways = way_ids(a, r);
  std::ranges::sort(ways);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{
                {0U, 4U}, {0U, 5U}, {2U, 4U}, {2U, 5U}, {4U, 5U}}),
            ways);
  expect_optimal(a, {.max_stretch_ = 1.0});
}

TEST(area_routing, interior_entry_behind_hole) {
  // The hole blocks the line from entry 5 (west) to the entry inside (east
  // of the hole).
  auto const a = make_area({{{0, 0}, {4000, 0}, {4000, 1500}, {4000, 3000},
                             {0, 3000}, {0, 1500}},
                            {{1500, 1000}, {2500, 1000}, {2500, 2100},
                             {1500, 2100}}},
                           {5U}, {{3000, 1500}});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{{5U, 6U, 7U, 10U}}),
            way_ids(a, r));
  expect_optimal(a, {.max_stretch_ = 1.0});
}

TEST(area_routing, interior_vertices_not_inside_are_ignored) {
  // On a vertex, on an edge, on a spike, outside, in the hole, on the
  // outline of the hole, far away, and with the id of an outline vertex.
  auto a = make_area({{{0, 0}, {2000, 0}, {2000, 1000}, {3000, 1000},
                       {2000, 1000}, {2000, 2000}, {0, 2000}},
                      {{500, 500}, {1000, 500}, {1000, 1000}, {500, 1000}}},
                     {0U, 5U},
                     {{0, 0},
                      {1000, 0},
                      {2500, 1000},
                      {2500, 500},
                      {750, 750},
                      {500, 750},
                      {1'000'000'000, 0}});
  a.interior_.push_back(
      {.id_ = 1U, .pos_ = a.rings_[0][1].pos_, .entry_ = true});
  a.interior_.push_back(
      {.id_ = 2U, .pos_ = point{kX0 + 250, kY0 + 250}, .entry_ = true});
  auto const with = compute_area_ways(a, {.max_stretch_ = 1.0});
  ASSERT_EQ(area_routing_status::kOk, with.status_);
  EXPECT_EQ(2U, with.n_entries_);
  EXPECT_EQ(0U, with.n_interior_entries_);

  auto without = a;
  without.interior_.clear();
  EXPECT_EQ(way_ids(without, compute_area_ways(without, {.max_stretch_ = 1.0})),
            way_ids(a, with));

  // Without any other entry, there is nothing to connect.
  auto const none = make_area({{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}},
                              {0U}, {{500, 0}});
  EXPECT_EQ(area_routing_status::kNoEntries, compute_area_ways(none).status_);
}

TEST(area_routing, links_replace_ways) {
  // Footway from corner 0 to the entry 4 in the center and on to entry 5.
  auto a = make_area({{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}},
                     {0U, 2U}, {{500, 500}, {700, 500}});
  auto const m = metric{a};
  auto const pos = [&](std::uint64_t const id) {
    return id < 4U ? a.rings_[0][id].pos_ : a.interior_[id - 4U].pos_;
  };
  a.links_.push_back({0U, 4U, m(pos(0U), pos(4U))});
  a.links_.push_back({4U, 5U, m(pos(4U), pos(5U))});
  // Not between entries: ignored.
  a.links_.push_back({1U, 5U, 1.0});
  a.links_.push_back({5U, 77U, 1.0});

  auto const r = compute_area_ways(a, {.max_stretch_ = 1.0});
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto ways = way_ids(a, r);
  std::ranges::sort(ways);
  // 0-4 and 4-5 exist, 0-2 runs through 4.
  EXPECT_EQ((std::vector<std::vector<std::uint64_t>>{
                {0U, 5U}, {2U, 4U}, {2U, 5U}}),
            ways);
  expect_optimal(a, {.max_stretch_ = 1.0});

  // A link that is a detour does not replace the way.
  a.links_[1].length_ *= 1.5;
  auto const detour = compute_area_ways(a, {.max_stretch_ = 1.2});
  EXPECT_TRUE(std::ranges::any_of(way_ids(a, detour), [](auto const& w) {
    return w == std::vector<std::uint64_t>{4U, 5U};
  }));
  expect_optimal(a, {.max_stretch_ = 1.2});
}

TEST(area_routing, interior_entries_count_for_limits) {
  auto const a = make_area({{{0, 0}, {1000, 0}, {1000, 1000}, {0, 1000}}},
                           {0U, 2U}, {{300, 500}, {700, 500}});
  EXPECT_EQ(area_routing_status::kTooLarge,
            compute_area_ways(a, {.max_entries_ = 3U}).status_);
  EXPECT_EQ(area_routing_status::kTooLarge,
            compute_area_ways(a, {.max_candidates_ = 3U}).status_);
  EXPECT_EQ(area_routing_status::kOk,
            compute_area_ways(a, {.max_candidates_ = 4U, .max_entries_ = 4U})
                .status_);

  auto const capped = compute_area_ways(
      a, {.max_stretch_ = 1.0, .max_ways_per_node_ = 1U});
  ASSERT_EQ(area_routing_status::kOk, capped.status_);
  auto n_ways = std::vector<unsigned>(6U);
  for (auto const& w : way_ids(a, capped)) {
    ++n_ways[w.front()];
    ++n_ways[w.back()];
  }
  EXPECT_TRUE(
      std::ranges::all_of(n_ways, [](auto const n) { return n <= 1U; }));
  EXPECT_NE(0U, capped.n_capped_);
}

TEST(area_routing, random_areas_with_interior_entries_match_reference) {
  auto rng = std::mt19937{1234U};
  auto real = std::uniform_real_distribution<double>{0.0, 1.0};
  auto n_interior = 0U;
  for (auto i = 0U; i != 300U; ++i) {
    auto const n_outer = 4U + rng() % 40U;
    auto const n_holes = n_outer < 12U ? 0U : rng() % 6U;
    // Every fifth area has no entry on the outline but the two fixed ones.
    auto a = random_area(rng, n_outer, n_holes, i % 5U == 0U ? 0.0 : 0.12,
                         1U + rng() % 8U);
    ASSERT_FALSE(a.interior_.empty());
    n_interior += a.interior_.size();
    SCOPED_TRACE(i);

    auto const check = [&]() {
      for (auto const stretch : {1.0, 1.02, 1.2}) {
        for (auto const ring_routable : {true, false}) {
          auto const opt =
              area_routing_options{.max_stretch_ = stretch,
                                   .max_ways_per_node_ = 1000U,
                                   .ring_routable_ = ring_routable};
          EXPECT_EQ(a.interior_.size(),
                    compute_area_ways(a, opt).n_interior_entries_);
          expect_optimal(a, opt);
        }
      }
    };
    check();

    // Existing connections between entries that are at least as long as
    // the shortest path: fewer ways, same guarantee.
    auto const f = flat{a};
    auto const ref = reference_distances(a);
    auto entries = std::vector<std::uint64_t>{};
    for (auto j = 0U; j != f.v_.size(); ++j) {
      if (f.v_[j].entry_) {
        entries.push_back(j);
      }
    }
    for (auto j = 0U; j != 1U + entries.size() / 2U; ++j) {
      auto const from = entries[rng() % entries.size()];
      auto const to = entries[rng() % entries.size()];
      if (from != to && ref[from][to] < kInf) {
        a.links_.push_back(
            {from, to, ref[from][to] * (1.0 + (j % 2U) * 0.3 * real(rng))});
      }
    }
    check();
  }
  EXPECT_GT(n_interior, 1000U);
}

namespace {

// Appends a line, vertex ids continue after the rings, the interior points
// and the lines added before. `entries` / `blocked`: indices on the line.
void add_line(area& a,
              ring_spec const& line,
              std::set<std::size_t> const& entries = {},
              std::set<std::size_t> const& blocked = {}) {
  auto id = a.interior_.size();
  for (auto const& r : a.rings_) {
    id += r.size();
  }
  for (auto const& l : a.lines_) {
    id += l.size();
  }
  auto& out = a.lines_.emplace_back();
  for (auto i = 0U; i != line.size(); ++i) {
    out.push_back({.id_ = id + i,
                   .pos_ = point{static_cast<std::int32_t>(kX0 + line[i].x_),
                                 static_cast<std::int32_t>(kY0 + line[i].y_)},
                   .entry_ = entries.contains(i),
                   .blocked_ = blocked.contains(i)});
  }
}

double dist(area const& a, std::vector<pt> const& path) {
  auto const m = metric{a};
  auto d = 0.0;
  for (auto i = 1U; i < path.size(); ++i) {
    d += m(point{static_cast<std::int32_t>(kX0 + path[i - 1U].x_),
                 static_cast<std::int32_t>(kY0 + path[i - 1U].y_)},
           point{static_cast<std::int32_t>(kX0 + path[i].x_),
                 static_cast<std::int32_t>(kY0 + path[i].y_)});
  }
  return d;
}

// 10000 x 10000 square, entries in the middle of the west side (5) and of
// the east side (2).
ring_spec west_east_square_ring() {
  return {{0, 0},         {10000, 0}, {10000, 5000},
          {10000, 10000}, {0, 10000}, {0, 5000}};
}

area west_east_square() {
  return make_area({west_east_square_ring()}, {2U, 5U});
}

constexpr auto kW = pt{0, 5000};
constexpr auto kE = pt{10000, 5000};

}  // namespace

TEST(area_routing, barrier_is_walked_around) {
  auto a = west_east_square();
  add_line(a, {{5000, 2000}, {5000, 8000}});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto const d = result_distances(a, r, true)[5][2];
  auto const expected = dist(a, {kW, {5000, 8000}, kE});
  EXPECT_GE(d, expected - 1e-3);
  EXPECT_LE(d, expected * 1.02 + 1e-3);
}

TEST(area_routing, barrier_sides_are_not_connected) {
  // Zig-zag: walking along the segment between the two bends would switch
  // from the west to the east side (W -> 7 -> 8 -> E: 11.1 vs 12.8 units).
  auto a = west_east_square();
  add_line(a, {{5000, 1000}, {4000, 4000}, {6000, 6000}, {5000, 9000}});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto const d = result_distances(a, r, true)[5][2];
  auto const expected = dist(a, {kW, {5000, 9000}, kE});
  EXPECT_GE(d, expected - 1e-3);
  EXPECT_LE(d, expected * 1.02 + 1e-3);
  EXPECT_LT(dist(a, {kW, {4000, 4000}, {6000, 6000}, kE}) * 1.1, expected);
}

TEST(area_routing, barrier_gate) {
  // Wall from the south to the north side with a gate (entry) in the middle.
  auto a = west_east_square();
  add_line(a, {{5000, 0}, {5000, 5000}, {5000, 10000}}, {1U});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto const d = result_distances(a, r, false);
  EXPECT_NEAR(dist(a, {kW, kE}), d[5][2], 1e-3);  // 5 -> gate (7) -> 2

  // Blocked gate: no way through, only along the outline.
  auto b = west_east_square();
  add_line(b, {{5000, 0}, {5000, 5000}, {5000, 10000}}, {1U}, {1U});
  auto const rb = compute_area_ways(b);
  ASSERT_EQ(area_routing_status::kOk, rb.status_);
  EXPECT_GE(result_distances(b, rb, false)[5][2], kInf / 2.0);
  EXPECT_NEAR(dist(b, {kW, {0, 10000}, {10000, 10000}, kE}),
              result_distances(b, rb, true)[5][2], 1e-3);
}

TEST(area_routing, closed_barrier) {
  // A fence around the center, the entry inside (6) cannot be reached.
  auto a = make_area({west_east_square_ring()}, {2U, 5U}, {{5000, 5000}});
  add_line(a, {{3000, 3000},
               {7000, 3000},
               {7000, 7000},
               {3000, 7000},
               {3000, 3000}});
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto const d = result_distances(a, r, false);
  EXPECT_NEAR(dist(a, {kW, {3000, 3000}, {7000, 3000}, kE}), d[5][2],
              dist(a, {kW, kE}) * 0.02);
  EXPECT_GE(d[5][6], kInf / 2.0);
  EXPECT_GE(d[2][6], kInf / 2.0);
}

TEST(area_routing, building_wall_is_not_routable) {
  // Building (hole) with two entrances (7, 8). Its wall is no way: the
  // entrances are connected by a generated way along it.
  auto rings = std::vector<ring_spec>{west_east_square_ring(),
                                      {{3000, 3000},
                                       {7000, 3000},
                                       {7000, 7000},
                                       {3000, 7000}}};
  auto a = make_area(rings, {2U, 5U, 6U, 7U});
  a.n_routable_rings_ = 1U;
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  EXPECT_NEAR(dist(a, {{3000, 3000}, {7000, 3000}}),
              result_distances(a, r, true)[6][7], 1e-3);

  // The inner ring of a multipolygon is a way: no generated way along it.
  auto b = make_area(rings, {2U, 5U, 6U, 7U});
  auto const rb = compute_area_ways(b);
  ASSERT_EQ(area_routing_status::kOk, rb.status_);
  for (auto const& w : way_ids(b, rb)) {
    EXPECT_FALSE((w.front() == 6U && w.back() == 7U) ||
                 (w.front() == 7U && w.back() == 6U));
  }
}

TEST(area_routing, overlapping_holes) {
  // Two buildings that overlap: no path through the overlap.
  auto a = make_area(
      {west_east_square_ring(),
       {{3000, 4000}, {6000, 4000}, {6000, 7000}, {3000, 7000}},
       {{5000, 3000}, {8000, 3000}, {8000, 6000}, {5000, 6000}}},
      {2U, 5U});
  a.n_routable_rings_ = 1U;
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto const d = result_distances(a, r, false)[5][2];
  auto const expected =
      std::min(dist(a, {kW, {3000, 7000}, {6000, 7000}, {8000, 6000}, kE}),
               dist(a, {kW, {3000, 4000}, {5000, 3000}, {8000, 3000}, kE}));
  EXPECT_GE(d, expected - 1e-3);
  EXPECT_LE(d, expected * 1.02 + 1e-3);
}

TEST(area_routing, building_crossing_the_outline) {
  // The building sticks out of the north side and blocks the straight line
  // between the entries 4 and 5 inside the area: around it in the south.
  auto a = make_area(
      {{{0, 0}, {10000, 0}, {10000, 10000}, {0, 10000}},
       {{4000, 7000}, {6000, 7000}, {6000, 12000}, {4000, 12000}}},
      {0U}, {{2000, 9000}, {8000, 9000}});
  a.n_routable_rings_ = 1U;
  auto const r = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, r.status_);
  auto const d = result_distances(a, r, true)[8][9];
  auto const expected =
      dist(a, {{2000, 9000}, {4000, 7000}, {6000, 7000}, {8000, 9000}});
  EXPECT_GE(d, expected - 1e-3);
  EXPECT_LE(d, expected * 1.02 + 1e-3);
}

