#include "osr/extract/area_routing.h"

#include <cmath>
#include <algorithm>
#include <limits>
#include <numbers>
#include <queue>
#include <tuple>

#include "osr/extract/point_in_polygon.h"
#include "osr/types.h"

namespace osr {

namespace {

using coord_t = std::int64_t;

// Coordinates are made relative to the first vertex. With an extent below
// 2^24 fixed-point units (1.6 degrees) all cross products fit into 2^50 and
// can be summed up over a ring without overflow.
constexpr auto kMaxExtent = coord_t{1} << 24U;

constexpr auto kInvalid = std::numeric_limits<std::uint32_t>::max();

// Key of one side of a line vertex that is no entry (see lines_): the sides
// must not share a search node, else paths would cross the line there.
constexpr auto kLineSideBit =
    (std::uint64_t{1U} << 62U) | (std::uint64_t{1U} << 61U);

struct vec2 {
  friend vec2 operator-(vec2 const a, vec2 const b) {
    return {a.x_ - b.x_, a.y_ - b.y_};
  }
  friend bool operator==(vec2 const&, vec2 const&) = default;
  coord_t x_, y_;
};

coord_t cross(vec2 const a, vec2 const b) { return a.x_ * b.y_ - a.y_ * b.x_; }

coord_t dot(vec2 const a, vec2 const b) { return a.x_ * b.x_ + a.y_ * b.y_; }

int sign(coord_t const x) { return x > 0 ? 1 : (x < 0 ? -1 : 0); }

// > 0: c is left of a->b, < 0: right, 0: collinear.
int orient(vec2 const a, vec2 const b, vec2 const c) {
  return sign(cross(b - a, c - a));
}

// A corner of a ring, stored so that the walkable side is on the left when
// walking prev -> corner -> next (outer ring counter-clockwise, holes
// clockwise).
struct corner {
  vec2 pos_;
  std::uint32_t prev_, next_;
  std::uint32_t node_;  // search graph node, shared by corners with equal id
  area_routing_result::vertex_ref ref_;
  bool entry_;
  bool reflex_;
  // Free-standing point inside the area: no outline (prev_ == next_ ==
  // itself), stored after all outline corners.
  bool free_;
  // End of an open line (prev_ == next_): walkable in every direction except
  // along the line, paths bend around it.
  bool end_{false};
};

// Is the direction d, seen from corner c, pointing into the walkable side?
bool points_inside(std::vector<corner> const& corners,
                   corner const& c,
                   vec2 const d) {
  if (c.free_) {
    return true;
  }
  auto const a = corners[c.next_].pos_ - c.pos_;
  auto const b = corners[c.prev_].pos_ - c.pos_;
  if (c.end_) {
    return cross(a, d) != 0 || dot(a, d) < 0;
  }
  if (c.prev_ == c.next_ && a == vec2{0, 0}) {
    return true;  // single vertex inside the area
  }
  auto const turn = cross(a, b);  // > 0: convex corner
  if (turn > 0) {
    return cross(a, d) > 0 && cross(d, b) > 0;
  } else if (turn < 0) {
    // Not within the (closed) blocked wedge between prev and next.
    return !(cross(b, d) >= 0 && cross(d, a) >= 0);
  } else {
    // Straight continuation. Spikes (a and b parallel) have no inside.
    return dot(a, b) < 0 && cross(a, d) > 0;
  }
}

// Does a line through the reflex corner c in direction d touch the corner
// without cutting it? Only then a shortest path can bend at c.
bool is_tangent(std::vector<corner> const& corners,
                corner const& c,
                vec2 const d) {
  auto const s1 = sign(cross(d, corners[c.next_].pos_ - c.pos_));
  auto const s2 = sign(cross(d, corners[c.prev_].pos_ - c.pos_));
  return s1 * s2 >= 0;
}

// Is the open segment p-q free of outline crossings and outline vertices?
// Only the first n_outline corners are part of the outline.
bool is_unobstructed(std::vector<corner> const& corners,
                     std::uint32_t const n_outline,
                     std::vector<int>& side,
                     vec2 const p,
                     vec2 const q) {
  auto const d = q - p;
  auto const len2 = dot(d, d);
  for (auto i = 0U; i != n_outline; ++i) {
    auto const& c = corners[i];
    auto const s = orient(p, q, c.pos_);
    side[i] = s;
    if (s == 0) {
      auto const t = dot(c.pos_ - p, d);
      if (t > 0 && t < len2) {
        return false;  // outline vertex on the segment
      }
    }
  }
  for (auto i = 0U; i != n_outline; ++i) {
    auto const j = corners[i].next_;
    if (side[i] * side[j] < 0 &&
        orient(corners[i].pos_, corners[j].pos_, p) *
                orient(corners[i].pos_, corners[j].pos_, q) <
            0) {
      return false;  // proper crossing
    }
  }
  return true;
}

struct edge {
  std::uint32_t to_;
  float dist_;
  bool ring_;
};

struct search_graph {
  void add(std::uint32_t const a,
           std::uint32_t const b,
           float const dist,
           bool const ring) {
    out_[a].push_back({b, dist, ring});
    out_[b].push_back({a, dist, ring});
  }

  std::vector<std::vector<edge>> out_;
};

struct dijkstra {
  struct label {
    friend bool operator>(label const& a, label const& b) {
      return std::tie(a.dist_, a.node_) > std::tie(b.dist_, b.node_);
    }
    double dist_;
    std::uint32_t node_;
  };

  // All distances from start.
  void run(search_graph const& g, std::uint32_t const start) {
    if (dist_.size() != g.out_.size()) {
      dist_.assign(g.out_.size(), std::numeric_limits<double>::max());
      pred_.assign(g.out_.size(), kInvalid);
    } else {
      for (auto const n : touched_) {
        dist_[n] = std::numeric_limits<double>::max();
        pred_[n] = kInvalid;
      }
    }
    touched_.clear();
    pq_ = {};
    dist_[start] = 0.0;
    touched_.push_back(start);
    pq_.push({0.0, start});
    while (!pq_.empty()) {
      auto const [d, u] = pq_.top();
      pq_.pop();
      if (d > dist_[u]) {
        continue;
      }
      for (auto const& e : g.out_[u]) {
        auto const nd = d + static_cast<double>(e.dist_);
        if (nd < dist_[e.to_]) {
          touched_.push_back(e.to_);
          dist_[e.to_] = nd;
          pred_[e.to_] = u;
          pq_.push({nd, e.to_});
        }
      }
    }
  }

  std::vector<double> dist_;
  std::vector<std::uint32_t> pred_;
  std::vector<std::uint32_t> touched_;
  std::priority_queue<label, std::vector<label>, std::greater<>> pq_;
};

}  // namespace

area_routing_result compute_area_ways(area const& a,
                                      area_routing_options const& opt) {
  auto result = area_routing_result{};
  auto const fail = [&](area_routing_status const s) {
    result.status_ = s;
    return result;
  };

  if (a.rings_.empty() || a.rings_.front().size() < 3U) {
    return fail(area_routing_status::kInvalid);
  }

  // Build corners with the walkable side on the left.
  auto const origin = a.rings_.front().front().pos_;
  auto corners = std::vector<corner>{};
  auto ids = std::vector<std::uint64_t>{};
  struct ring_range {
    std::uint32_t base_, n_;
    bool routable_;  // part of the routing graph (see ring_routable_)
  };
  auto ring_ranges = std::vector<ring_range>{};
  for (auto r = 0U; r != a.rings_.size(); ++r) {
    auto const& ring = a.rings_[r];
    auto const n = static_cast<std::uint32_t>(ring.size());
    auto const is_point = r != 0U && n == 1U && ring.front().entry_;
    if (n < 3U && !is_point) {
      if (r == 0U) {
        return fail(area_routing_status::kInvalid);
      }
      continue;
    }

    auto const base = static_cast<std::uint32_t>(corners.size());
    auto doubled_area = coord_t{0};
    for (auto i = 0U; i != n; ++i) {
      auto const pos =
          vec2{static_cast<coord_t>(ring[i].pos_.lat_) - origin.lat_,
               static_cast<coord_t>(ring[i].pos_.lng_) - origin.lng_};
      if (std::abs(pos.x_) >= kMaxExtent || std::abs(pos.y_) >= kMaxExtent) {
        return fail(area_routing_status::kTooLarge);
      }
      corners.push_back({.pos_ = pos,
                         .prev_ = base + (i + n - 1U) % n,
                         .next_ = base + (i + 1U) % n,
                         .node_ = kInvalid,
                         .ref_ = {r, i},
                         .entry_ = ring[i].entry_,
                         .reflex_ = false,
                         .free_ = false});
      ids.push_back(ring[i].id_);
    }
    for (auto i = base; i != base + n; ++i) {
      doubled_area += cross(corners[i].pos_, corners[corners[i].next_].pos_);
    }
    if (is_point) {
      // Its own neighbor: no outline segment, never reflex.
      ring_ranges.push_back({base, n, false});
      continue;
    }
    if (doubled_area == 0) {
      if (r == 0U) {
        return fail(area_routing_status::kInvalid);
      }
      corners.resize(base);
      ids.resize(base);
      continue;
    }

    // Outer ring counter-clockwise, holes clockwise.
    if ((doubled_area > 0) != (r == 0U)) {
      for (auto i = base; i != base + n; ++i) {
        std::swap(corners[i].prev_, corners[i].next_);
      }
    }
    ring_ranges.push_back({base, n, r < a.n_routable_rings_});
  }

  // Lines (barriers), walkable on both sides: an open line becomes a ring
  // around it (there and back) whose ends are spikes, a closed line two
  // rings, one per side.
  auto const n_rings = static_cast<std::uint32_t>(a.rings_.size());
  auto const same_pos = [](area::vertex const& x, area::vertex const& y) {
    return x.pos_.lat_ == y.pos_.lat_ && x.pos_.lng_ == y.pos_.lng_;
  };
  // An end of a line is only walked around if it is free: strictly inside
  // the area and not on another line. A fence that ends at the outline, at
  // a wall or at another fence leaves no gap.
  using line_position = std::tuple<std::int32_t, std::int32_t, std::uint32_t>;
  auto line_positions = std::vector<line_position>{};
  for (auto l = 0U; l != a.lines_.size(); ++l) {
    for (auto const& v : a.lines_[l]) {
      line_positions.emplace_back(v.pos_.lat_, v.pos_.lng_, l);
    }
  }
  std::ranges::sort(line_positions);
  auto const is_free_end = [&](area::vertex const& v, std::uint32_t const l) {
    auto const [lo, hi] = std::ranges::equal_range(
        line_positions, std::pair{v.pos_.lat_, v.pos_.lng_}, {},
        [](line_position const& x) {
          return std::pair{std::get<0>(x), std::get<1>(x)};
        });
    return std::all_of(
               lo, hi,
               [&](line_position const& x) { return std::get<2>(x) == l; }) &&
           locate_in_rings(a.rings_, v.pos_, &area::vertex::pos_) ==
               point_location::kInside;
  };

  auto seq = std::vector<std::uint32_t>{};
  auto line_ids = std::vector<std::uint64_t>{};
  for (auto l = 0U; l != a.lines_.size(); ++l) {
    auto const& line = a.lines_[l];
    seq.clear();
    for (auto i = 0U; i != line.size(); ++i) {
      if (seq.empty() || !same_pos(line[i], line[seq.back()])) {
        seq.push_back(i);
      }
    }
    auto const closed =
        seq.size() > 3U && same_pos(line[seq.front()], line[seq.back()]);
    if (closed) {
      seq.pop_back();
    }
    if (seq.size() < 2U) {
      continue;
    }

    auto const add_ring = [&](bool const is_open) {
      auto const base = static_cast<std::uint32_t>(corners.size());
      auto const n = static_cast<std::uint32_t>(seq.size());
      auto const last_end = (n + 2U) / 2U - 1U;  // index of the far end
      for (auto k = 0U; k != n; ++k) {
        auto const& v = line[seq[k]];
        auto const pos = vec2{static_cast<coord_t>(v.pos_.lat_) - origin.lat_,
                              static_cast<coord_t>(v.pos_.lng_) - origin.lng_};
        if (std::abs(pos.x_) >= kMaxExtent || std::abs(pos.y_) >= kMaxExtent) {
          return false;
        }
        corners.push_back({.pos_ = pos,
                           .prev_ = base + (k + n - 1U) % n,
                           .next_ = base + (k + 1U) % n,
                           .node_ = kInvalid,
                           .ref_ = {n_rings + l, seq[k]},
                           .entry_ = v.entry_,
                           .reflex_ = false,
                           .free_ = false,
                           .end_ = is_open && (k == 0U || k == last_end) &&
                                   is_free_end(v, l)});
        ids.push_back(v.entry_ && !v.blocked_
                          ? v.id_
                          : kLineSideBit | (corners.size() - 1U));
        line_ids.push_back(v.id_);
      }
      ring_ranges.push_back({base, n, false});
      return true;
    };
    if (closed) {
      auto ok = add_ring(false);
      std::ranges::reverse(seq);
      ok = ok && add_ring(false);
      if (!ok) {
        return fail(area_routing_status::kTooLarge);
      }
    } else {
      for (auto i = seq.size() - 1U; i-- > 1U;) {
        seq.push_back(seq[i]);
      }
      if (!add_ring(true)) {
        return fail(area_routing_status::kTooLarge);
      }
    }
  }

  // Entries that are not on the outline: point corners behind the outline
  // corners. They have to be strictly inside (exact test), otherwise a line
  // of sight could run along or outside of the outline.
  auto const n_outline = static_cast<std::uint32_t>(corners.size());
  if (!a.interior_.empty()) {
    auto outline_ids = ids;
    outline_ids.insert(end(outline_ids), begin(line_ids), end(line_ids));
    std::ranges::sort(outline_ids);
    for (auto i = 0U; i != a.interior_.size(); ++i) {
      auto const& v = a.interior_[i];
      auto const pos = vec2{static_cast<coord_t>(v.pos_.lat_) - origin.lat_,
                            static_cast<coord_t>(v.pos_.lng_) - origin.lng_};
      if (std::abs(pos.x_) >= kMaxExtent || std::abs(pos.y_) >= kMaxExtent ||
          std::ranges::binary_search(outline_ids, v.id_) ||
          locate_in_rings(a.rings_, v.pos_, &area::vertex::pos_) !=
              point_location::kInside) {
        continue;
      }
      auto const self = static_cast<std::uint32_t>(corners.size());
      corners.push_back(
          {.pos_ = pos,
           .prev_ = self,
           .next_ = self,
           .node_ = kInvalid,
           .ref_ = {area_routing_result::vertex_ref::kInterior, i},
           .entry_ = true,
           .reflex_ = false,
           .free_ = true});
      ids.push_back(v.id_);
    }
  }

  for (auto& c : corners) {
    c.reflex_ = c.end_ || cross(c.pos_ - corners[c.prev_].pos_,
                                corners[c.next_].pos_ - c.pos_) < 0;
  }

  // Search graph nodes: one per distinct vertex id.
  auto n_nodes = 0U;
  auto node_corner = std::vector<std::uint32_t>{};  // representative corner
  auto node_ids = std::vector<std::uint64_t>{};  // ascending
  {
    auto order = std::vector<std::uint32_t>(corners.size());
    for (auto i = 0U; i != order.size(); ++i) {
      order[i] = i;
    }
    std::ranges::stable_sort(order, {}, [&](auto const i) { return ids[i]; });
    for (auto i = 0U; i != order.size(); ++i) {
      if (i == 0U || ids[order[i]] != ids[order[i - 1U]]) {
        node_corner.push_back(order[i]);
        node_ids.push_back(ids[order[i]]);
        ++n_nodes;
      }
      corners[order[i]].node_ = n_nodes - 1U;
    }
  }

  // Nothing passes a blocked vertex: it is not an entry.
  auto node_is_entry = std::vector<bool>(n_nodes);
  for (auto const& c : corners) {
    if (c.entry_) {
      node_is_entry[c.node_] = true;
    }
  }
  for (auto const& c : corners) {
    if (get_vertex(a, c.ref_).blocked_) {
      node_is_entry[c.node_] = false;
    }
  }
  auto entries = std::vector<std::uint32_t>{};
  for (auto i = 0U; i != n_nodes; ++i) {
    if (node_is_entry[i]) {
      entries.push_back(i);
      result.n_interior_entries_ += node_corner[i] >= n_outline ? 1U : 0U;
    }
  }
  result.n_entries_ = static_cast<std::uint32_t>(entries.size());
  if (entries.size() < 2U) {
    return fail(area_routing_status::kNoEntries);
  }
  if (entries.size() > opt.max_entries_) {
    return fail(area_routing_status::kTooLarge);
  }

  auto candidates = std::vector<std::uint32_t>{};
  for (auto i = 0U; i != corners.size(); ++i) {
    if (corners[i].reflex_ || node_is_entry[corners[i].node_]) {
      candidates.push_back(i);
    }
  }
  if (candidates.size() > opt.max_candidates_) {
    return fail(area_routing_status::kTooLarge);
  }

  // Metric scaling. Fixed-point coordinates are 1e-7 degrees.
  constexpr auto kMetersPerUnit = 6371000.0 * std::numbers::pi / 180.0 * 1e-7;
  auto const x_scale =
      kMetersPerUnit * std::cos(origin.as_latlng().lat_ * std::numbers::pi /
                                180.0);
  auto const y_scale = kMetersPerUnit;
  auto const distance = [&](vec2 const p, vec2 const q) {
    auto const d = q - p;
    return static_cast<float>(std::hypot(static_cast<double>(d.x_) * x_scale,
                                         static_cast<double>(d.y_) * y_scale));
  };

  auto g = search_graph{};
  g.out_.resize(n_nodes);

  // Outline.
  for (auto const& c : corners) {
    auto const& next = corners[c.next_];
    if (c.node_ != next.node_) {
      g.add(c.node_, next.node_, distance(c.pos_, next.pos_), true);
    }
  }

  // Tangent lines of sight.
  auto side = std::vector<int>(corners.size());
  for (auto i = 0U; i != candidates.size(); ++i) {
    auto const& ci = corners[candidates[i]];
    auto const i_is_entry = node_is_entry[ci.node_];
    for (auto j = i + 1U; j != candidates.size(); ++j) {
      auto const& cj = corners[candidates[j]];
      if (ci.node_ == cj.node_ || ci.pos_ == cj.pos_) {
        continue;
      }
      auto const d = cj.pos_ - ci.pos_;
      auto const rd = ci.pos_ - cj.pos_;
      if (!points_inside(corners, ci, d) || !points_inside(corners, cj, rd)) {
        continue;
      }
      // Paths end at entries in any direction, but only bend at tangents.
      auto const j_is_entry = node_is_entry[cj.node_];
      if ((!i_is_entry && !is_tangent(corners, ci, d)) ||
          (!j_is_entry && !is_tangent(corners, cj, rd))) {
        continue;
      }
      if (!is_unobstructed(corners, n_outline, side, ci.pos_, cj.pos_)) {
        continue;
      }
      g.add(ci.node_, cj.node_, distance(ci.pos_, cj.pos_), false);
      ++result.n_visibility_edges_;
    }
  }

  // Shortest paths from every entry to all others.
  auto const k = entries.size();
  auto entry_idx = std::vector<std::uint32_t>(n_nodes, kInvalid);
  for (auto i = 0U; i != k; ++i) {
    entry_idx[entries[i]] = i;
  }

  struct entry_pair {
    double dist_;
    std::uint32_t from_, to_;
  };
  auto pairs = std::vector<entry_pair>{};
  auto trees = std::vector<std::vector<std::uint32_t>>(k);
  auto d = dijkstra{};
  for (auto si = 0U; si + 1U < k; ++si) {
    d.run(g, entries[si]);
    for (auto ti = si + 1U; ti != k; ++ti) {
      if (d.pred_[entries[ti]] != kInvalid) {  // else: invalid geometry
        pairs.push_back({d.dist_[entries[ti]], si, ti});
      }
    }
    trees[si] = d.pred_;
  }
  std::ranges::sort(pairs, [](entry_pair const& x, entry_pair const& y) {
    return std::tie(x.dist_, x.from_, x.to_) <
           std::tie(y.dist_, y.from_, y.to_);
  });

  // Greedy spanner: going through the pairs from short to long, the
  // shortest path of a pair is only added if the paths added so far (and
  // the outline) do not already connect it within the allowed stretch.
  // Paths share their segments: where they meet or part at a corner, the
  // corner becomes a junction. This drops all paths that follow the
  // outline, and the parts of paths that are covered by others already.
  // h: outline between consecutive entries, links, added segments.
  auto h = search_graph{};
  h.out_.resize(n_nodes);
  if (opt.ring_routable_) {
    auto const is_blocked = [&](std::uint32_t const i) {
      return get_vertex(a, corners[i].ref_).blocked_;
    };
    for (auto const& [base, n, routable] : ring_ranges) {
      if (n < 2U || !routable) {
        continue;
      }
      // Walk the ring once, starting at a blocked vertex if there is one:
      // consecutive entries are connected unless a blocked vertex lies
      // between them (a blocked entry is not connected at all).
      auto start = base;
      for (auto i = base; i != base + n; ++i) {
        if (is_blocked(i)) {
          start = i;
          break;
        }
      }
      auto const has_block = is_blocked(start);
      auto first = kInvalid, prev = kInvalid;
      auto first_at = 0.0, prev_at = 0.0, length = 0.0;
      for (auto step = 0U; step != n; ++step) {
        auto const i = base + (start - base + step) % n;
        if (is_blocked(i)) {
          prev = kInvalid;
        } else if (auto const e = entry_idx[corners[i].node_]; e != kInvalid) {
          if (prev != kInvalid && prev != e) {
            h.add(entries[prev], entries[e], static_cast<float>(length - prev_at),
                  true);
          }
          if (first == kInvalid) {
            first = e;
            first_at = length;
          }
          prev = e;
          prev_at = length;
        }
        auto const next = i + 1U == base + n ? base : i + 1U;
        length += static_cast<double>(
            distance(corners[i].pos_, corners[next].pos_));
      }
      // Closing the ring: from the last entry back to the first one.
      if (!has_block && prev != kInvalid && first != kInvalid &&
          prev != first) {
        h.add(entries[prev], entries[first],
              static_cast<float>(length - prev_at + first_at), true);
      }
    }
  }

  // Connections between entries that exist already.
  auto const find_entry = [&](std::uint64_t const id) {
    auto const it = std::ranges::lower_bound(node_ids, id);
    return it != end(node_ids) && *it == id
               ? entry_idx[static_cast<std::size_t>(it - begin(node_ids))]
               : kInvalid;
  };
  for (auto const& l : a.links_) {
    auto const from = find_entry(l.from_);
    auto const to = find_entry(l.to_);
    if (from != kInvalid && to != kInvalid && from != to) {
      h.add(entries[from], entries[to], static_cast<float>(l.length_), true);
    }
  }

  auto const is_added = [&](std::uint32_t const x, std::uint32_t const y) {
    return std::ranges::any_of(
        h.out_[x], [&](edge const& e) { return e.to_ == y && !e.ring_; });
  };
  auto degree = std::vector<std::uint32_t>(n_nodes);
  auto path = std::vector<std::uint32_t>{};
  auto in_h = dijkstra{};
  if (opt.all_visibility_edges_) {
    for (auto x = 0U; x != n_nodes; ++x) {
      for (auto const& e : g.out_[x]) {
        if (x < e.to_ && !is_added(x, e.to_)) {
          ++degree[x];
          ++degree[e.to_];
          h.add(x, e.to_, e.dist_, false);
        }
      }
    }
    pairs.clear();
  }
  // Distances in h only get shorter as segments are added: distances from
  // an entry computed earlier are upper bounds. They are only computed
  // again when the bound is not good enough (Farshi and Gudmundsson).
  auto upper = std::vector<std::vector<double>>(k);
  for (auto const& p : pairs) {
    // Tolerance: edge lengths are stored as float.
    auto const bound = p.dist_ * opt.max_stretch_ * (1.0 + 1e-5) + 1e-3;
    auto& from = upper[p.from_];
    auto const& to = upper[p.to_];
    if ((!from.empty() && from[entries[p.to_]] <= bound) ||
        (!to.empty() && to[entries[p.from_]] <= bound)) {
      continue;
    }
    in_h.run(h, entries[p.from_]);
    from = in_h.dist_;
    if (from[entries[p.to_]] <= bound) {
      continue;
    }

    // New segments of the path, and the degrees they lead to.
    path.clear();
    for (auto n = entries[p.to_]; n != kInvalid; n = trees[p.from_][n]) {
      path.push_back(n);
    }
    auto over = false;
    for (auto i = 0U; i + 1U < path.size(); ++i) {
      if (!is_added(path[i], path[i + 1U])) {
        over = over || degree[path[i]] + 1U > opt.max_ways_per_node_ ||
               degree[path[i + 1U]] + 1U > opt.max_ways_per_node_;
      }
    }
    if (over) {
      ++result.n_capped_;
      continue;
    }
    for (auto i = 0U; i + 1U < path.size(); ++i) {
      auto const x = path[i], y = path[i + 1U];
      if (!is_added(x, y)) {
        ++degree[x];
        ++degree[y];
        h.add(x, y,
              distance(corners[node_corner[x]].pos_,
                       corners[node_corner[y]].pos_),
              false);
      }
    }
  }

  // Ways: chains of added segments between entries and junctions (corners
  // where more than two segments meet).
  auto const is_end = [&](std::uint32_t const n) {
    return node_is_entry[n] || degree[n] != 2U;
  };
  auto done = std::vector<std::pair<std::uint32_t, std::uint32_t>>{};
  auto const is_done = [&](std::uint32_t const x, std::uint32_t const y) {
    return std::ranges::binary_search(done, std::pair{std::min(x, y),
                                                      std::max(x, y)});
  };
  for (auto start = 0U; start != n_nodes; ++start) {
    if (degree[start] == 0U || !is_end(start)) {
      continue;
    }
    auto next = std::vector<std::uint32_t>{};
    for (auto const& e : h.out_[start]) {
      if (!e.ring_) {
        next.push_back(e.to_);
      }
    }
    std::ranges::sort(next);
    for (auto const first : next) {
      if (is_done(start, first)) {
        continue;
      }
      auto& way = result.ways_.emplace_back();
      way.push_back(corners[node_corner[start]].ref_);
      auto prev = start, curr = first;
      auto chain = std::vector<std::pair<std::uint32_t, std::uint32_t>>{};
      while (true) {
        chain.emplace_back(std::min(prev, curr), std::max(prev, curr));
        way.push_back(corners[node_corner[curr]].ref_);
        if (is_end(curr)) {
          break;
        }
        for (auto const& e : h.out_[curr]) {
          if (!e.ring_ && e.to_ != prev) {
            prev = curr;
            curr = e.to_;
            break;
          }
        }
      }
      result.end_ids_.emplace_back(node_ids[start], node_ids[curr]);
      for (auto const& c : chain) {
        done.insert(std::ranges::upper_bound(done, c), c);
      }
    }
  }
  for (auto n = 0U; n != n_nodes; ++n) {
    if (degree[n] != 0U && !node_is_entry[n] && degree[n] != 2U) {
      result.junctions_.push_back(node_ids[n]);
    }
  }

  return result;
}

}  // namespace osr
