#include "osr/extract/area_ways.h"

#include <algorithm>
#include <array>
#include <limits>
#include <mutex>
#include <tuple>

#include "rtree.h"

#include "geo/latlng.h"

#include "utl/parallel_for.h"
#include "utl/verify.h"

#include "osr/extract/point_in_polygon.h"
#include "osr/ways.h"

namespace osr {

namespace {

// Same limit as in compute_area_ways: larger areas are skipped there. Keeps
// all coordinate differences small enough for the exact predicates.
constexpr auto kMaxExtent = std::int64_t{1} << 24U;

struct box {
  bool contains(point const p) const {
    return p.lat_ >= min_.lat_ && p.lat_ <= max_.lat_ &&
           p.lng_ >= min_.lng_ && p.lng_ <= max_.lng_;
  }

  void extend(point const p) {
    min_ = {std::min(min_.lat_, p.lat_), std::min(min_.lng_, p.lng_)};
    max_ = {std::max(max_.lat_, p.lat_), std::max(max_.lng_, p.lng_)};
  }

  std::array<double, 2> min() const {
    return {static_cast<double>(min_.lat_), static_cast<double>(min_.lng_)};
  }

  std::array<double, 2> max() const {
    return {static_cast<double>(max_.lat_), static_cast<double>(max_.lng_)};
  }

  point min_{std::numeric_limits<std::int32_t>::max(),
             std::numeric_limits<std::int32_t>::max()};
  point max_{std::numeric_limits<std::int32_t>::min(),
             std::numeric_limits<std::int32_t>::min()};
};

// Ways without level are on the ground.
level_t ground(level_t const l) {
  return l == kNoLevel ? level_t{0.F} : l;
}

}  // namespace

// Entries that are not on the outline: nodes of walkable ways strictly inside
// an area. A walker can leave the way there and cross the area freely.
//
// Not every node inside an area becomes an entry:
//  - Ends of ways (a footway that ends in the middle of a plaza) and nodes
//    that are routing nodes anyway (two ways meet, barriers) are entries.
//    The latter are optional (interior_entries): they do not cost routing
//    nodes, but most of the additional ways.
//  - All other nodes are not. They are only the shape of a way that enters
//    the area through one of the entries above or through an outline node,
//    and leaving the way later is never shorter than leaving it right there.
//    Making them entries would turn every one of them into a routing node.
//  - The first node of a closed way (a platform or a loop inside the area)
//    is neither: it is a routing node only because it is counted twice, and
//    an arbitrary place to connect the loop to the area. If another way
//    starts there, the node is found as the end of that way.
//  - Nodes on an outline segment that are not outline nodes are not mapped
//    as connected to the area. Where a way crosses the outline without a
//    common node, nothing is connected either - but if it ends inside, its
//    end is an entry like any other.
//
// What lies inside the outline on the map is not necessarily on the area:
// ways on another level, in a tunnel or on a bridge are skipped, and so are
// steps and elevators (one of their ends is on another level).
//
// One r-tree lookup per way with the bounding boxes of all areas, then the
// exact point-in-polygon test for the nodes of the few ways that are close to
// an area.
area_ways::interior area_ways::find_interior(
    std::vector<bool> const& is_outline) const {
  auto const n_way_areas = way_areas_.size();
  auto const n_areas = n_way_areas + relation_areas_.size();
  auto const is_off_ground = [&](way_idx_t const way) {
    return to_idx(way) < off_ground_.size() && off_ground_[to_idx(way)];
  };

  auto result = interior{};
  if (n_areas == 0U || interior_ == interior_entries::kNone) {
    return result;
  }

  // Areas that can have interior entries, with level and bounding box.
  struct indexed_area {
    box box_;
    level_t level_;
  };
  auto areas = std::vector<indexed_area>(n_areas);
  auto* const index = rtree_new();
  utl::verify(index != nullptr, "area rtree creation failed");
  for (auto i = std::size_t{0U}; i != n_areas; ++i) {
    auto& a = areas[i];
    auto p = way_properties{};
    if (i < n_way_areas) {
      auto const way = way_areas_[i];
      p = w_.r_->way_properties_[way];
      if (is_off_ground(way)) {
        continue;
      }
      for (auto const pos : w_.way_polylines_[way]) {
        a.box_.extend(pos);
      }
    } else {
      auto const& rel = relation_areas_[i - n_way_areas];
      p = rel.p_;
      if (rel.off_ground_ || rel.rings_.empty()) {
        continue;
      }
      for (auto const& [node, pos] : rel.rings_.front()) {
        a.box_.extend(pos);
      }
    }
    a.level_ = ground(p.from_level());
    if (p.from_level_ != p.to_level_ ||
        a.box_.min_.lat_ > a.box_.max_.lat_ ||
        std::int64_t{a.box_.max_.lat_} - a.box_.min_.lat_ >= kMaxExtent ||
        std::int64_t{a.box_.max_.lng_} - a.box_.min_.lng_ >= kMaxExtent) {
      continue;
    }
    rtree_insert(index, a.box_.min().data(), a.box_.max().data(),
                 reinterpret_cast<void*>(i));
  }

  auto const locate = [&](std::size_t const i, point const p) {
    if (i < n_way_areas) {
      return locate_in_ring(w_.way_polylines_[way_areas_[i]], p);
    } else {
      return locate_in_rings(relation_areas_[i - n_way_areas].rings_, p,
                             [](auto const& x) { return x.second; });
    }
  };

  // Proper crossings of the way with the outline segments of area i: the
  // way does not use the outline nodes there.
  auto const count_crossings = [&](std::size_t const i,
                                   auto const& polyline) {
    using wide = __int128;
    auto const orient = [](point const a, point const b, point const c) {
      auto const x = static_cast<wide>(std::int64_t{b.lat_} - a.lat_) *
                         (std::int64_t{c.lng_} - a.lng_) -
                     static_cast<wide>(std::int64_t{b.lng_} - a.lng_) *
                         (std::int64_t{c.lat_} - a.lat_);
      return x > 0 ? 1 : (x < 0 ? -1 : 0);
    };
    auto n = std::uint64_t{0U};
    auto const check = [&](point const a, point const b) {
      for (auto j = 1U; j < polyline.size(); ++j) {
        auto const p = point{polyline[j - 1U]};
        auto const q = point{polyline[j]};
        if (orient(a, b, p) * orient(a, b, q) < 0 &&
            orient(p, q, a) * orient(p, q, b) < 0) {
          ++n;
        }
      }
    };
    if (i < n_way_areas) {
      auto const ring = w_.way_polylines_[way_areas_[i]];
      for (auto j = 1U; j < ring.size(); ++j) {
        check(point{ring[j - 1U]}, point{ring[j]});
      }
    } else {
      for (auto const& r : relation_areas_[i - n_way_areas].rings_) {
        for (auto j = 0U; j != r.size(); ++j) {
          check(r[j].second, r[(j + 1U) % r.size()].second);
        }
      }
    }
    return n;
  };

  struct scratch {
    std::uint64_t n_crossings_{0U};
    std::vector<std::size_t> areas_;
    std::vector<point_location> location_;
    std::vector<bool> is_entry_;
    std::vector<interior::entry> entries_;
    std::vector<interior::link> links_;
  };
  auto mutex = std::mutex{};
  utl::parallel_for_run_threadlocal<scratch>(
      w_.n_ways(), [&](scratch& s, std::size_t const idx) {
        auto const way = way_idx_t{static_cast<way_idx_t::value_t>(idx)};
        auto const p = w_.r_->way_properties_[way];
        if (is_outline[idx] || is_off_ground(way) ||
            !(p.is_foot_accessible() || p.is_bike_accessible()) ||
            p.is_steps() || p.is_elevator() || p.from_level_ != p.to_level_) {
          return;
        }

        auto const polyline = w_.way_polylines_[way];
        if (polyline.size() < 2U) {
          return;
        }
        auto way_box = box{};
        for (auto const pos : polyline) {
          way_box.extend(pos);
        }

        s.areas_.clear();
        rtree_search(
            index, way_box.min().data(), way_box.max().data(),
            [](auto, auto, void const* item, void* udata) {
              static_cast<std::vector<std::size_t>*>(udata)->push_back(
                  reinterpret_cast<std::size_t>(item));
              return true;
            },
            &s.areas_);
        if (s.areas_.empty()) {
          return;
        }

        auto const nodes = w_.way_osm_nodes_[way];
        auto const n = polyline.size();
        auto const is_closed = nodes.front() == nodes.back();
        s.entries_.clear();
        s.links_.clear();
        for (auto const i : s.areas_) {
          auto const& a = areas[i];
          if (a.level_ != ground(p.from_level())) {
            continue;
          }

          auto any = false;
          s.location_.assign(n, point_location::kOutside);
          s.is_entry_.assign(n, false);
          for (auto j = 0U; j != n; ++j) {
            auto const pos = point{polyline[j]};
            if (!a.box_.contains(pos)) {
              continue;
            }
            s.location_[j] = locate(i, pos);
            if (s.location_[j] != point_location::kInside) {
              continue;
            }
            auto const is_end = j == 0U || j + 1U == n;
            auto const is_shared =
                w_.node_way_counter_.is_multi(to_idx(nodes[j]));
            if (interior_ == interior_entries::kAllNodes ||
                ((is_end ? !is_closed : is_shared) &&
                 (!is_shared ||
                  interior_ == interior_entries::kRoutingNodes))) {
              s.is_entry_[j] = true;
              any = true;
            }
          }
          s.n_crossings_ += count_crossings(i, polyline);
          if (!any) {
            continue;
          }

          // The way connects each entry to the previous and next entry on
          // it, there is no need to generate a way for these pairs. Nodes on
          // the outline can be entries as well (if they are outline nodes).
          auto last = n;
          auto length = 0.0;
          for (auto j = 0U; j != n; ++j) {
            if (j != 0U) {
              length += geo::distance(point{polyline[j - 1U]}.as_latlng(),
                                      point{polyline[j]}.as_latlng());
            }
            if (s.is_entry_[j]) {
              s.entries_.push_back({static_cast<std::uint32_t>(i), nodes[j],
                                    point{polyline[j]}});
            } else if (s.location_[j] != point_location::kBoundary) {
              continue;
            }
            if (last != n && (s.is_entry_[last] || s.is_entry_[j]) &&
                nodes[last] != nodes[j] && p.is_foot_accessible()) {
              s.links_.push_back(
                  {static_cast<std::uint32_t>(i),
                   std::min(nodes[last], nodes[j]),
                   std::max(nodes[last], nodes[j]), length});
            }
            last = j;
            length = 0.0;
          }
        }

        if (s.n_crossings_ != 0U) {
          auto const lock = std::scoped_lock{mutex};
          result.n_crossings_ += s.n_crossings_;
          s.n_crossings_ = 0U;
        }
        if (!s.entries_.empty()) {
          auto const lock = std::scoped_lock{mutex};
          result.entries_.insert(end(result.entries_), begin(s.entries_),
                                 end(s.entries_));
          result.links_.insert(end(result.links_), begin(s.links_),
                               end(s.links_));
        }
      });
  rtree_free(index);

  // Deterministic order, independent of the threads. A node is found once
  // for every way it is on.
  std::ranges::sort(result.entries_, {}, [](interior::entry const& x) {
    return std::tuple{x.area_, x.node_};
  });
  auto const [first, last] = std::ranges::unique(
      result.entries_, {},
      [](interior::entry const& x) { return std::tuple{x.area_, x.node_}; });
  result.entries_.erase(first, last);
  std::ranges::sort(result.links_, {}, [](interior::link const& x) {
    return std::tuple{x.area_, x.from_, x.to_, x.length_};
  });
  return result;
}

}  // namespace osr
