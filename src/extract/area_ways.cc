#include "osr/extract/area_ways.h"

#include <cstdlib>
#include <cstring>
#include <algorithm>
#include <fstream>
#include <ranges>
#include <sstream>

#include "rtree.h"

#include "utl/parallel_for.h"
#include "utl/verify.h"

#include "osr/extract/area_merge.h"
#include "osr/extract/point_in_polygon.h"
#include "osr/ways.h"

namespace osr {

namespace {

osm_way_idx_t to_area_way_idx(std::uint64_t const encoded_origin,
                              bool const from_relation) {
  constexpr auto kIdMask = kAreaWayNegativeIdBit - 1U;
  return osm_way_idx_t{
      kAreaWayBit | (from_relation ? kAreaWayFromRelationBit : 0U) |
      (is_negative_encoded_osm_id(encoded_origin) ? kAreaWayNegativeIdBit
                                                  : 0U) |
      (encoded_origin & kIdMask)};
}

// A way crossing an area is walked (or cycled) in any direction: it keeps
// the access and the level of the area, everything that describes the
// outline as a street (oneway, steps, vehicle access, ...) does not apply.
way_properties to_area_way_properties(way_properties const& area) {
  auto p = way_properties{};
  p.is_foot_accessible_ = area.is_foot_accessible_;
  p.is_bike_accessible_ = area.is_bike_accessible_;
  p.motor_vehicle_no_ = area.motor_vehicle_no_;
  p.from_level_ = area.from_level_;
  p.to_level_ = area.to_level_;
  return p;
}

std::uint64_t segment_hash(osm_node_idx_t const a, osm_node_idx_t const b) {
  auto const lo = std::min(to_idx(a), to_idx(b));
  auto const hi = std::max(to_idx(a), to_idx(b));
  auto x = (lo * 0x9E3779B97F4A7C15ULL) ^ hi;
  x = (x ^ (x >> 30U)) * 0xBF58476D1CE4E5B9ULL;
  x = (x ^ (x >> 27U)) * 0x94D049BB133111EBULL;
  return x ^ (x >> 31U);
}

// Areas are only merged if the generated ways get the same properties from
// each of them: a way crossing the union must not open one area for a profile
// that only the other one allows (or drop it where it was allowed). This
// includes the levels. Not compared: motor_vehicle_no_, which only changes
// the costs (the generated ways get it if all merged areas have it).
std::uint64_t to_merge_class(way_properties const& area) {
  auto p = to_area_way_properties(area);
  p.motor_vehicle_no_ = false;
  auto c = std::uint64_t{0U};
  std::memcpy(&c, &p, sizeof(p));
  return c;
}

// Ways without level are on the ground.
level_t ground(level_t const l) { return l == kNoLevel ? level_t{0.F} : l; }

std::array<double, 2> to_array(point const p) {
  return {static_cast<double>(p.lat_), static_cast<double>(p.lng_)};
}

// An area as collected from OSM.
struct input {
  osm_way_idx_t id_;
  way_properties p_;
  string_idx_t name_;
  bool mixed_levels_;
  std::int8_t layer_;
};

struct generated {
  osm_way_idx_t id_;
  way_properties p_;
  string_idx_t name_;
  bool mixed_levels_{false};
  std::int8_t layer_{0};
  std::uint32_t n_holes_{0U}, n_lines_{0U};
  area_routing_status status_{area_routing_status::kOk};
  std::uint64_t n_visibility_edges_{0U};
  std::uint64_t n_capped_{0U};
  std::uint32_t n_entries_{0U}, n_interior_entries_{0U};
  bool interior_dropped_{false};
  std::vector<area_ways::ring> ways_;
  std::vector<std::uint64_t> junctions_;  // vertex ids, sorted
};

// What is processed at once: one area or a group of merged areas. A group
// that turns out to be too large is processed as separate areas.
struct unit {
  std::vector<generated> out_;
  bool too_large_{false};
  std::uint64_t n_on_removed_{0U};
};

}  // namespace

area_ways::area_ways(ways& w,
                     area_routing_options const& opt,
                     interior_entries const interior)
    : w_{w}, opt_{opt}, interior_{interior} {}

area_ways::~area_ways() {
  if (area_index_ != nullptr) {
    rtree_free(static_cast<rtree*>(area_index_));
  }
}

void area_ways::add_way(way_idx_t const way, std::int8_t const layer) {
  way_areas_.push_back(way);
  way_area_layers_.push_back(layer);
}

void area_ways::add_relation(std::int64_t const relation_id,
                             std::vector<ring> rings,
                             std::vector<way_idx_t> member_ways,
                             way_properties const p,
                             bool const off_ground,
                             std::int8_t const layer) {
  relation_areas_.push_back({relation_id, std::move(rings),
                             std::move(member_ways), p, off_ground, layer});
}

void area_ways::index_areas() {
  auto* const index = rtree_new();
  utl::verify(index != nullptr, "area rtree creation failed");
  area_index_ = index;
  auto const insert = [&](auto const& positions) {
    auto min = point{std::numeric_limits<std::int32_t>::max(),
                     std::numeric_limits<std::int32_t>::max()};
    auto max = point{std::numeric_limits<std::int32_t>::min(),
                     std::numeric_limits<std::int32_t>::min()};
    for (point const p : positions) {
      min = {std::min(min.lat_, p.lat_), std::min(min.lng_, p.lng_)};
      max = {std::max(max.lat_, p.lat_), std::max(max.lng_, p.lng_)};
    }
    if (min.lat_ <= max.lat_) {
      rtree_insert(index, to_array(min).data(), to_array(max).data(), nullptr);
    }
  };
  for (auto const way : way_areas_) {
    insert(w_.way_polylines_[way]);
  }
  for (auto const& rel : relation_areas_) {
    if (!rel.rings_.empty()) {
      insert(rel.rings_.front() |
             std::views::transform([](auto const& x) { return x.second; }));
    }
  }
}

bool area_ways::is_near_area(point const min, point const max) const {
  if (area_index_ == nullptr) {
    return false;
  }
  auto found = false;
  rtree_search(
      static_cast<rtree const*>(area_index_), to_array(min).data(),
      to_array(max).data(),
      [](double const*, double const*, void const*, void* udata) {
        *static_cast<bool*>(udata) = true;
        return false;
      },
      &found);
  return found;
}

void area_ways::add_obstacles(std::vector<obstacle>&& obstacles) {
  auto const l = std::scoped_lock{obstacle_mutex_};
  obstacles_.insert(end(obstacles_), std::make_move_iterator(begin(obstacles)),
                    std::make_move_iterator(end(obstacles)));
}

void area_ways::set_off_ground(way_idx_t const way) {
  off_ground_.resize(
      std::max(off_ground_.size(), static_cast<std::size_t>(to_idx(way)) + 1U));
  off_ground_[to_idx(way)] = true;
}

void area_ways::set_connector(way_idx_t const way) {
  connector_.resize(
      std::max(connector_.size(), static_cast<std::size_t>(to_idx(way)) + 1U));
  connector_[to_idx(way)] = true;
}

void area_ways::add_elevator_node(osm_node_idx_t const n) {
  elevator_nodes_.push_back(n);
}

void area_ways::add_barrier(std::span<osm_node_idx_t const> const nodes) {
  auto const l = std::scoped_lock{barrier_mutex_};
  for (auto i = 1U; i < nodes.size(); ++i) {
    barrier_segments_.push_back(segment_hash(nodes[i - 1U], nodes[i]));
  }
}

void area_ways::add_blocked_node(osm_node_idx_t const n) {
  blocked_nodes_.push_back(n);
}

area_ways::stats area_ways::build() {
  // How the nodes of the area outlines are used:
  //  - outline_: number of area outlines the node is part of
  //  - other_: number of other walkable ways - these make a node an entry
  // Nodes that are only shared between outlines are not entries: adjacent
  // areas are connected through their outlines.
  struct node_use {
    std::uint8_t outline_{0U}, other_{0U}, generated_{0U};
  };
  auto const inc = [](std::uint8_t& x) {
    x = static_cast<std::uint8_t>(std::min(254, x + 1));
  };

  auto use = hash_map<osm_node_idx_t, node_use>{};
  auto is_outline = std::vector<bool>(w_.n_ways());
  auto const add_outline = [&](way_idx_t const way) {
    if (is_outline[to_idx(way)]) {
      return;
    }
    is_outline[to_idx(way)] = true;
    auto const nodes = w_.way_osm_nodes_[way];
    auto const closed = nodes.size() > 1U && nodes.front() == nodes.back();
    for (auto i = 0U; i != nodes.size() - (closed ? 1U : 0U); ++i) {
      inc(use[nodes[i]].outline_);
    }
  };
  for (auto const way : way_areas_) {
    add_outline(way);
  }
  for (auto const& rel : relation_areas_) {
    for (auto const way : rel.member_ways_) {
      add_outline(way);
    }
  }

  // Entries that are not on the outline. In `use` for the node budget, their
  // other_ ways are counted below.
  std::ranges::sort(elevator_nodes_);
  auto const inside = find_interior(is_outline);
  auto s = stats{};
  s.n_crossings_ = inside.n_crossings_;
  for (auto const& e : inside.entries_) {
    use[e.node_];
  }

  // Obstacles: deterministic order, nodes can be entries (building
  // entrances, gates).
  std::ranges::sort(obstacles_, {}, &obstacle::id_);
  auto* const obstacle_index = rtree_new();
  utl::verify(obstacle_index != nullptr, "obstacle rtree creation failed");
  for (auto i = std::size_t{0U}; i != obstacles_.size(); ++i) {
    auto min = obstacles_[i].nodes_.front().second;
    auto max = min;
    for (auto const& [n, p] : obstacles_[i].nodes_) {
      min = {std::min(min.lat_, p.lat_), std::min(min.lng_, p.lng_)};
      max = {std::max(max.lat_, p.lat_), std::max(max.lng_, p.lng_)};
      use[n];
    }
    rtree_insert(obstacle_index, to_array(min).data(), to_array(max).data(),
                 reinterpret_cast<void*>(i));
  }

  for (auto way = way_idx_t{0U}; way != w_.n_ways(); ++way) {
    if (is_outline[to_idx(way)]) {
      continue;
    }
    // A way that cannot be walked on (rails, motorway, ...) and runs along
    // an outline is handled like a barrier: no merging across it.
    auto const p = w_.r_->way_properties_[way];
    auto const walkable = p.is_foot_accessible() || p.is_bike_accessible();
    auto prev = osm_node_idx_t::invalid();
    for (auto const n : w_.way_osm_nodes_[way]) {
      auto const it = use.find(n);
      if (it == end(use)) {
        prev = osm_node_idx_t::invalid();
        continue;
      }
      if (walkable) {
        inc(it->second.other_);
      } else if (prev != osm_node_idx_t::invalid()) {
        barrier_segments_.push_back(segment_hash(prev, n));
      }
      prev = n;
    }
  }
  std::ranges::sort(barrier_segments_);
  std::ranges::sort(blocked_nodes_);
  auto const is_blocked = [&](osm_node_idx_t const n) {
    return std::ranges::binary_search(blocked_nodes_, n);
  };

  auto const is_entry = [&](osm_node_idx_t const n) {
    auto const it = use.find(n);
    return it != end(use) && it->second.other_ != 0U;
  };

  // EVAL-ONLY: OSR_AREA_DUMP=<file> writes the input of every area.
  auto dump_mutex = std::mutex{};
  auto dump = std::unique_ptr<std::ofstream>{};
  if (auto const* f = std::getenv("OSR_AREA_DUMP"); f != nullptr) {
    dump = std::make_unique<std::ofstream>(f);
  }
  auto const dump_area = [&](area const& a, osm_way_idx_t const id) {
    auto out = std::ostringstream{};
    out.precision(12);
    out << "A " << to_idx(id) << ' ' << a.rings_.size() << '\n';
    auto const vertex = [&](area::vertex const& v) {
      out << ' ' << v.id_ << ' ' << v.pos_.lat_ << ' ' << v.pos_.lng_ << ' '
          << (v.entry_ ? 1 : 0) + (v.blocked_ ? 2 : 0);
    };
    for (auto const& r : a.rings_) {
      out << "R " << r.size();
      for (auto const& v : r) {
        vertex(v);
      }
      out << '\n';
    }
    out << "I " << a.interior_.size();
    for (auto const& v : a.interior_) {
      vertex(v);
    }
    out << "\nL " << a.links_.size();
    for (auto const& l : a.links_) {
      out << ' ' << l.from_ << ' ' << l.to_ << ' ' << l.length_;
    }
    out << "\nX " << std::min(a.n_routable_rings_, a.rings_.size()) << ' '
        << a.lines_.size() << '\n';
    for (auto const& l : a.lines_) {
      out << "R " << l.size();
      for (auto const& v : l) {
        vertex(v);
      }
      out << '\n';
    }
    auto const lock = std::scoped_lock{dump_mutex};
    *dump << out.str();
  };

  // Buildings that stand on the area (a node strictly inside, on the level
  // of the area) become holes, barrier segments with an end strictly inside
  // (same level and layer) lines. The outline and the holes of the area
  // decide what is inside: overlapping obstacles are not merged, each blocks
  // on its own.
  auto const add_obstacles = [&](area& a, generated& g) {
    a.n_routable_rings_ = a.rings_.size();
    if (obstacles_.empty() || a.rings_.empty() || a.rings_.front().empty()) {
      return;
    }
    auto const n_rings = a.rings_.size();
    auto const rings = std::span{a.rings_.data(), n_rings};
    auto const is_inside = [&](point const p) {
      return locate_in_rings(rings, p, &area::vertex::pos_) ==
             point_location::kInside;
    };

    auto min = a.rings_.front().front().pos_;
    auto max = min;
    for (auto const& v : a.rings_.front()) {
      min = {std::min(min.lat_, v.pos_.lat_), std::min(min.lng_, v.pos_.lng_)};
      max = {std::max(max.lat_, v.pos_.lat_), std::max(max.lng_, v.pos_.lng_)};
    }
    auto candidates = std::vector<std::size_t>{};
    rtree_search(
        obstacle_index, to_array(min).data(), to_array(max).data(),
        [](double const*, double const*, void const* item, void* udata) {
          static_cast<std::vector<std::size_t>*>(udata)->push_back(
              reinterpret_cast<std::size_t>(item));
          return true;
        },
        &candidates);
    std::ranges::sort(candidates);

    auto const level = ground(g.p_.from_level());
    auto const to_vertex = [&](std::pair<osm_node_idx_t, point> const& x) {
      return area::vertex{.id_ = to_idx(x.first),
                          .pos_ = x.second,
                          .entry_ = is_entry(x.first),
                          .blocked_ = is_blocked(x.first)};
    };
    auto rings_to_add = std::vector<area::ring>{};
    for (auto const i : candidates) {
      auto const& o = obstacles_[i];
      auto const on_level =
          o.has_level_ ? ((o.levels_ >> to_idx(level)) & 1U) != 0U
                       : !o.on_pillars_ && level == level_t{0.F};
      if (!on_level) {
        continue;
      }
      auto const& nodes = o.nodes_;
      if (o.kind_ == obstacle::kind::kBuilding) {
        auto const closed =
            nodes.size() > 3U && nodes.front().first == nodes.back().first;
        if (!closed || std::ranges::none_of(nodes, [&](auto const& x) {
              return is_inside(x.second);
            })) {
          continue;
        }
        auto& ring = rings_to_add.emplace_back();
        for (auto j = 0U; j + 1U < nodes.size(); ++j) {
          ring.push_back(to_vertex(nodes[j]));
        }
        ++g.n_holes_;
      } else if (o.layer_ == g.layer_) {
        // Runs of segments with at least one end inside.
        auto in = std::vector<bool>(nodes.size());
        for (auto j = 0U; j != nodes.size(); ++j) {
          in[j] = is_inside(nodes[j].second);
        }
        auto line = area::ring{};
        for (auto j = 1U; j < nodes.size(); ++j) {
          if (in[j - 1U] || in[j]) {
            if (line.empty()) {
              line.push_back(to_vertex(nodes[j - 1U]));
            }
            line.push_back(to_vertex(nodes[j]));
          } else if (!line.empty()) {
            a.lines_.push_back(std::move(line));
            line = {};
            ++g.n_lines_;
          }
        }
        if (!line.empty()) {
          a.lines_.push_back(std::move(line));
          ++g.n_lines_;
        }
      }
    }
    for (auto& r : rings_to_add) {
      a.rings_.push_back(std::move(r));
    }
  };

  auto const run = [&](area& a, generated& g) {
    add_obstacles(a, g);
    if (dump) {
      dump_area(a, g.id_);
    }
    auto r = compute_area_ways(a, opt_);
    if (r.status_ == area_routing_status::kTooLarge && !a.interior_.empty()) {
      // The limits must not cost an area that works without them.
      a.interior_.clear();
      a.links_.clear();
      r = compute_area_ways(a, opt_);
      g.interior_dropped_ = true;
    }
    g.n_entries_ = r.n_entries_;
    g.n_interior_entries_ = r.n_interior_entries_;
    g.status_ = r.status_;
    g.n_visibility_edges_ = r.n_visibility_edges_;
    g.n_capped_ = r.n_capped_;
    g.junctions_ = std::move(r.junctions_);
    for (auto i = 0U; i != r.ways_.size(); ++i) {
      auto& out = g.ways_.emplace_back();
      for (auto const& v : r.ways_[i]) {
        auto const& x = get_vertex(a, v);
        out.emplace_back(osm_node_idx_t{x.id_}, x.pos_);
      }
      out.front().first = osm_node_idx_t{r.end_ids_[i].first};
      out.back().first = osm_node_idx_t{r.end_ids_[i].second};
    }
  };

  // Collect the areas.
  auto const n_way_areas = way_areas_.size();
  auto const n_areas = n_way_areas + relation_areas_.size();
  auto in = std::vector<input>(n_areas);
  auto areas = std::vector<area>(n_areas);
  utl::parallel_for_run(n_areas, [&](std::size_t const i) {
    auto& a = areas[i];
    auto& x = in[i];
    if (i < n_way_areas) {
      auto const way = way_areas_[i];
      auto const nodes = w_.way_osm_nodes_[way];
      auto const polyline = w_.way_polylines_[way];
      x.id_ = to_area_way_idx(to_idx(w_.way_osm_idx_[way]), false);
      x.p_ = w_.r_->way_properties_[way];
      x.name_ = w_.way_names_[way];
      x.mixed_levels_ = x.p_.from_level_ != x.p_.to_level_;
      x.layer_ = way_area_layers_[i];

      // Stored ways are closed: skip the repeated first node.
      auto& ring = a.rings_.emplace_back();
      for (auto j = 0U; j + 1U < nodes.size(); ++j) {
        ring.push_back({.id_ = to_idx(nodes[j]),
                        .pos_ = polyline[j],
                        .entry_ = is_entry(nodes[j]),
                        .blocked_ = is_blocked(nodes[j])});
      }
    } else {
      auto const& rel = relation_areas_[i - n_way_areas];
      x.id_ = to_area_way_idx(encode_osm_id(rel.id_), true);
      x.p_ = rel.p_;
      x.name_ = string_idx_t::invalid();
      x.layer_ = rel.layer_;
      // Member ways can have their own tags (e.g. steps). A generated way
      // must not connect levels that the outline keeps apart.
      x.mixed_levels_ =
          x.p_.from_level_ != x.p_.to_level_ ||
          std::ranges::any_of(rel.member_ways_, [&](way_idx_t const way) {
            auto const p = w_.r_->way_properties_[way];
            return p.from_level_ != x.p_.from_level_ ||
                   p.to_level_ != x.p_.to_level_;
          });
      for (auto const& r : rel.rings_) {
        auto& ring = a.rings_.emplace_back();
        for (auto const& [node, pos] : r) {
          ring.push_back({.id_ = to_idx(node),
                          .pos_ = pos,
                          .entry_ = is_entry(node),
                          .blocked_ = is_blocked(node)});
        }
      }
    }
    // Entries that are not on the outline.
    auto const idx = static_cast<std::uint32_t>(i);
    for (auto const& e : std::ranges::equal_range(
             inside.entries_, idx, {}, &interior::entry::area_)) {
      a.interior_.push_back(
          {.id_ = to_idx(e.node_), .pos_ = e.pos_, .entry_ = true});
    }
    for (auto const& l : std::ranges::equal_range(inside.links_, idx, {},
                                                  &interior::link::area_)) {
      a.links_.push_back({to_idx(l.from_), to_idx(l.to_), l.length_});
    }
  });

  // Merge adjacent areas: closed ways and multipolygons alike.
  auto classes = std::vector<std::uint64_t>(n_areas);
  for (auto i = 0U; i != n_areas; ++i) {
    classes[i] = in[i].mixed_levels_ ? kNoAreaMerge : to_merge_class(in[i].p_);
  }
  auto merged = merge_areas(
      areas, classes,
      {.is_blocked_ = [&](std::uint64_t const a, std::uint64_t const b) {
        return std::ranges::binary_search(
            barrier_segments_,
            segment_hash(osm_node_idx_t{a}, osm_node_idx_t{b}));
      }});

  // Where two areas only touch in a node, this node is the only way from one
  // into the other: an entry of both. It has to be a routing node: used by
  // two stored outlines (members of a relation are not always stored).
  auto const mark_touching = [&](area& a) {
    for (auto& ring : a.rings_) {
      for (auto& v : ring) {
        if (v.entry_ || !std::ranges::binary_search(merged.touch_, v.id_)) {
          continue;
        }
        auto const it = use.find(osm_node_idx_t{v.id_});
        v.entry_ = it != end(use) && it->second.outline_ >= 2U;
      }
    }
  };

  auto const run_area = [&](std::uint32_t const i, generated& g) {
    g.id_ = in[i].id_;
    g.p_ = in[i].p_;
    g.name_ = in[i].name_;
    g.mixed_levels_ = in[i].mixed_levels_;
    g.layer_ = in[i].layer_;
    if (!g.mixed_levels_) {
      mark_touching(areas[i]);
      run(areas[i], g);
    }
  };

  auto const run_group = [&](merged_area& m, unit& u) {
    // Smallest id, common name of the named areas. The properties that
    // matter are the same for all members (see to_merge_class).
    auto g = generated{};
    g.id_ = in[m.members_.front()].id_;
    g.p_ = in[m.members_.front()].p_;
    g.layer_ = in[m.members_.front()].layer_;
    g.name_ = string_idx_t::invalid();
    auto same_name = true;
    for (auto const i : m.members_) {
      g.id_ = std::min(g.id_, in[i].id_);
      g.p_.in_route_ = g.p_.in_route_ && in[i].p_.in_route_;
      g.p_.motor_vehicle_no_ =
          g.p_.motor_vehicle_no_ && in[i].p_.motor_vehicle_no_;
      if (in[i].name_ != string_idx_t::invalid()) {
        same_name = same_name && (g.name_ == string_idx_t::invalid() ||
                                  g.name_ == in[i].name_);
        g.name_ = in[i].name_;
      }
    }
    if (!same_name) {
      g.name_ = string_idx_t::invalid();
    }

    // Entries inside of the members are inside of the union.
    for (auto const i : m.members_) {
      m.area_.interior_.insert(end(m.area_.interior_),
                               begin(areas[i].interior_),
                               end(areas[i].interior_));
      m.area_.links_.insert(end(m.area_.links_), begin(areas[i].links_),
                            end(areas[i].links_));
    }
    mark_touching(m.area_);
    run(m.area_, g);
    if (g.status_ == area_routing_status::kTooLarge ||
        g.status_ == area_routing_status::kInvalid) {
      u.too_large_ = true;
      for (auto const i : m.members_) {
        run_area(i, u.out_.emplace_back());
      }
      return;
    }

    // The removed segments are still part of the outline ways: no need for
    // a generated way that only follows them.
    u.n_on_removed_ = std::erase_if(g.ways_, [&](ring const& way) {
      auto const is_junction = [&](osm_node_idx_t const n) {
        return std::ranges::binary_search(g.junctions_, to_idx(n));
      };
      if (is_junction(way.front().first) || is_junction(way.back().first)) {
        return false;  // the outline does not reach junctions
      }
      for (auto i = 1U; i < way.size(); ++i) {
        auto const a = to_idx(way[i - 1U].first);
        auto const b = to_idx(way[i].first);
        if (!std::ranges::binary_search(
                m.removed_, std::pair{std::min(a, b), std::max(a, b)})) {
          return false;
        }
      }
      return true;
    });
    u.out_.emplace_back(std::move(g));
  };

  // Groups first, then the areas that are not part of a group.
  auto single = std::vector<std::uint32_t>{};
  for (auto i = 0U; i != n_areas; ++i) {
    if (merged.group_[i] == area_merge_result::kNoGroup) {
      single.push_back(i);
    }
  }
  auto const n_groups = merged.merged_.size();
  auto units = std::vector<unit>(n_groups + single.size());
  utl::parallel_for_run(units.size(), [&](std::size_t const i) {
    if (i < n_groups) {
      run_group(merged.merged_[i], units[i]);
    } else {
      run_area(single[i - n_groups], units[i].out_.emplace_back());
    }
  });

  s.n_areas_ = n_areas;
  s.n_removed_segments_ = merged.n_removed_;
  s.n_incompatible_segments_ = merged.n_incompatible_;
  s.n_blocked_segments_ = merged.n_blocked_;
  s.n_conflict_segments_ = merged.n_conflicts_;
  s.n_failed_groups_ = merged.n_failed_groups_;
  s.n_touch_nodes_ = merged.touch_.size();

  auto out = std::vector<generated const*>{};
  for (auto i = 0U; i != units.size(); ++i) {
    auto const& u = units[i];
    for (auto const& g : u.out_) {
      out.push_back(&g);
    }
    s.n_on_removed_ += u.n_on_removed_;
    if (u.too_large_) {
      ++s.n_too_large_groups_;
    } else if (i < n_groups) {
      auto const size = merged.merged_[i].members_.size();
      ++s.n_merged_groups_;
      s.n_merged_areas_ += size;
      ++s.group_sizes_[size <= 4U   ? size - 2U
                       : size <= 8U ? 3U
                                    : (size <= 16U ? 4U : 5U)];
    }
  }

  // Append in id order to keep way_osm_idx_ sorted.
  std::ranges::stable_sort(out, {}, [](auto const* g) { return g->id_; });

  for (auto const* gp : out) {
    auto const& g = *gp;
    if (g.mixed_levels_) {
      ++s.n_mixed_levels_;
      continue;
    }
    s.n_interior_entries_ += g.n_interior_entries_;
    s.n_building_holes_ += g.n_holes_;
    s.n_barrier_lines_ += g.n_lines_;
    s.n_with_interior_ += g.n_interior_entries_ != 0U ? 1U : 0U;
    s.n_interior_dropped_ += g.interior_dropped_ ? 1U : 0U;
    switch (g.status_) {
      case area_routing_status::kNoEntries: ++s.n_no_entries_; continue;
      case area_routing_status::kTooLarge: ++s.n_too_large_; continue;
      case area_routing_status::kInvalid: ++s.n_invalid_; continue;
      case area_routing_status::kOk: break;
    }
    ++s.n_ok_;
    s.n_ok_by_interior_ +=
        g.n_entries_ - g.n_interior_entries_ < 2U ? 1U : 0U;
    s.n_visibility_edges_ += g.n_visibility_edges_;
    s.n_capped_ += g.n_capped_;

    auto const p = to_area_way_properties(g.p_);
    auto const name = g.name_;

    // Junctions: new routing nodes, one per vertex id of this area.
    auto junction_ids = hash_map<std::uint64_t, osm_node_idx_t>{};
    auto const to_node = [&](std::pair<osm_node_idx_t, point> const& v) {
      if (!std::ranges::binary_search(g.junctions_, to_idx(v.first))) {
        return v.first;
      }
      auto const [it, added] = junction_ids.emplace(
          to_idx(v.first), osm_node_idx_t{kAreaNodeBit | junctions_.size()});
      if (added) {
        junctions_.push_back({it->second, v.second, p});
        w_.node_way_counter_.increment(to_idx(it->second));
        w_.node_way_counter_.increment(to_idx(it->second));
        ++s.n_junctions_;
      }
      return it->second;
    };

    for (auto const& way : g.ways_) {
      // ways::connect_ways drops ways at nodes with more than kMaxWaysPerNode
      // ways - which could hit OSM ways as well. Stay below. Junctions have
      // at most max_ways_per_node_ ways.
      auto const is_junction = [&](osm_node_idx_t const n) {
        return std::ranges::binary_search(g.junctions_, to_idx(n));
      };
      auto const is_full = [&](osm_node_idx_t const n) {
        if (is_junction(n)) {
          return false;
        }
        auto const& x = use.at(n);
        return x.outline_ + x.other_ + x.generated_ + 2 >= kMaxWaysPerNode;
      };
      if (is_full(way.front().first) || is_full(way.back().first)) {
        ++s.n_over_node_budget_;
        continue;
      }
      for (auto const n : {way.front().first, way.back().first}) {
        if (!is_junction(n)) {
          inc(use.at(n).generated_);
        }
      }

      auto const way_idx = way_idx_t{w_.way_osm_idx_.size()};
      w_.way_osm_idx_.push_back(g.id_);
      w_.r_->way_properties_.emplace_back(p);
      w_.way_polylines_.emplace_back(
          way | std::views::transform([](auto const& x) { return x.second; }));
      // Only the ends are connected to the routing graph. A bend that is a
      // routing node (e.g. the first node of a closed hole outline) would
      // otherwise become a junction - outside of the per node way budget.
      auto nodes = std::vector<osm_node_idx_t>(way.size(),
                                               osm_node_idx_t::invalid());
      nodes.front() = to_node(way.front());
      nodes.back() = to_node(way.back());
      w_.way_osm_nodes_.emplace_back(nodes);
      w_.way_names_.emplace_back(name);
      w_.way_has_conditional_access_no_.resize(to_idx(way_idx) + 1U);

      // Entries on the outline are routing nodes already (outline + other
      // way). The end of a way inside the area is not: it was counted once
      // and becomes a routing node now, which only extends its way by one
      // edge. Not done for entries without generated way.
      for (auto const n : {nodes.front(), nodes.back()}) {
        if (!w_.node_way_counter_.is_multi(to_idx(n))) {
          ++s.n_new_routing_nodes_;
          do {
            w_.node_way_counter_.increment(to_idx(n));
          } while (!w_.node_way_counter_.is_multi(to_idx(n)));
        }
      }

      ++s.n_ways_;
      s.n_way_nodes_ += way.size();
    }
  }
  rtree_free(obstacle_index);
  return s;
}

void area_ways::set_junction_nodes() {
  for (auto const& j : junctions_) {
    auto const n = w_.find_node_idx(j.id_);
    if (!n.has_value()) {
      continue;  // all of its ways were dropped
    }
    auto p = node_properties{};
    std::memset(&p, 0, sizeof(node_properties));
    p.is_foot_accessible_ = j.p_.is_foot_accessible_;
    p.is_bike_accessible_ = j.p_.is_bike_accessible_;
    p.from_level_ = j.p_.from_level_;
    p.to_level_ = j.p_.to_level_;
    w_.r_->node_properties_[*n] = p;
    w_.r_->node_positions_[*n] = j.pos_;
  }
}

}  // namespace osr
