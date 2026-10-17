#include "osr/extract/area_merge.h"

#include <algorithm>
#include <numeric>
#include <tuple>

namespace osr {

namespace {

// Cross products of coordinate differences do not fit into 64 bit.
using wide_t = __int128;

constexpr auto kInvalid = std::numeric_limits<std::uint32_t>::max();

struct vec2 {
  friend vec2 operator-(vec2 const a, vec2 const b) {
    return {a.x_ - b.x_, a.y_ - b.y_};
  }
  std::int64_t x_, y_;
};

vec2 to_vec(point const p) { return {p.lat_, p.lng_}; }

wide_t cross(vec2 const a, vec2 const b) {
  return static_cast<wide_t>(a.x_) * b.y_ - static_cast<wide_t>(a.y_) * b.x_;
}

wide_t dot(vec2 const a, vec2 const b) {
  return static_cast<wide_t>(a.x_) * b.x_ + static_cast<wide_t>(a.y_) * b.y_;
}

// > 0: counter-clockwise, < 0: clockwise, 0: degenerate.
template <typename Ring>
int winding(Ring const& ring) {
  auto doubled_area = wide_t{0};
  auto const origin = to_vec(ring[0].pos_);
  for (auto i = 1U; i + 1U < ring.size(); ++i) {
    doubled_area += cross(to_vec(ring[i].pos_) - origin,
                          to_vec(ring[i + 1U].pos_) - origin);
  }
  return doubled_area > 0 ? 1 : (doubled_area < 0 ? -1 : 0);
}

// Where the direction d comes when rotating clockwise, starting after r:
// 0 = within the first half turn, 1 = opposite to r, 2 = within the second
// half turn, 3 = same direction as r, -1 = no direction.
int clockwise_rank(vec2 const r, vec2 const d) {
  auto const c = cross(r, d);
  if (c != 0) {
    return c < 0 ? 0 : 2;
  }
  auto const t = dot(r, d);
  return t < 0 ? 1 : (t > 0 ? 3 : -1);
}

// A directed outline segment with the walkable side on the left:
// rings_[ring_][from_] -> rings_[ring_][to_] of an input area.
struct segment {
  enum class state : std::uint8_t {
    kKeep,
    kShared,  // to be removed if the group can be merged
    kRemoved
  };

  std::uint32_t area_, ring_, from_, to_;
  state state_{state::kKeep};
};

struct merger {
  merger(std::vector<area> const& areas,
         std::vector<std::uint64_t> const& classes,
         area_merge_options const& opt)
      : areas_{areas}, classes_{classes}, opt_{opt} {}

  area::vertex const& from(segment const& s) const {
    return areas_[s.area_].rings_[s.ring_][s.from_];
  }

  area::vertex const& to(segment const& s) const {
    return areas_[s.area_].rings_[s.ring_][s.to_];
  }

  std::pair<std::uint64_t, std::uint64_t> key(segment const& s) const {
    return {std::min(from(s).id_, to(s).id_),
            std::max(from(s).id_, to(s).id_)};
  }

  std::uint32_t find(std::uint32_t x) {
    while (parent_[x] != x) {
      parent_[x] = parent_[parent_[x]];
      x = parent_[x];
    }
    return x;
  }

  // The smaller index becomes the root: deterministic.
  void unite(std::uint32_t const a, std::uint32_t const b) {
    auto const x = find(a);
    auto const y = find(b);
    parent_[std::max(x, y)] = std::min(x, y);
  }

  void collect_segments() {
    area_segments_.resize(areas_.size());
    for (auto a = 0U; a != areas_.size(); ++a) {
      auto const begin = static_cast<std::uint32_t>(segments_.size());
      area_segments_[a] = {begin, begin};
      if (classes_[a] == kNoAreaMerge) {
        continue;
      }
      auto const& rings = areas_[a].rings_;
      for (auto r = 0U; r != rings.size(); ++r) {
        auto const n = static_cast<std::uint32_t>(rings[r].size());
        auto const w = n < 3U ? 0 : winding(rings[r]);
        if (w == 0) {
          if (r == 0U) {
            break;  // no valid outer ring
          }
          continue;
        }
        // Outer ring counter-clockwise, holes clockwise.
        auto const forward = (w > 0) == (r == 0U);
        for (auto i = 0U; i != n; ++i) {
          auto const j = forward ? (i + 1U) % n : (i + n - 1U) % n;
          if (rings[r][i].id_ != rings[r][j].id_) {
            segments_.push_back({.area_ = a, .ring_ = r, .from_ = i, .to_ = j});
          }
        }
      }
      area_segments_[a].second = static_cast<std::uint32_t>(segments_.size());
    }

    // Both directions of a segment next to each other.
    by_key_.resize(segments_.size());
    std::iota(begin(by_key_), end(by_key_), 0U);
    std::ranges::sort(by_key_, [&](auto const x, auto const y) {
      return std::pair{key(segments_[x]), x} < std::pair{key(segments_[y]), y};
    });
  }

  // Calls fn(first, last) for every range of by_key_ with more than one
  // segment between the same two vertices.
  template <typename Fn>
  void for_each_shared(Fn&& fn) const {
    for (auto b = 0U; b != by_key_.size();) {
      auto e = b + 1U;
      while (e != by_key_.size() &&
             key(segments_[by_key_[e]]) == key(segments_[by_key_[b]])) {
        ++e;
      }
      if (e - b > 1U) {
        fn(b, e);
      }
      b = e;
    }
  }

  void find_groups() {
    parent_.resize(areas_.size());
    std::iota(begin(parent_), end(parent_), 0U);
    for_each_shared([&](std::uint32_t const b, std::uint32_t const e) {
      auto& s1 = segments_[by_key_[b]];
      auto& s2 = segments_[by_key_[b + 1U]];
      if (e - b == 2U && s1.area_ != s2.area_ &&
          from(s1).id_ == to(s2).id_) {
        // Two areas on opposite sides.
        if (classes_[s1.area_] != classes_[s2.area_]) {
          ++result_.n_incompatible_;
        } else if (opt_.is_blocked_ &&
                   opt_.is_blocked_(from(s1).id_, to(s1).id_)) {
          ++result_.n_blocked_;
        } else {
          s1.state_ = s2.state_ = segment::state::kShared;
          unite(s1.area_, s2.area_);
        }
      } else if (std::any_of(begin(by_key_) + b, begin(by_key_) + e,
                             [&](auto const x) {
                               return segments_[x].area_ != s1.area_;
                             })) {
        ++result_.n_conflicts_;
      }
    });
  }

  // Links the segments of the group that are not shared to rings.
  bool link(std::vector<std::uint32_t> const& members, merged_area& m) {
    // remaining_: in input order, by_from_: indices into remaining_.
    remaining_.clear();
    for (auto const a : members) {
      for (auto s = area_segments_[a].first; s != area_segments_[a].second;
           ++s) {
        if (segments_[s].state_ == segment::state::kKeep) {
          remaining_.push_back(s);
        }
      }
    }
    auto const seg = [&](std::uint32_t const i) -> segment const& {
      return segments_[remaining_[i]];
    };
    by_from_.resize(remaining_.size());
    std::iota(begin(by_from_), end(by_from_), 0U);
    std::ranges::sort(by_from_, [&](auto const x, auto const y) {
      return std::tuple{from(seg(x)).id_, to(seg(x)).id_, x} <
             std::tuple{from(seg(y)).id_, to(seg(y)).id_, y};
    });

    // The segment that follows s on the outline of the union: the first one
    // that leaves the end of s when rotating clockwise, starting at s. The
    // wedge in between is walkable (it is left of both).
    auto const next = [&](segment const& s) {
      auto const node = to(s).id_;
      auto const [b, e] = std::ranges::equal_range(
          by_from_, node, {}, [&](auto const x) { return from(seg(x)).id_; });
      if (e - b == 1) {
        return *b;
      }
      auto const pos = to_vec(to(s).pos_);
      auto const back = to_vec(from(s).pos_) - pos;
      auto best = kInvalid;
      auto best_rank = 4;
      auto best_dir = vec2{};
      for (auto it = b; it != e; ++it) {
        auto const dir = to_vec(to(seg(*it)).pos_) - pos;
        auto const rank = clockwise_rank(back, dir);
        if (rank == -1) {
          return kInvalid;
        }
        auto const c = rank == best_rank ? cross(dir, best_dir) : wide_t{0};
        if (rank == best_rank && c == 0) {
          return kInvalid;  // two segments in the same direction: overlap
        }
        if (rank < best_rank || (rank == best_rank && c < 0)) {
          best = *it;
          best_rank = rank;
          best_dir = dir;
        }
      }
      return best;
    };

    auto visited = std::vector<bool>(remaining_.size());
    auto n_outer = 0U;
    m.area_.rings_.emplace_back();  // outer ring
    for (auto start = 0U; start != remaining_.size(); ++start) {
      if (visited[start]) {
        continue;
      }
      auto ring = area::ring{};
      auto i = start;
      do {
        visited[i] = true;
        ring.push_back(from(seg(i)));
        i = next(seg(i));
        if (i == kInvalid || (visited[i] && i != start)) {
          return false;  // not a ring
        }
      } while (i != start);

      auto const w = ring.size() < 3U ? 0 : winding(ring);
      if (w == 0) {
        return false;  // e.g. a wall inside of the area
      } else if (w > 0) {
        ++n_outer;
        m.area_.rings_[0] = std::move(ring);
      } else {
        m.area_.rings_.emplace_back(std::move(ring));
      }
    }
    return n_outer == 1U;
  }

  void merge(std::vector<std::uint32_t> const& members) {
    auto m = merged_area{};
    auto const ok = link(members, m);

    // Entries on removed segments that are not on the outline anymore.
    auto on_outline = std::vector<std::uint64_t>{};
    auto inside = std::vector<area::vertex>{};
    if (ok) {
      for (auto const& ring : m.area_.rings_) {
        for (auto const& v : ring) {
          on_outline.push_back(v.id_);
        }
      }
      std::ranges::sort(on_outline);
    }

    for (auto const a : members) {
      for (auto i = area_segments_[a].first; i != area_segments_[a].second;
           ++i) {
        auto& s = segments_[i];
        if (s.state_ != segment::state::kShared) {
          continue;
        }
        s.state_ = ok ? segment::state::kRemoved : segment::state::kKeep;
        if (!ok) {
          continue;
        }
        if (from(s).id_ < to(s).id_) {
          m.removed_.emplace_back(from(s).id_, to(s).id_);
        }
        if (from(s).entry_ &&
            !std::ranges::binary_search(on_outline, from(s).id_)) {
          inside.push_back(from(s));
        }
      }
    }

    if (!ok) {
      ++result_.n_failed_groups_;
      result_.n_failed_areas_ += members.size();
      return;
    }

    std::ranges::sort(inside, {}, &area::vertex::id_);
    for (auto i = 0U; i != inside.size(); ++i) {
      if (i == 0U || inside[i].id_ != inside[i - 1U].id_) {
        m.area_.rings_.push_back({inside[i]});
      }
    }
    std::ranges::sort(m.removed_);
    result_.n_removed_ += m.removed_.size();
    m.members_ = members;
    for (auto const a : members) {
      result_.group_[a] = static_cast<std::uint32_t>(result_.merged_.size());
    }
    result_.merged_.emplace_back(std::move(m));
  }

  void merge_groups() {
    result_.group_.assign(areas_.size(), area_merge_result::kNoGroup);

    // Members in ascending order, groups ordered by their first member (the
    // root, see unite).
    auto members = std::vector<std::vector<std::uint32_t>>(areas_.size());
    for (auto a = 0U; a != areas_.size(); ++a) {
      members[find(a)].push_back(a);
    }
    for (auto const& group : members) {
      if (group.size() > 1U) {
        merge(group);
      }
    }
  }

  void find_touching() {
    // Vertices on a segment that is still used by more than one area.
    auto shared = std::vector<std::uint64_t>{};
    for_each_shared([&](std::uint32_t const b, std::uint32_t const e) {
      auto first = kInvalid;
      for (auto i = b; i != e; ++i) {
        auto const& s = segments_[by_key_[i]];
        if (s.state_ == segment::state::kRemoved) {
          continue;
        }
        if (first == kInvalid) {
          first = s.area_;
        } else if (s.area_ != first) {
          shared.push_back(from(s).id_);
          shared.push_back(to(s).id_);
          return;
        }
      }
    });
    std::ranges::sort(shared);

    // Vertices on the remaining outlines of more than one area.
    auto uses = std::vector<std::pair<std::uint64_t, std::uint32_t>>{};
    for (auto const& s : segments_) {
      if (s.state_ != segment::state::kRemoved) {
        uses.emplace_back(from(s).id_, s.area_);
      }
    }
    std::ranges::sort(uses);
    for (auto i = 1U; i < uses.size(); ++i) {
      auto const id = uses[i].first;
      if (uses[i - 1U].first == id && uses[i - 1U].second != uses[i].second &&
          (result_.touch_.empty() || result_.touch_.back() != id) &&
          !std::ranges::binary_search(shared, id)) {
        result_.touch_.push_back(id);
      }
    }
  }

  std::vector<area> const& areas_;
  std::vector<std::uint64_t> const& classes_;
  area_merge_options const& opt_;

  std::vector<segment> segments_;
  std::vector<std::pair<std::uint32_t, std::uint32_t>> area_segments_;
  std::vector<std::uint32_t> by_key_;
  std::vector<std::uint32_t> parent_;
  std::vector<std::uint32_t> remaining_, by_from_;

  area_merge_result result_;
};

}  // namespace

area_merge_result merge_areas(std::vector<area> const& areas,
                              std::vector<std::uint64_t> const& classes,
                              area_merge_options const& opt) {
  auto m = merger{areas, classes, opt};
  m.collect_segments();
  m.find_groups();
  m.merge_groups();
  m.find_touching();
  return std::move(m.result_);
}

}  // namespace osr
