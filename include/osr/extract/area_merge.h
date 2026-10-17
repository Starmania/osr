#pragma once

#include <cstdint>
#include <functional>
#include <limits>
#include <utility>
#include <vector>

#include "osr/extract/area_routing.h"

namespace osr {

// Areas with this class are never merged (and do not make other areas touch).
constexpr auto const kNoAreaMerge = std::numeric_limits<std::uint64_t>::max();

struct area_merge_options {
  // Called with the vertex ids of an outline segment shared by two areas:
  // something that must not be crossed (e.g. a wall) is mapped on it. The
  // areas are not merged across this segment.
  std::function<bool(std::uint64_t, std::uint64_t)> is_blocked_;
};

// The union of two or more areas.
struct merged_area {
  // rings_[0] is the outer ring, followed by the holes of the union, followed
  // by one single vertex ring for every entry that is not part of the outline
  // of the union anymore (it was on a removed segment). Outer ring
  // counter-clockwise, holes clockwise.
  area area_;

  // Indices of the merged input areas, ascending.
  std::vector<std::uint32_t> members_;

  // The removed outline segments (smaller vertex id first), sorted.
  std::vector<std::pair<std::uint64_t, std::uint64_t>> removed_;
};

struct area_merge_result {
  static constexpr auto const kNoGroup =
      std::numeric_limits<std::uint32_t>::max();

  // Ordered by the smallest member.
  std::vector<merged_area> merged_;

  // For each input area: index into merged_, kNoGroup = not merged.
  std::vector<std::uint32_t> group_;

  // Vertices where outlines of different areas meet without sharing a segment
  // (after merging). Crossing from one area to the other is only possible
  // there. Sorted vertex ids.
  std::vector<std::uint64_t> touch_;

  // Statistics: shared segments that were removed / were not used to merge
  // because the classes differ / is_blocked_ returned true / the areas are on
  // the same side of the segment or more than two outlines use it.
  std::uint64_t n_removed_{0U};
  std::uint64_t n_incompatible_{0U};
  std::uint64_t n_blocked_{0U};
  std::uint64_t n_conflicts_{0U};

  // Groups (and their areas) that are not merged because the remaining
  // segments do not form one outer ring with holes (overlapping or otherwise
  // invalid input).
  std::uint64_t n_failed_groups_{0U};
  std::uint64_t n_failed_areas_{0U};
};

// Merges areas that share outline segments.
//
// Two areas are merged if they have the same class and a segment (two
// consecutive vertices with the same ids) in common that they are on
// opposite sides of. No general polygon union is required: the shared
// segments are removed and the remaining segments are linked to rings again.
// Where more than two segments remain at a vertex, the ring continues along
// the walkable wedge - so rings may touch but never cross.
//
// Not merged: areas that only touch in a vertex (see touch_), areas that
// overlap, and areas where one outline has a vertex on the shared part that
// the other does not have (no common segment - the two ends are in touch_).
//
// A shared segment that is not removed between two areas of the same group
// (is_blocked_, or used by more than two outlines) stays in the union as a
// wall with both sides.
area_merge_result merge_areas(std::vector<area> const&,
                              std::vector<std::uint64_t> const& classes,
                              area_merge_options const& = {});

}  // namespace osr
