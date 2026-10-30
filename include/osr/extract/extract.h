#pragma once

#include <cstdint>
#include <filesystem>

#include "osr/extract/walkable_areas.h"

namespace osr {

// Ways generated across walkable areas (pedestrian squares, platforms, ...),
// see area_ways.h. Off by default: the extract is the same as without areas.
struct area_routing_config {
  enum class approach : std::uint8_t {
    kNone,
    // Shortest paths between the access points of an area, only as many as
    // needed to stay within max_stretch_ (greedy spanner).
    kSpanner,
    // All lines of sight between the corners and access points (the
    // visibility graph). Baseline for comparisons.
    kVisibilityGraph
  };

  approach approach_{approach::kNone};

  // kSpanner: a pair of access points does not get its own path if the
  // generated ways connect it with a path at most this factor longer.
  double max_stretch_{1.02};

  // kSpanner: upper bound for the generated ways at one node.
  std::uint32_t max_ways_per_node_{24U};
};

// `on_walkable_areas`: optional, see walkable_areas_hook. Without it, walkable
// areas are not collected and the graph is the same as before.
void extract(bool with_platforms,
             std::filesystem::path const& in,
             std::filesystem::path const& out,
             std::filesystem::path const& elevation_dir,
             walkable_areas_hook const& on_walkable_areas = {},
             area_routing_config const& = {});

}  // namespace osr
