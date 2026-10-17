#pragma once

#include <cstdint>
#include <limits>
#include <vector>

#include "osr/point.h"

namespace osr {

// A walkable polygon: one outer ring plus its holes.
//
// Rings are open (the first vertex is not repeated at the end) and may be given
// in any winding order. Vertices that represent the same OSM node (a ring
// touching itself or another ring) carry the same id_.
struct area {
  struct vertex {
    std::uint64_t id_;
    point pos_;
    // Entry: the vertex connects the area to the rest of the routing graph.
    bool entry_;
    // The outline cannot be walked through this vertex (e.g. a node with
    // access=private): with ring_routable_, the outline does not connect the
    // entries on both sides of it. Generated ways can still bend here.
    bool blocked_{false};
  };

  using ring = std::vector<vertex>;

  // rings_[0] is the outer ring, all others are holes. A hole that consists
  // of a single entry is a point inside the area that can be walked to and
  // from in any direction (not an obstacle).
  std::vector<ring> rings_;

  // Entries that are not on the outline (e.g. the end of a footway that
  // leads into the area): free-standing points that can be reached from all
  // directions and do not block anything. entry_ is not used. Vertices that
  // are not strictly inside the area (on the outline, in a hole) or that
  // have the id of an outline vertex are ignored.
  std::vector<vertex> interior_;

  // A connection between two entries (ids) that is part of the routing graph
  // already, like the outline if ring_routable_ is set: the footway between
  // an outline node and its end inside the area. Pairs that are connected
  // well enough this way do not get a generated way. Links between vertices
  // that are not entries are ignored.
  struct link {
    std::uint64_t from_, to_;
    double length_;  // in meters
  };
  std::vector<link> links_;
};

struct area_routing_options {
  // Areas with more candidate corners (reflex corners + entries) or more
  // entries are skipped.
  std::uint32_t max_candidates_{4096U};
  std::uint32_t max_entries_{512U};

  // A pair of entries does not get its own way if the other generated ways
  // connect it with a path that is at most this factor longer than the
  // shortest path. 1.0 = always optimal.
  double max_stretch_{1.02};

  // Upper bound for the number of generated ways that end at one entry or
  // junction (ways::connect_ways drops ways at nodes with too many ways).
  std::uint32_t max_ways_per_node_{24U};

  // Instead of the selected shortest paths, generate the whole search graph:
  // all tangent lines of sight plus the outline between the corners (the
  // visibility graph). For comparison only.
  bool all_visibility_edges_{false};

  // The outline of the area is already part of the routing graph (except at
  // blocked vertices): paths that only follow the outline are not generated.
  bool ring_routable_{true};
};

enum class area_routing_status : std::uint8_t {
  kOk,
  kNoEntries,  // less than two entries, nothing to connect
  kTooLarge,  // extent or number of candidates over the limit
  kInvalid  // degenerate geometry
};

struct area_routing_result {
  // A vertex of the input: rings_[ring_][idx_], or interior_[idx_] if
  // ring_ == kInterior.
  struct vertex_ref {
    static constexpr auto kInterior =
        std::numeric_limits<std::uint32_t>::max();
    std::uint32_t ring_, idx_;
  };

  area_routing_status status_{area_routing_status::kOk};

  // Every way ends at entries or junctions: corners where more than two ways
  // meet, which are not entries. Only the ends have to be connected to the
  // routing graph, the vertices in between are geometry. A junction is a
  // new routing node (the same one for all ways that end at the same vertex
  // id), not the outline node: the outline does not pass through it.
  std::vector<std::vector<vertex_ref>> ways_;

  // Vertex ids of the junctions, sorted.
  std::vector<std::uint64_t> junctions_;

  // Statistics: line-of-sight edges that survived all filters, entry pairs
  // that did not get a path because of max_ways_per_node_.
  std::uint64_t n_visibility_edges_{0U};
  std::uint64_t n_capped_{0U};

  // Statistics: number of entries (set for kNoEntries and kTooLarge as
  // well), and how many of them are not on the outline.
  std::uint32_t n_entries_{0U};
  std::uint32_t n_interior_entries_{0U};
};

inline area::vertex const& get_vertex(
    area const& a, area_routing_result::vertex_ref const v) {
  return v.ring_ == area_routing_result::vertex_ref::kInterior
             ? a.interior_[v.idx_]
             : a.rings_[v.ring_][v.idx_];
}

// Computes the ways required to cross the area on (nearly) shortest paths
// between all pairs of entries.
//
// The shortest path between two points inside a polygon only bends at reflex
// corners, and only along lines that are tangent to the polygon there. So the
// search graph consists of the outline plus all tangent lines of sight between
// reflex corners and entries. Entries inside the area are points of this
// graph without outline. All geometric predicates are exact (integer).
// Out of the shortest paths between the entries, a greedy spanner selects the
// ones that are needed to stay within max_stretch_. Selected paths share
// their common segments (junctions where they meet).
area_routing_result compute_area_ways(area const&,
                                      area_routing_options const& = {});

}  // namespace osr
