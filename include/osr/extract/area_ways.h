#pragma once

#include <cstdint>
#include <array>
#include <mutex>
#include <span>
#include <utility>
#include <vector>

#include "osr/extract/area_routing.h"
#include "osr/point.h"
#include "osr/types.h"
#include "osr/ways.h"

namespace osr {

// Collects walkable areas during the extract and adds the ways that cross
// them (see compute_area_ways) to the ways, as if they were mapped in OSM.
//
// Has to run after all OSM ways are stored (entries of an area are the outline
// nodes that are also used by a way which is not an area outline, and the
// ends and routing nodes of ways inside the area) and before
// ways::connect_ways().
struct area_ways {
  struct stats {
    std::uint64_t n_areas_{0U};
    std::uint64_t n_no_entries_{0U};
    std::uint64_t n_too_large_{0U};
    std::uint64_t n_invalid_{0U};
    std::uint64_t n_mixed_levels_{0U};
    std::uint64_t n_ok_{0U};
    std::uint64_t n_ways_{0U};
    std::uint64_t n_way_nodes_{0U};
    std::uint64_t n_visibility_edges_{0U};
    std::uint64_t n_capped_{0U};
    std::uint64_t n_over_node_budget_{0U};
    std::uint64_t n_junctions_{0U};
    // Entries that are not on the outline: how many, in how many areas, how
    // many areas have two entries only because of them, how many areas had
    // to do without them (limits), and how many of them were not routing
    // nodes before.
    std::uint64_t n_interior_entries_{0U};
    std::uint64_t n_with_interior_{0U};
    std::uint64_t n_ok_by_interior_{0U};
    std::uint64_t n_interior_dropped_{0U};
    std::uint64_t n_new_routing_nodes_{0U};
    // Places where a walkable way crosses an area outline without a common
    // node (not connected to the area there).
    std::uint64_t n_crossings_{0U};
    // Obstacles: buildings standing on areas (holes), barrier lines.
    std::uint64_t n_building_holes_{0U};
    std::uint64_t n_barrier_lines_{0U};

    // Merging of adjacent areas (see merge_areas). Merged areas count as one
    // in n_no_entries_, n_too_large_, n_invalid_ and n_ok_.
    std::uint64_t n_merged_groups_{0U};
    std::uint64_t n_merged_areas_{0U};
    // Number of groups with 2, 3, 4, 5-8, 9-16, more areas.
    std::array<std::uint64_t, 6U> group_sizes_{};
    std::uint64_t n_removed_segments_{0U};
    // Shared segments not used for merging: different properties / levels,
    // barrier or not walkable way on the segment, overlapping areas.
    std::uint64_t n_incompatible_segments_{0U};
    std::uint64_t n_blocked_segments_{0U};
    std::uint64_t n_conflict_segments_{0U};
    // Groups processed as separate areas: the union is not a polygon / is
    // over the size limits.
    std::uint64_t n_failed_groups_{0U};
    std::uint64_t n_too_large_groups_{0U};
    // Nodes where areas only touch: entries of both areas.
    std::uint64_t n_touch_nodes_{0U};
    // Generated ways dropped because they only follow removed segments,
    // which are still part of the outline ways.
    std::uint64_t n_on_removed_{0U};
  };

  using ring = std::vector<std::pair<osm_node_idx_t, point>>;

  // Which nodes of walkable ways inside an area are entries, in addition to
  // the entries on the outline (see area_ways_interior.cc).
  enum class interior_entries : std::uint8_t {
    kNone,
    // Ends of ways that are not connected to any other way.
    kDeadEnds,
    // All ends of ways and all nodes that are routing nodes already.
    kRoutingNodes,
    // Every node of a way that lies inside the area.
    kAllNodes
  };

  // A building or a barrier that may stand on an area (see add_obstacles).
  struct obstacle {
    enum class kind : std::uint8_t { kBuilding, kBarrier };
    osm_way_idx_t id_;
    kind kind_;
    // The levels it is on (level tag), the ground level without.
    level_bits_t levels_;
    bool has_level_;
    // building:min_level > 0 or min_height > 2.5m: does not block the ground
    // below, only its own levels.
    bool on_pillars_;
    std::int8_t layer_;
    std::vector<std::pair<osm_node_idx_t, point>> nodes_;
  };

  explicit area_ways(
      ways& w,
      area_routing_options const& opt = {},
      interior_entries const interior = interior_entries::kAllNodes);
  ~area_ways();
  area_ways(area_ways const&) = delete;
  area_ways& operator=(area_ways const&) = delete;

  // A closed way that is an area. Its outline is the stored way.
  void add_way(way_idx_t, std::int8_t layer = 0);

  // One outer ring of a multipolygon relation with its holes. Rings are open
  // (first node not repeated). `member_ways` are the stored ways the outline
  // is made of. The generated ways get the properties of the relation.
  void add_relation(std::int64_t relation_id,
                    std::vector<ring> rings,
                    std::vector<way_idx_t> member_ways,
                    way_properties,
                    bool off_ground = false,
                    std::int8_t layer = 0);

  // After all areas are added: is the box (min, max) close to an area?
  // Obstacles elsewhere are not needed.
  void index_areas();
  bool is_near_area(point min, point max) const;

  // Buildings and barriers close to areas. Buildings with a node strictly
  // inside an area (on its level) are holes, barrier segments with an end
  // strictly inside (same level and layer) are lines that cannot be crossed.
  // Can be called concurrently.
  void add_obstacles(std::vector<obstacle>&&);
  std::size_t n_obstacles() const { return obstacles_.size(); }

  // The way is not on the ground (tunnel, bridge, layer): what is inside an
  // area on the map is not on the area then, and the other way round.
  void set_off_ground(way_idx_t);

  // Steps, ramps and elevators connect levels: their nodes inside an area
  // are entries whatever their level (see area_ways_interior.cc).
  void set_connector(way_idx_t);

  // A node with highway=elevator: serves all levels.
  void add_elevator_node(osm_node_idx_t);

  // A way that is not stored but must not be crossed (wall, fence, ...):
  // areas are not merged across it. Can be called concurrently.
  void add_barrier(std::span<osm_node_idx_t const>);

  // A node that cannot be walked through (e.g. access=private): the outline
  // does not connect the entries on both sides of it.
  void add_blocked_node(osm_node_idx_t);

  stats build();

  // Positions and properties of the junctions (see compute_area_ways): call
  // after the node properties of the OSM nodes are set.
  void set_junction_nodes();

private:
  struct relation_area {
    std::int64_t id_;
    std::vector<ring> rings_;
    std::vector<way_idx_t> member_ways_;
    way_properties p_;
    bool off_ground_;
    std::int8_t layer_;
  };

  // Entries that are not on the outline, see area_ways_interior.cc. area_ is
  // the index of the area: way areas first, then the relation areas.
  struct interior {
    struct entry {
      std::uint32_t area_;
      osm_node_idx_t node_;
      point pos_;
    };
    struct link {
      std::uint32_t area_;
      osm_node_idx_t from_, to_;
      double length_;
    };
    std::vector<entry> entries_;  // sorted by area_ and node_, unique
    std::vector<link> links_;  // sorted by area_
    // Walkable ways crossing an area outline without a common node.
    std::uint64_t n_crossings_{0U};
  };
  interior find_interior(std::vector<bool> const& is_outline) const;

  ways& w_;
  area_routing_options opt_;
  interior_entries interior_;
  std::vector<way_idx_t> way_areas_;
  std::vector<std::int8_t> way_area_layers_;
  std::vector<relation_area> relation_areas_;

  // R-tree over the bounding boxes of the areas (index_areas).
  void* area_index_{nullptr};

  std::mutex obstacle_mutex_;
  std::vector<obstacle> obstacles_;
  std::vector<bool> off_ground_;
  std::vector<bool> connector_;
  std::vector<osm_node_idx_t> elevator_nodes_;

  // Hashes of the barrier segments (both nodes).
  std::mutex barrier_mutex_;
  std::vector<std::uint64_t> barrier_segments_;

  std::vector<osm_node_idx_t> blocked_nodes_;

  struct junction {
    osm_node_idx_t id_;
    point pos_;
    way_properties p_;
  };
  std::vector<junction> junctions_;
};

}  // namespace osr
