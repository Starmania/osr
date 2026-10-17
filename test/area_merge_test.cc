#include "gtest/gtest.h"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string_view>
#include <vector>

#include "fmt/format.h"

#include "osr/extract/area_merge.h"
#include "osr/extract/area_routing.h"
#include "osr/extract/extract.h"
#include "osr/lookup.h"
#include "osr/routing/parameters.h"
#include "osr/routing/path.h"
#include "osr/routing/route.h"
#include "osr/ways.h"
#include "xml_to_pbf.h"

namespace fs = std::filesystem;
using namespace osr;

namespace {

constexpr auto kX0 = 20'000'000;
constexpr auto kY0 = 480'000'000;

struct vertex_spec {
  std::uint64_t id_;
  std::int32_t x_, y_;
  bool entry_{false};
};

using ring_spec = std::vector<vertex_spec>;
using ids = std::vector<std::uint64_t>;
using segments = std::vector<std::pair<std::uint64_t, std::uint64_t>>;

area make_area(std::vector<ring_spec> const& rings) {
  auto a = area{};
  for (auto const& r : rings) {
    auto& ring = a.rings_.emplace_back();
    for (auto const& v : r) {
      ring.push_back({.id_ = v.id_,
                      .pos_ = point{kX0 + v.x_, kY0 + v.y_},
                      .entry_ = v.entry_});
    }
  }
  return a;
}

// Vertex ids of a ring, starting at the smallest id.
ids ring_ids(area::ring const& ring) {
  auto v = ids{};
  for (auto const& x : ring) {
    v.push_back(x.id_);
  }
  std::ranges::rotate(v, std::ranges::min_element(v));
  return v;
}

area_merge_result merge(std::vector<area> const& areas,
                        area_merge_options const& opt = {}) {
  return merge_areas(areas, std::vector<std::uint64_t>(areas.size(), 0U), opt);
}

// 100 x 100 squares on a grid, counter-clockwise from the south west corner:
//
//   8 --- 7 --- 10
//   |  c  |  d  |
//   4 --- 3 --- 6
//   |  a  |  b  |
//   1 --- 2 --- 5
constexpr auto kV1 = vertex_spec{1U, 0, 0};
constexpr auto kV2 = vertex_spec{2U, 100, 0};
constexpr auto kV3 = vertex_spec{3U, 100, 100};
constexpr auto kV4 = vertex_spec{4U, 0, 100};
constexpr auto kV5 = vertex_spec{5U, 200, 0};
constexpr auto kV6 = vertex_spec{6U, 200, 100};
constexpr auto kV7 = vertex_spec{7U, 100, 200};
constexpr auto kV8 = vertex_spec{8U, 0, 200};
constexpr auto kV10 = vertex_spec{10U, 200, 200};

area square_a() { return make_area({{kV1, kV2, kV3, kV4}}); }
area square_b() { return make_area({{kV2, kV5, kV6, kV3}}); }
area square_c() { return make_area({{kV4, kV3, kV7, kV8}}); }
area square_d() { return make_area({{kV3, kV6, kV10, kV7}}); }

}  // namespace

TEST(area_merge, two_squares_sharing_an_edge) {
  auto const r = merge({square_a(), square_b()});
  ASSERT_EQ(1U, r.merged_.size());
  auto const& m = r.merged_[0];
  EXPECT_EQ((std::vector<std::uint32_t>{0U, 1U}), m.members_);
  EXPECT_EQ((std::vector<std::uint32_t>{0U, 0U}), r.group_);
  ASSERT_EQ(1U, m.area_.rings_.size());
  EXPECT_EQ((ids{1U, 2U, 5U, 6U, 3U, 4U}), ring_ids(m.area_.rings_[0]));
  EXPECT_EQ((segments{{2U, 3U}}), m.removed_);
  EXPECT_EQ(1U, r.n_removed_);
  EXPECT_TRUE(r.touch_.empty());
}

TEST(area_merge, winding_order_does_not_matter) {
  // Both outlines run 2 -> 3: the squares are still on opposite sides.
  auto const r =
      merge({square_a(), make_area({{kV2, kV3, kV6, kV5}})});
  ASSERT_EQ(1U, r.merged_.size());
  ASSERT_EQ(1U, r.merged_[0].area_.rings_.size());
  EXPECT_EQ((ids{1U, 2U, 5U, 6U, 3U, 4U}),
            ring_ids(r.merged_[0].area_.rings_[0]));
}

TEST(area_merge, three_squares_in_an_l) {
  // The not merged area in between does not change the result.
  auto const far = make_area({{{20U, 500, 500}, {21U, 600, 500},
                               {22U, 600, 600}}});
  auto const r = merge({square_c(), far, square_b(), square_a()});
  ASSERT_EQ(1U, r.merged_.size());
  auto const& m = r.merged_[0];
  EXPECT_EQ((std::vector<std::uint32_t>{0U, 2U, 3U}), m.members_);
  EXPECT_EQ((std::vector<std::uint32_t>{0U, area_merge_result::kNoGroup, 0U,
                                        0U}),
            r.group_);
  ASSERT_EQ(1U, m.area_.rings_.size());
  // 3 is the reflex corner of the L.
  EXPECT_EQ((ids{1U, 2U, 5U, 6U, 3U, 7U, 8U, 4U}),
            ring_ids(m.area_.rings_[0]));
  EXPECT_EQ((segments{{2U, 3U}, {3U, 4U}}), m.removed_);
  EXPECT_TRUE(r.touch_.empty());
}

TEST(area_merge, four_squares_around_a_node) {
  // The node in the middle is not part of the outline anymore. As an entry,
  // it is kept as a point inside the area.
  for (auto const entry : {false, true}) {
    auto v3 = kV3;
    v3.entry_ = entry;
    auto const r = merge({make_area({{kV1, kV2, v3, kV4}}),
                          make_area({{kV2, kV5, kV6, v3}}),
                          make_area({{kV4, v3, kV7, kV8}}),
                          make_area({{v3, kV6, kV10, kV7}})});
    ASSERT_EQ(1U, r.merged_.size());
    auto const& m = r.merged_[0];
    ASSERT_EQ(entry ? 2U : 1U, m.area_.rings_.size());
    EXPECT_EQ((ids{1U, 2U, 5U, 6U, 10U, 7U, 8U, 4U}),
              ring_ids(m.area_.rings_[0]));
    EXPECT_EQ((segments{{2U, 3U}, {3U, 4U}, {3U, 6U}, {3U, 7U}}), m.removed_);
    if (entry) {
      EXPECT_EQ((ids{3U}), ring_ids(m.area_.rings_[1]));
      EXPECT_TRUE(m.area_.rings_[1][0].entry_);
    }
  }
}

TEST(area_merge, hole_filled_by_another_area) {
  auto const hole = ring_spec{
      {11U, 100, 100}, {12U, 200, 100}, {13U, 200, 200}, {14U, 100, 200}};
  auto const outer = ring_spec{
      {1U, 0, 0}, {2U, 300, 0}, {3U, 300, 300}, {4U, 0, 300}};
  auto const r = merge({make_area({outer, hole}), make_area({hole})});
  ASSERT_EQ(1U, r.merged_.size());
  auto const& m = r.merged_[0];
  ASSERT_EQ(1U, m.area_.rings_.size());
  EXPECT_EQ((ids{1U, 2U, 3U, 4U}), ring_ids(m.area_.rings_[0]));
  EXPECT_EQ(4U, m.removed_.size());
}

TEST(area_merge, hole_partly_filled_by_another_area) {
  // The western half of the hole is filled: the hole shrinks.
  auto const hole = ring_spec{{11U, 100, 100}, {15U, 150, 100},
                              {12U, 200, 100}, {13U, 200, 200},
                              {16U, 150, 200}, {14U, 100, 200}};
  auto const outer = ring_spec{
      {1U, 0, 0}, {2U, 300, 0}, {3U, 300, 300}, {4U, 0, 300}};
  auto const fill = ring_spec{
      {11U, 100, 100}, {15U, 150, 100}, {16U, 150, 200}, {14U, 100, 200}};
  auto const r = merge({make_area({outer, hole}), make_area({fill})});
  ASSERT_EQ(1U, r.merged_.size());
  auto const& m = r.merged_[0];
  ASSERT_EQ(2U, m.area_.rings_.size());
  EXPECT_EQ((ids{1U, 2U, 3U, 4U}), ring_ids(m.area_.rings_[0]));
  // Clockwise.
  EXPECT_EQ((ids{12U, 15U, 16U, 13U}), ring_ids(m.area_.rings_[1]));
}

TEST(area_merge, ring_of_squares_has_a_hole) {
  // 3 x 3 squares without the one in the middle. Vertex ids: 1 + x + 4 * y.
  auto areas = std::vector<area>{};
  for (auto y = 0; y != 3; ++y) {
    for (auto x = 0; x != 3; ++x) {
      if (x == 1 && y == 1) {
        continue;
      }
      auto const v = [](int const vx, int const vy) {
        return vertex_spec{static_cast<std::uint64_t>(1 + vx + 4 * vy),
                           100 * vx, 100 * vy};
      };
      areas.push_back(
          make_area({{v(x, y), v(x + 1, y), v(x + 1, y + 1), v(x, y + 1)}}));
    }
  }
  auto const r = merge(areas);
  ASSERT_EQ(1U, r.merged_.size());
  auto const& m = r.merged_[0];
  EXPECT_EQ(8U, m.members_.size());
  ASSERT_EQ(2U, m.area_.rings_.size());
  EXPECT_EQ((ids{1U, 2U, 3U, 4U, 8U, 12U, 16U, 15U, 14U, 13U, 9U, 5U}),
            ring_ids(m.area_.rings_[0]));
  EXPECT_EQ((ids{6U, 10U, 11U, 7U}), ring_ids(m.area_.rings_[1]));
  EXPECT_EQ(8U, m.removed_.size());
  EXPECT_TRUE(r.touch_.empty());
}

TEST(area_merge, touching_in_one_node_is_not_merged) {
  auto const r = merge({square_a(), square_d()});
  EXPECT_TRUE(r.merged_.empty());
  EXPECT_EQ((std::vector<std::uint32_t>{area_merge_result::kNoGroup,
                                        area_merge_result::kNoGroup}),
            r.group_);
  EXPECT_EQ((ids{3U}), r.touch_);
}

TEST(area_merge, squares_touching_within_a_group) {
  // a and d only touch in 3, but both are merged with b: 3 is not a place
  // where areas only touch, the outline of the union passes it once.
  auto const r = merge({square_a(), square_b(), square_d()});
  ASSERT_EQ(1U, r.merged_.size());
  ASSERT_EQ(1U, r.merged_[0].area_.rings_.size());
  EXPECT_EQ((ids{1U, 2U, 5U, 6U, 10U, 7U, 3U, 4U}),
            ring_ids(r.merged_[0].area_.rings_[0]));
  EXPECT_TRUE(r.touch_.empty());
}

TEST(area_merge, extra_node_on_one_side_is_not_merged) {
  // b has a node between 3 and 2 that a does not have: no common segment.
  // Decision: no geometric matching, valid OSM data shares the nodes. The
  // areas are connected where they touch (2 and 3).
  auto const r = merge(
      {square_a(), make_area({{kV2, kV5, kV6, kV3, {20U, 100, 50}}})});
  EXPECT_TRUE(r.merged_.empty());
  EXPECT_EQ((ids{2U, 3U}), r.touch_);
}

TEST(area_merge, extra_node_on_both_sides_is_merged) {
  auto const mid = vertex_spec{20U, 100, 50};
  auto const r = merge({make_area({{kV1, kV2, mid, kV3, kV4}}),
                        make_area({{kV2, kV5, kV6, kV3, mid}})});
  ASSERT_EQ(1U, r.merged_.size());
  EXPECT_EQ((ids{1U, 2U, 5U, 6U, 3U, 4U}),
            ring_ids(r.merged_[0].area_.rings_[0]));
  EXPECT_EQ((segments{{2U, 20U}, {3U, 20U}}), r.merged_[0].removed_);
}

TEST(area_merge, different_classes_are_not_merged) {
  auto const areas = std::vector<area>{square_a(), square_b(), square_c()};
  auto const r = merge_areas(areas, {7U, 8U, 7U});
  ASSERT_EQ(1U, r.merged_.size());
  EXPECT_EQ((std::vector<std::uint32_t>{0U, 2U}), r.merged_[0].members_);
  EXPECT_EQ(1U, r.n_incompatible_);
  // Shared segments are no touch nodes (3 is on the segment 2 - 3).
  EXPECT_TRUE(r.touch_.empty());

  auto const none = merge_areas(areas, {7U, kNoAreaMerge, kNoAreaMerge});
  EXPECT_TRUE(none.merged_.empty());
  EXPECT_EQ(0U, none.n_incompatible_);
}

TEST(area_merge, blocked_segment_is_not_merged) {
  auto const opt = area_merge_options{
      .is_blocked_ = [](std::uint64_t const a, std::uint64_t const b) {
        return std::min(a, b) == 2U && std::max(a, b) == 3U;
      }};
  auto const r = merge({square_a(), square_b()}, opt);
  EXPECT_TRUE(r.merged_.empty());
  EXPECT_EQ(1U, r.n_blocked_);
  EXPECT_TRUE(r.touch_.empty());

  // Merged through c and d: the blocked segment stays as a wall with two
  // sides, the outline runs 2 -> 3 -> 2.
  auto const g =
      merge({square_a(), square_b(), square_c(), square_d()}, opt);
  ASSERT_EQ(1U, g.merged_.size());
  ASSERT_EQ(1U, g.merged_[0].area_.rings_.size());
  EXPECT_EQ((ids{1U, 2U, 3U, 2U, 5U, 6U, 10U, 7U, 8U, 4U}),
            ring_ids(g.merged_[0].area_.rings_[0]));
  EXPECT_EQ((segments{{3U, 4U}, {3U, 6U}, {3U, 7U}}), g.merged_[0].removed_);
}

TEST(area_merge, overlapping_areas_are_not_merged) {
  // Same side of the segment 2 - 3.
  auto const r = merge(
      {square_a(), make_area({{kV2, kV3, {30U, 50, 150}, {31U, 50, 50}}})});
  EXPECT_TRUE(r.merged_.empty());
  EXPECT_EQ(1U, r.n_conflicts_);

  // Three areas on one segment.
  auto const three = merge({square_a(), square_b(), square_b()});
  EXPECT_TRUE(three.merged_.empty());
  EXPECT_LE(1U, three.n_conflicts_);
}

TEST(area_merge, wall_inside_of_the_union_is_not_merged) {
  // b and c fill the hole of a, the segment between them is blocked: what
  // is left of it is not a ring.
  auto const outer = ring_spec{
      {1U, 0, 0}, {2U, 300, 0}, {3U, 300, 300}, {4U, 0, 300}};
  auto const hole = ring_spec{{11U, 100, 100}, {15U, 150, 100},
                              {12U, 200, 100}, {13U, 200, 200},
                              {16U, 150, 200}, {14U, 100, 200}};
  auto const west = ring_spec{
      {11U, 100, 100}, {15U, 150, 100}, {16U, 150, 200}, {14U, 100, 200}};
  auto const east = ring_spec{
      {15U, 150, 100}, {12U, 200, 100}, {13U, 200, 200}, {16U, 150, 200}};
  auto const r = merge(
      {make_area({outer, hole}), make_area({west}), make_area({east})},
      {.is_blocked_ = [](std::uint64_t const a, std::uint64_t const b) {
        return std::min(a, b) == 15U && std::max(a, b) == 16U;
      }});
  EXPECT_TRUE(r.merged_.empty());
  EXPECT_EQ(1U, r.n_failed_groups_);
  EXPECT_EQ(3U, r.n_failed_areas_);
  EXPECT_EQ(0U, r.n_removed_);
}

TEST(area_merge, merged_area_is_crossed) {
  // Entries: 1, 10 and the node in the middle of four squares.
  auto v1 = kV1;
  auto v3 = kV3;
  auto v10 = kV10;
  v1.entry_ = v3.entry_ = v10.entry_ = true;
  auto const r = merge({make_area({{v1, kV2, v3, kV4}}),
                        make_area({{kV2, kV5, kV6, v3}}),
                        make_area({{kV4, v3, kV7, kV8}}),
                        make_area({{v3, kV6, v10, kV7}})});
  ASSERT_EQ(1U, r.merged_.size());
  auto const& a = r.merged_[0].area_;
  auto const ways = compute_area_ways(a);
  ASSERT_EQ(area_routing_status::kOk, ways.status_);

  // 1 - 3 and 3 - 10: the diagonal 1 - 10 passes through the entry 3.
  auto found = std::vector<ids>{};
  for (auto const& w : ways.ways_) {
    auto& way = found.emplace_back();
    for (auto const& v : w) {
      way.push_back(a.rings_[v.ring_][v.idx_].id_);
    }
    if (way.front() > way.back()) {
      std::ranges::reverse(way);
    }
  }
  std::ranges::sort(found);
  EXPECT_EQ((std::vector<ids>{{1U, 3U}, {3U, 10U}}), found);
}

namespace {

// Three 100m x 100m squares, footways leading away from 1, 6 and 7:
//
//         8 --- 7
//         |  n  |
//   4 --- 3 --- 6
//   |  w  9  e  |
//   1 --- 2 --- 5
constexpr auto const kNodes = R"osm(
  <node id="1" version="1" lat="48.0000" lon="9.00000"/>
  <node id="2" version="1" lat="48.0000" lon="9.00134"/>
  <node id="3" version="1" lat="48.0009" lon="9.00134"/>
  <node id="4" version="1" lat="48.0009" lon="9.00000"/>
  <node id="5" version="1" lat="48.0000" lon="9.00268"/>
  <node id="6" version="1" lat="48.0009" lon="9.00268"/>
  <node id="7" version="1" lat="48.0018" lon="9.00268"/>
  <node id="8" version="1" lat="48.0018" lon="9.00134"/>
  <node id="9" version="1" lat="48.00045" lon="9.00134"/>
  <node id="10" version="1" lat="47.9998" lon="8.99970"/>
  <node id="11" version="1" lat="48.0011" lon="9.00298"/>
  <node id="12" version="1" lat="48.0020" lon="9.00298"/>
  <way id="10" version="1">
    <nd ref="10"/><nd ref="1"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="11" version="1">
    <nd ref="6"/><nd ref="11"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="12" version="1">
    <nd ref="7"/><nd ref="12"/>
    <tag k="highway" v="footway"/>
  </way>
)osm";

constexpr auto const kWest = R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="9"/><nd ref="3"/><nd ref="4"/>
    <nd ref="1"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
  </way>)osm";

constexpr auto const kEast = R"osm(
  <way id="101" version="1">
    <nd ref="2"/><nd ref="5"/><nd ref="6"/><nd ref="3"/><nd ref="9"/>
    <nd ref="2"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
  </way>)osm";

constexpr auto const kNorth = R"osm(
  <way id="102" version="1">
    <nd ref="3"/><nd ref="6"/><nd ref="7"/><nd ref="8"/><nd ref="3"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
  </way>)osm";

constexpr auto const kAt10 = location{47.9998, 8.9997, kNoLevel};
constexpr auto const kAt11 = location{48.0011, 9.00298, kNoLevel};
constexpr auto const kAt12 = location{48.0020, 9.00298, kNoLevel};

struct graph {
  graph(std::string_view const name,
        std::initializer_list<std::string_view> const areas) {
    auto const osm =
        fmt::format("<?xml version=\"1.0\"?>\n<osm version=\"0.6\">{}{}</osm>",
                    kNodes, fmt::join(areas, "\n"));
    auto const osm_path = test::write_osm_pbf(name, osm);
    dir_ = fs::temp_directory_path() / fmt::format("{}_dir", name);

    auto ec = std::error_code{};
    fs::remove_all(dir_, ec);
    fs::create_directories(dir_, ec);

    extract(false, osm_path.generic_string(), dir_, {});
    w_ = std::make_unique<ways>(dir_, cista::mmap::protection::READ);
    l_ = std::make_unique<lookup>(*w_, dir_, cista::mmap::protection::READ);
  }

  path route(location const& from, location const& to) const {
    auto const p =
        osr::route(get_parameters(search_profile::kFoot), *w_, *l_,
                   search_profile::kFoot, from, to, 10'000U,
                   direction::kForward, 50.0, nullptr, nullptr, nullptr,
                   routing_algorithm::kDijkstra);
    EXPECT_TRUE(p.has_value());
    return p.value_or(path{});
  }

  unsigned n_area_ways() const {
    return static_cast<unsigned>(
        std::ranges::count_if(w_->way_osm_idx_, is_area_way));
  }

  bool has_area_way(std::uint64_t const origin) const {
    return w_->find_way(osm_way_idx_t{kAreaWayBit | origin}).has_value();
  }

  bool uses_area_way(path const& p) const {
    return std::ranges::any_of(p.segments_, [&](auto const& s) {
      return s.way_ != way_idx_t::invalid() && s.way_ < w_->n_ways() &&
             is_area_way(w_->way_osm_idx_[s.way_]);
    });
  }

  fs::path dir_;
  std::unique_ptr<ways> w_;
  std::unique_ptr<lookup> l_;
};

}  // namespace

TEST(area_merge, adjacent_areas_are_crossed_diagonally) {
  auto const g = graph{"osr_area_merge_adjacent", {kWest, kEast}};

  // One way 1 - 6, named after the area with the smaller id. Each area on
  // its own has only one entry.
  EXPECT_EQ(1U, g.n_area_ways());
  EXPECT_TRUE(g.has_area_way(100U));
  EXPECT_FALSE(g.has_area_way(101U));

  // 2 x 31m to the corners + 224m diagonal (along the outlines: 300m).
  auto const diagonal = g.route(kAt10, kAt11);
  EXPECT_TRUE(g.uses_area_way(diagonal));
  EXPECT_NEAR(286.0, diagonal.dist_, 4.0);
}

TEST(area_merge, areas_separated_by_a_wall_are_not_merged) {
  auto const g = graph{"osr_area_merge_wall", {kWest, kEast, R"osm(
  <way id="200" version="1">
    <nd ref="2"/><nd ref="9"/><nd ref="3"/>
    <tag k="barrier" v="wall"/>
  </way>)osm"}};
  EXPECT_EQ(0U, g.n_area_ways());
  EXPECT_NEAR(362.0, g.route(kAt10, kAt11).dist_, 4.0);
}

TEST(area_merge, areas_separated_by_rails_are_not_merged) {
  auto const g = graph{"osr_area_merge_rails", {kWest, kEast, R"osm(
  <way id="200" version="1">
    <nd ref="2"/><nd ref="9"/><nd ref="3"/>
    <tag k="railway" v="rail"/>
  </way>)osm"}};
  EXPECT_EQ(0U, g.n_area_ways());
}

TEST(area_merge, kerb_between_areas_does_not_prevent_merging) {
  auto const g = graph{"osr_area_merge_kerb", {kWest, kEast, R"osm(
  <way id="200" version="1">
    <nd ref="2"/><nd ref="9"/><nd ref="3"/>
    <tag k="barrier" v="kerb"/>
  </way>)osm"}};
  EXPECT_EQ(1U, g.n_area_ways());
}

TEST(area_merge, areas_on_different_levels_are_not_merged) {
  auto const g = graph{"osr_area_merge_level", {kWest, R"osm(
  <way id="101" version="1">
    <nd ref="2"/><nd ref="5"/><nd ref="6"/><nd ref="3"/><nd ref="9"/>
    <nd ref="2"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
    <tag k="level" v="1"/>
  </way>)osm"}};
  EXPECT_EQ(0U, g.n_area_ways());
}

TEST(area_merge, multipolygon_is_merged_with_closed_way) {
  // Ways have to be sorted by id.
  auto const g = graph{"osr_area_merge_multipolygon", {R"osm(
  <way id="99" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="9"/><nd ref="3"/>
  </way>
  <way id="100" version="1">
    <nd ref="3"/><nd ref="4"/><nd ref="1"/>
  </way>)osm", kEast, R"osm(
  <relation id="1000" version="1">
    <member type="way" ref="99" role="outer"/>
    <member type="way" ref="100" role="outer"/>
    <tag k="type" v="multipolygon"/>
    <tag k="highway" v="pedestrian"/>
  </relation>)osm"}};

  // Ids of ways generated from closed ways are smaller.
  EXPECT_EQ(1U, g.n_area_ways());
  EXPECT_TRUE(g.has_area_way(101U));

  auto const diagonal = g.route(kAt10, kAt11);
  EXPECT_TRUE(g.uses_area_way(diagonal));
  EXPECT_NEAR(286.0, diagonal.dist_, 4.0);
}

TEST(area_merge, areas_touching_in_a_node_are_crossed_through_it) {
  auto const g = graph{"osr_area_merge_touch", {kWest, kNorth}};

  // 1 - 3 and 3 - 7.
  EXPECT_EQ(2U, g.n_area_ways());
  EXPECT_TRUE(g.has_area_way(100U));
  EXPECT_TRUE(g.has_area_way(102U));

  // 2 x 31m to the corners + 2 x 141m (along the outlines: 400m).
  auto const diagonal = g.route(kAt10, kAt12);
  EXPECT_TRUE(g.uses_area_way(diagonal));
  EXPECT_NEAR(345.0, diagonal.dist_, 4.0);
}

TEST(area_merge, entry_on_removed_segment_stays_connected) {
  // A footway ends at 9, in the middle of the removed segment 2 - 3.
  auto const g = graph{"osr_area_merge_inside", {kWest, kEast, R"osm(
  <node id="20" version="1" lat="48.00045" lon="9.00067"/>
  <way id="200" version="1">
    <nd ref="20"/><nd ref="9"/>
    <tag k="highway" v="footway"/>
  </way>)osm"}};

  // 50m on the footway + 112m from 9 to 6 + 31m (along the outline: 150m).
  auto const p = g.route(location{48.00045, 9.00067, kNoLevel}, kAt11);
  EXPECT_TRUE(g.uses_area_way(p));
  EXPECT_NEAR(193.0, p.dist_, 4.0);
}
