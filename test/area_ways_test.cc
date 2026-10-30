#include "gtest/gtest.h"

#include <algorithm>
#include <filesystem>
#include <memory>
#include <string_view>

#include "fmt/format.h"

#include "osr/extract/extract.h"
#include "osr/lookup.h"
#include "osr/routing/parameters.h"
#include "osr/routing/path.h"
#include "osr/routing/route.h"
#include "osr/ways.h"
#include "xml_to_pbf.h"

namespace fs = std::filesystem;

namespace {

// 100m x 100m square (nodes 1-4, counter-clockwise from the south west
// corner) with footways leading away from the south west and the north east
// corner. Nodes 30-33 are inside the square (30: center, 31: half way
// between 5 and the center), 15 is on the outline between 1 and 2.
constexpr auto const kNodes = R"osm(
  <node id="1" version="1" lat="48.0000" lon="9.00000"/>
  <node id="2" version="1" lat="48.0000" lon="9.00134"/>
  <node id="3" version="1" lat="48.0009" lon="9.00134"/>
  <node id="4" version="1" lat="48.0009" lon="9.00000"/>
  <node id="5" version="1" lat="48.00045" lon="9.00000"/>
  <node id="6" version="1" lat="48.00045" lon="9.00134"/>
  <node id="10" version="1" lat="47.9998" lon="8.99970"/>
  <node id="11" version="1" lat="48.0011" lon="9.00164"/>
  <node id="12" version="1" lat="48.00045" lon="8.99970"/>
  <node id="13" version="1" lat="48.00045" lon="9.00164"/>
  <node id="14" version="1" lat="47.9998" lon="9.00050"/>
  <node id="15" version="1" lat="48.0000" lon="9.00050"/>
  <node id="16" version="1" lat="47.9998" lon="9.00067"/>
  <node id="17" version="1" lat="48.0000" lon="9.00090">
    <tag k="access" v="private"/>
    <tag k="entrance" v="yes"/>
  </node>
  <node id="20" version="1" lat="48.00027" lon="9.00040"/>
  <node id="21" version="1" lat="48.00027" lon="9.00094"/>
  <node id="22" version="1" lat="48.00063" lon="9.00094"/>
  <node id="23" version="1" lat="48.00063" lon="9.00040"/>
  <node id="30" version="1" lat="48.00045" lon="9.00067"/>
  <node id="31" version="1" lat="48.00045" lon="9.000335"/>
  <node id="32" version="1" lat="48.00030" lon="9.00030"/>
  <node id="33" version="1" lat="48.00060" lon="9.00030"/>
  <node id="40" version="1" lat="48.0020000" lon="9.0000000"/>
  <node id="41" version="1" lat="48.0020000" lon="9.0013425"/>
  <node id="42" version="1" lat="48.0020898" lon="9.0013425"/>
  <node id="43" version="1" lat="48.0022695" lon="9.0013425"/>
  <node id="44" version="1" lat="48.0023593" lon="9.0013425"/>
  <node id="45" version="1" lat="48.0023593" lon="9.0005370"/>
  <node id="46" version="1" lat="48.0028983" lon="9.0005370"/>
  <node id="47" version="1" lat="48.0028983" lon="9.0004027"/>
  <node id="48" version="1" lat="48.0028983" lon="9.0001342"/>
  <node id="49" version="1" lat="48.0028983" lon="9.0000000"/>
  <node id="50" version="1" lat="48.0020898" lon="9.0016110"/>
  <node id="51" version="1" lat="48.0022695" lon="9.0016110"/>
  <node id="52" version="1" lat="48.0030779" lon="9.0001342"/>
  <node id="53" version="1" lat="48.0030779" lon="9.0004027"/>
  <node id="60" version="1" lat="48.0000" lon="9.00067"/>
  <node id="61" version="1" lat="48.0009" lon="9.00067"/>
)osm";

constexpr auto const kFootways = R"osm(
  <way id="10" version="1">
    <nd ref="10"/><nd ref="1"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="11" version="1">
    <nd ref="3"/><nd ref="11"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="12" version="1">
    <nd ref="12"/><nd ref="5"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="13" version="1">
    <nd ref="6"/><nd ref="13"/>
    <tag k="highway" v="footway"/>
  </way>
)osm";

constexpr auto const kSouthWest = osr::location{47.9998, 8.9997, osr::kNoLevel};
constexpr auto const kNorthEast =
    osr::location{48.0011, 9.00164, osr::kNoLevel};
constexpr auto const kWest = osr::location{48.00045, 8.9997, osr::kNoLevel};
constexpr auto const kEast = osr::location{48.00045, 9.00164, osr::kNoLevel};
constexpr auto const kCenter = osr::location{48.00045, 9.00067, osr::kNoLevel};
constexpr auto const kQuarter =
    osr::location{48.00045, 9.000335, osr::kNoLevel};

constexpr auto const kSquare = R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="6"/><nd ref="3"/><nd ref="4"/>
    <nd ref="5"/><nd ref="1"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
  </way>)osm";

struct graph {
  graph(std::string_view const name,
        std::string_view const ways,
        std::string_view const more_ways = "") {
    auto const osm = fmt::format(
        "<?xml version=\"1.0\"?>\n<osm version=\"0.6\">{}{}{}{}</osm>", kNodes,
        kFootways, more_ways, ways);
    auto const osm_path = osr::test::write_osm_pbf(name, osm);
    dir_ = fs::temp_directory_path() / fmt::format("{}_dir", name);

    auto ec = std::error_code{};
    fs::remove_all(dir_, ec);
    fs::create_directories(dir_, ec);

    osr::extract(
        false, osm_path.generic_string(), dir_, {}, {},
        {.approach_ = osr::area_routing_config::approach::kSpanner});
    w_ = std::make_unique<osr::ways>(dir_, cista::mmap::protection::READ);
    l_ = std::make_unique<osr::lookup>(*w_, dir_,
                                       cista::mmap::protection::READ);
  }

  osr::path route(osr::location const& from, osr::location const& to) const {
    auto const p = osr::route(
        osr::get_parameters(osr::search_profile::kFoot), *w_, *l_,
        osr::search_profile::kFoot, from, to, std::chrono::seconds{10'000},
        osr::direction::kForward, 50.0, nullptr, nullptr, nullptr,
        osr::routing_algorithm::kDijkstra);
    EXPECT_TRUE(p.has_value());
    return p.value_or(osr::path{});
  }

  bool is_routing_node(std::uint64_t const osm_node) const {
    return w_->find_node_idx(osr::osm_node_idx_t{osm_node}).has_value();
  }

  unsigned n_area_ways() const {
    return static_cast<unsigned>(
        std::ranges::count_if(w_->way_osm_idx_, osr::is_area_way));
  }

  // Does a generated way end at the node?
  bool has_area_way_at(std::uint64_t const osm_node) const {
    auto const n = w_->find_node_idx(osr::osm_node_idx_t{osm_node});
    return n.has_value() &&
           std::ranges::any_of(w_->r_->node_ways_[*n], [&](auto const way) {
             return osr::is_area_way(w_->way_osm_idx_[way]);
           });
  }

  bool uses_area_way(osr::path const& p) const {
    return std::ranges::any_of(p.segments_, [&](auto const& s) {
      return s.way_ != osr::way_idx_t::invalid() &&
             s.way_ < w_->n_ways() &&
             osr::is_area_way(w_->way_osm_idx_[s.way_]);
    });
  }

  fs::path dir_;
  std::unique_ptr<osr::ways> w_;
  std::unique_ptr<osr::lookup> l_;
};

}  // namespace

TEST(area_ways, closed_way_is_crossed_diagonally) {
  auto const g = graph{"osr_area_ways_square", R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="6"/><nd ref="3"/><nd ref="4"/>
    <nd ref="5"/><nd ref="1"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
  </way>)osm"};

  // Diagonals 1-3 and 5-6. The other pairs are close enough via the outline
  // or one of the diagonals.
  EXPECT_LE(g.n_area_ways(), 4U);

  // 2 x 31m to the corners + 141m diagonal (along the outline: 200m).
  auto const diagonal = g.route(kSouthWest, kNorthEast);
  EXPECT_TRUE(g.uses_area_way(diagonal));
  EXPECT_NEAR(204.0, diagonal.dist_, 4.0);

  // 2 x 22m + 100m straight across (along the outline: 200m).
  auto const straight = g.route(kWest, kEast);
  EXPECT_TRUE(g.uses_area_way(straight));
  EXPECT_NEAR(145.0, straight.dist_, 4.0);

  // The generated ways are found by id, OSM ways still are as well.
  EXPECT_TRUE(g.w_->find_way(osr::osm_way_idx_t{100U}).has_value());
  EXPECT_TRUE(g.w_->find_way(osr::osm_way_idx_t{13U}).has_value());
}

TEST(area_ways, blocked_outline_node_is_not_passed) {
  // Node 17 on the south side is a private entrance: the outline cannot be
  // walked from 1 to 2, the area has to be crossed instead.
  auto const g = graph{"osr_area_ways_blocked", R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="17"/><nd ref="2"/><nd ref="6"/><nd ref="3"/>
    <nd ref="4"/><nd ref="5"/><nd ref="1"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
  </way>
  <way id="20" version="1">
    <nd ref="2"/><nd ref="16"/>
    <tag k="highway" v="footway"/>
  </way>)osm"};

  // 31m to 1, 100m along the south side, 55m from 2 to 16.
  auto const p =
      g.route(kSouthWest, osr::location{47.9998, 9.00067, osr::kNoLevel});
  EXPECT_TRUE(g.uses_area_way(p));
  EXPECT_NEAR(186.0, p.dist_, 6.0);
}

TEST(area_ways, paths_meet_at_a_junction) {
  // L-shaped area (100m x 100m, 60m x 60m cut out in the north east), two
  // entries at the east end (42, 43) and two at the north end (47, 48). All
  // paths between them bend at the reflex corner 45, which becomes a routing
  // node of its own: four ways meet there.
  auto const g = graph{"osr_area_ways_junction", R"osm(
  <way id="400" version="1">
    <nd ref="40"/><nd ref="41"/><nd ref="42"/><nd ref="43"/><nd ref="44"/>
    <nd ref="45"/><nd ref="46"/><nd ref="47"/><nd ref="48"/><nd ref="49"/>
    <nd ref="40"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="area" v="yes"/>
  </way>
  <way id="401" version="1">
    <nd ref="42"/><nd ref="50"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="402" version="1">
    <nd ref="43"/><nd ref="51"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="403" version="1">
    <nd ref="48"/><nd ref="52"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="404" version="1">
    <nd ref="47"/><nd ref="53"/>
    <tag k="highway" v="footway"/>
  </way>)osm"};

  EXPECT_EQ(4U, g.n_area_ways());
  EXPECT_EQ(1, std::ranges::count_if(g.w_->node_to_osm_, osr::is_area_node));
  EXPECT_FALSE(g.is_routing_node(45U));  // the outline does not join

  // 20m + 67m to the corner + 61m to 47 + 20m.
  auto const p =
      g.route(osr::location{48.0020898, 9.0016110, osr::kNoLevel},
              osr::location{48.0030779, 9.0004027, osr::kNoLevel});
  EXPECT_TRUE(g.uses_area_way(p));
  EXPECT_NEAR(168.0, p.dist_, 4.0);
}

TEST(area_ways, area_without_tag_is_not_crossed) {
  auto const g = graph{"osr_area_ways_no_area", R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="6"/><nd ref="3"/><nd ref="4"/>
    <nd ref="5"/><nd ref="1"/>
    <tag k="highway" v="pedestrian"/>
  </way>)osm"};
  EXPECT_EQ(0U, g.n_area_ways());
  EXPECT_NEAR(263.0, g.route(kSouthWest, kNorthEast).dist_, 4.0);
}

TEST(area_ways, platforms_and_other_highways_are_areas) {
  // Same square with the other tags that make an area.
  auto const crossed = [](std::string_view const name,
                          std::string_view const tags) {
    auto const g = graph{name, fmt::format(R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="6"/><nd ref="3"/><nd ref="4"/>
    <nd ref="5"/><nd ref="1"/>
    {}
    <tag k="area" v="yes"/>
  </way>)osm",
                                           tags)};
    return g.n_area_ways() != 0U &&
           std::abs(g.route(kWest, kEast).dist_ - 145.0) < 4.0;
  };
  // railway=platform alone is not walkable in osr (no foot access without
  // public_transport=platform), neither is its outline.
  EXPECT_FALSE(crossed("osr_area_ways_railway_platform",
                       R"(<tag k="railway" v="platform"/>)"));
  EXPECT_TRUE(crossed("osr_area_ways_railway_pt_platform",
                      R"(<tag k="railway" v="platform"/>
    <tag k="public_transport" v="platform"/>)"));
  EXPECT_TRUE(crossed("osr_area_ways_pt_platform",
                      R"(<tag k="public_transport" v="platform"/>)"));
  EXPECT_TRUE(
      crossed("osr_area_ways_service", R"(<tag k="highway" v="service"/>)"));
  EXPECT_TRUE(crossed("osr_area_ways_corridor",
                      R"(<tag k="highway" v="corridor"/>)"));
  EXPECT_FALSE(crossed("osr_area_ways_residential",
                       R"(<tag k="highway" v="residential"/>)"));
}

TEST(area_ways, station_is_no_area) {
  // A station covers all of its levels.
  auto const g = graph{"osr_area_ways_station", R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="6"/><nd ref="3"/><nd ref="4"/>
    <nd ref="5"/><nd ref="1"/>
    <tag k="highway" v="pedestrian"/>
    <tag k="public_transport" v="station"/>
    <tag k="area" v="yes"/>
  </way>)osm"};
  EXPECT_EQ(0U, g.n_area_ways());
}

TEST(area_ways, multipolygon_with_hole) {
  // Hole: 40m x 40m in the center, blocks the straight line from 5 to 6.
  auto const g = graph{"osr_area_ways_multipolygon", R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="6"/><nd ref="3"/>
  </way>
  <way id="101" version="1">
    <nd ref="3"/><nd ref="4"/><nd ref="5"/><nd ref="1"/>
  </way>
  <way id="102" version="1">
    <nd ref="20"/><nd ref="21"/><nd ref="22"/><nd ref="23"/><nd ref="20"/>
  </way>
  <relation id="1000" version="1">
    <member type="way" ref="100" role="outer"/>
    <member type="way" ref="101" role="outer"/>
    <member type="way" ref="102" role="inner"/>
    <tag k="type" v="multipolygon"/>
    <tag k="highway" v="pedestrian"/>
  </relation>)osm"};

  EXPECT_GE(g.n_area_ways(), 1U);

  // 2 x 22m + around the hole: 2 x 36m + 40m (along the outline: 200m).
  auto const around = g.route(kWest, kEast);
  EXPECT_TRUE(g.uses_area_way(around));
  EXPECT_NEAR(157.0, around.dist_, 4.0);
}

TEST(area_ways, footway_ending_inside_is_an_entry) {
  // Footway from the west side to the center of the square, where it ends.
  auto const g = graph{"osr_area_ways_dead_end", kSquare, R"osm(
  <way id="20" version="1">
    <nd ref="5"/><nd ref="31"/><nd ref="30"/>
    <tag k="highway" v="footway"/>
  </way>)osm"};

  // The end of the footway became a routing node, and so did the node on
  // the way to it (every node inside the area is an access point).
  EXPECT_TRUE(g.is_routing_node(30U));
  EXPECT_TRUE(g.is_routing_node(31U));

  // From the center straight to the north east corner: 71m + 31m (back to
  // the outline first and on from there: 50m + 112m + 31m).
  auto const to_corner = g.route(kCenter, kNorthEast);
  EXPECT_TRUE(g.uses_area_way(to_corner));
  EXPECT_NEAR(102.0, to_corner.dist_, 4.0);
  EXPECT_NEAR(102.0, g.route(kNorthEast, kCenter).dist_, 4.0);

  // Straight on to the east side: 50m + 22m.
  EXPECT_NEAR(72.0, g.route(kCenter, kEast).dist_, 4.0);

  // The footway itself connects its end to the outline: 50m + 22m, and
  // crossing the square from west to east did not get longer.
  EXPECT_NEAR(72.0, g.route(kCenter, kWest).dist_, 4.0);
  EXPECT_NEAR(145.0, g.route(kWest, kEast).dist_, 4.0);
  EXPECT_NEAR(204.0, g.route(kSouthWest, kNorthEast).dist_, 4.0);
}

TEST(area_ways, crossing_footway_can_be_left_where_ways_meet) {
  // Two footways cross the square from west to east and meet in the center.
  auto const g = graph{"osr_area_ways_crossing", kSquare, R"osm(
  <way id="20" version="1">
    <nd ref="5"/><nd ref="31"/><nd ref="30"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="21" version="1">
    <nd ref="30"/><nd ref="6"/>
    <tag k="highway" v="footway"/>
  </way>)osm"};
  EXPECT_TRUE(g.is_routing_node(30U));
  EXPECT_TRUE(g.is_routing_node(31U));

  // Starting on the footway, 25m west of the center: leave it right there
  // towards the north east corner: 90m + 31m (via the center: 25m + 71m +
  // 31m, to the outline and on from there: 25m + 112m + 31m).
  auto const leave = g.route(kQuarter, kNorthEast);
  EXPECT_TRUE(g.uses_area_way(leave));
  EXPECT_NEAR(121.0, leave.dist_, 4.0);

  // The footways are not duplicated: only ways from the center to the four
  // corners are added to those of the square without footways.
  auto const plain = graph{"osr_area_ways_crossing_plain", kSquare};
  EXPECT_LE(g.n_area_ways(), plain.n_area_ways() + 4U);
  EXPECT_NEAR(145.0, g.route(kWest, kEast).dist_, 4.0);
}

TEST(area_ways, crossing_footway_can_be_left_at_every_node) {
  // One footway from west to east. Its nodes inside the square are access
  // points: a route that starts on the footway leaves it where it starts.
  auto const g = graph{"osr_area_ways_shape", kSquare, R"osm(
  <way id="20" version="1">
    <nd ref="5"/><nd ref="31"/><nd ref="30"/><nd ref="6"/>
    <tag k="highway" v="footway"/>
  </way>)osm"};
  EXPECT_TRUE(g.is_routing_node(30U));
  EXPECT_TRUE(g.is_routing_node(31U));

  // Straight to the north east corner: 71m + 31m (along the footway to the
  // east side and the outline: 50m + 50m + 31m).
  auto const leave = g.route(kCenter, kNorthEast);
  EXPECT_TRUE(g.uses_area_way(leave));
  EXPECT_NEAR(102.0, leave.dist_, 4.0);
}

TEST(area_ways, ways_not_on_the_area_are_not_connected) {
  // Not connected to the square: a footway that ends on the outline, but
  // not in one of its nodes (15), a tunnel that ends below the center, a
  // footway on another level and one in a lower layer.
  auto const g = graph{"osr_area_ways_not_connected", kSquare, R"osm(
  <way id="20" version="1">
    <nd ref="14"/><nd ref="15"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="21" version="1">
    <nd ref="16"/><nd ref="30"/>
    <tag k="highway" v="footway"/>
    <tag k="tunnel" v="yes"/>
  </way>
  <way id="22" version="1">
    <nd ref="16"/><nd ref="32"/>
    <tag k="highway" v="footway"/>
    <tag k="level" v="-1"/>
  </way>
  <way id="24" version="1">
    <nd ref="16"/><nd ref="31"/>
    <tag k="highway" v="footway"/>
    <tag k="layer" v="-1"/>
  </way>)osm"};
  auto const plain = graph{"osr_area_ways_not_connected_plain", kSquare};
  EXPECT_EQ(plain.n_area_ways(), g.n_area_ways());
  for (auto const n : {15U, 30U, 31U, 32U}) {
    EXPECT_FALSE(g.is_routing_node(n)) << n;
  }
  EXPECT_NEAR(204.0, g.route(kSouthWest, kNorthEast).dist_, 4.0);
}

TEST(area_ways, entry_inside_multipolygon_but_not_in_its_hole) {
  // Same multipolygon as above. Footway 20 ends between the west side and
  // the hole, footway 21 ends in the hole.
  auto const g = graph{"osr_area_ways_multipolygon_inside", R"osm(
  <way id="100" version="1">
    <nd ref="1"/><nd ref="2"/><nd ref="6"/><nd ref="3"/>
  </way>
  <way id="101" version="1">
    <nd ref="3"/><nd ref="4"/><nd ref="5"/><nd ref="1"/>
  </way>
  <way id="102" version="1">
    <nd ref="20"/><nd ref="21"/><nd ref="22"/><nd ref="23"/><nd ref="20"/>
  </way>
  <relation id="1000" version="1">
    <member type="way" ref="100" role="outer"/>
    <member type="way" ref="101" role="outer"/>
    <member type="way" ref="102" role="inner"/>
    <tag k="type" v="multipolygon"/>
    <tag k="highway" v="pedestrian"/>
  </relation>)osm", R"osm(
  <way id="20" version="1">
    <nd ref="5"/><nd ref="31"/>
    <tag k="highway" v="footway"/>
  </way>
  <way id="21" version="1">
    <nd ref="16"/><nd ref="30"/>
    <tag k="highway" v="footway"/>
  </way>)osm"};
  EXPECT_TRUE(g.is_routing_node(31U));
  EXPECT_FALSE(g.is_routing_node(30U));

  // From the end of footway 20 around the hole to the east side: 23m to
  // the corner of the hole, 40m along it, 36m to node 6, 22m (via node 5:
  // 25m more).
  auto const around = g.route(kQuarter, kEast);
  EXPECT_TRUE(g.uses_area_way(around));
  EXPECT_NEAR(121.0, around.dist_, 4.0);
}

TEST(area_ways, building_on_the_square_is_a_hole) {
  // Same hole as in multipolygon_with_hole, mapped as a building.
  auto const building = [](std::string_view const name,
                           std::string_view const tags) {
    return graph{name, kSquare, fmt::format(R"osm(
  <way id="102" version="1">
    <nd ref="20"/><nd ref="21"/><nd ref="22"/><nd ref="23"/><nd ref="20"/>
    {}
  </way>)osm",
                                            tags)};
  };

  // 2 x 22m + around the building: 2 x 36m + 40m.
  auto const g = building("osr_area_ways_building",
                          R"(<tag k="building" v="yes"/>)");
  auto const around = g.route(kWest, kEast);
  EXPECT_TRUE(g.uses_area_way(around));
  EXPECT_NEAR(157.0, around.dist_, 4.0);

  // Not in the way: a roof, a building on pillars, on another level.
  for (auto const& [name, tags] :
       {std::pair{"osr_area_ways_roof", R"(<tag k="building" v="roof"/>)"},
        std::pair{"osr_area_ways_pillars",
                  R"(<tag k="building" v="yes"/>
                     <tag k="building:min_level" v="1"/>)"},
        std::pair{"osr_area_ways_level",
                  R"(<tag k="building" v="yes"/><tag k="level" v="1"/>)"}}) {
    EXPECT_NEAR(145.0, building(name, tags).route(kWest, kEast).dist_, 4.0)
        << name;
  }
}

TEST(area_ways, fence_across_the_square) {
  // Fence from the south to the north side (not connected to the outline).
  constexpr auto kFence = R"osm(
  <way id="200" version="1">
    <nd ref="60"/><nd ref="30"/><nd ref="61"/>
    <tag k="barrier" v="fence"/>
    {}
  </way>)osm";

  // Along the outline: 2 x 22m + 50m + 100m + 50m.
  auto const g =
      graph{"osr_area_ways_fence", kSquare, fmt::format(kFence, "")};
  EXPECT_NEAR(244.0, g.route(kWest, kEast).dist_, 6.0);

  // A footway that passes the fence in the center: straight across.
  auto const gate = graph{"osr_area_ways_fence_gate", kSquare,
                          fmt::format(kFence, "") + R"osm(
  <way id="201" version="1">
    <nd ref="31"/><nd ref="30"/>
    <tag k="highway" v="footway"/>
  </way>)osm"};
  EXPECT_NEAR(145.0, gate.route(kWest, kEast).dist_, 4.0);

  // A fence on a bridge does not block the square.
  auto const bridge =
      graph{"osr_area_ways_fence_bridge", kSquare,
            fmt::format(kFence, R"(<tag k="layer" v="1"/>)")};
  EXPECT_NEAR(145.0, bridge.route(kWest, kEast).dist_, 4.0);
}

TEST(area_ways, steps_and_elevators_are_entries_on_any_level) {
  // Steps from the square (33) down to a tunnel (30, below the center), and
  // a tunnel footway with an elevator (32) that ends in the square.
  auto const g = graph{"osr_area_ways_steps", kSquare, R"osm(
  <node id="70" version="1" lat="48.00060" lon="9.00067">
    <tag k="highway" v="elevator"/>
  </node>
  <way id="20" version="1">
    <nd ref="33"/><nd ref="30"/>
    <tag k="highway" v="steps"/>
    <tag k="level" v="-1;0"/>
  </way>
  <way id="21" version="1">
    <nd ref="30"/><nd ref="16"/>
    <tag k="highway" v="footway"/>
    <tag k="tunnel" v="yes"/>
    <tag k="level" v="-1"/>
  </way>
  <way id="22" version="1">
    <nd ref="16"/><nd ref="70"/><nd ref="32"/>
    <tag k="highway" v="footway"/>
    <tag k="tunnel" v="yes"/>
    <tag k="level" v="-1"/>
  </way>)osm"};
  EXPECT_TRUE(g.has_area_way_at(33U));  // top of the steps
  EXPECT_FALSE(g.has_area_way_at(30U));  // bottom, joined by the tunnel
  EXPECT_TRUE(g.has_area_way_at(70U));  // the elevator
  EXPECT_FALSE(g.has_area_way_at(32U));  // the tunnel itself
}
