#include "gtest/gtest.h"

#include <cmath>
#include <algorithm>
#include <random>
#include <vector>

#include "osr/extract/point_in_polygon.h"

using namespace osr;

namespace {

constexpr auto kX0 = 20'000'000;  // 2.0 E
constexpr auto kY0 = 480'000'000;  // 48.0 N

point pt(std::int32_t const x, std::int32_t const y) {
  return point{kX0 + x, kY0 + y};
}

using ring = std::vector<point>;

constexpr auto kIn = point_location::kInside;
constexpr auto kOn = point_location::kBoundary;
constexpr auto kOut = point_location::kOutside;

}  // namespace

TEST(point_in_polygon, square) {
  auto const r = ring{pt(0, 0), pt(10, 0), pt(10, 10), pt(0, 10)};
  EXPECT_EQ(kIn, locate_in_ring(r, pt(5, 5)));
  EXPECT_EQ(kIn, locate_in_ring(r, pt(1, 9)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(11, 5)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(-1, 5)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(5, 11)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(5, -1)));
}

TEST(point_in_polygon, on_vertex_and_on_edge) {
  auto const r = ring{pt(0, 0), pt(10, 0), pt(10, 10), pt(0, 10)};
  for (auto const& v : r) {
    EXPECT_EQ(kOn, locate_in_ring(r, v));
  }
  EXPECT_EQ(kOn, locate_in_ring(r, pt(5, 0)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(10, 5)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(5, 10)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(0, 5)));

  // On a slanted edge, and right next to it.
  auto const t = ring{pt(0, 0), pt(9, 0), pt(0, 6)};
  EXPECT_EQ(kOn, locate_in_ring(t, pt(3, 4)));
  EXPECT_EQ(kOn, locate_in_ring(t, pt(6, 2)));
  EXPECT_EQ(kIn, locate_in_ring(t, pt(5, 2)));
  EXPECT_EQ(kOut, locate_in_ring(t, pt(7, 2)));
  EXPECT_EQ(kOut, locate_in_ring(t, pt(4, 4)));
}

TEST(point_in_polygon, winding_order_and_closed_rings) {
  auto const ccw = ring{pt(0, 0), pt(10, 0), pt(10, 10), pt(0, 10)};
  auto const cw = ring{pt(0, 10), pt(10, 10), pt(10, 0), pt(0, 0)};
  auto const closed =
      ring{pt(0, 0), pt(10, 0), pt(10, 10), pt(0, 10), pt(0, 0)};
  for (auto x = -2; x != 13; ++x) {
    for (auto y = -2; y != 13; ++y) {
      EXPECT_EQ(locate_in_ring(ccw, pt(x, y)), locate_in_ring(cw, pt(x, y)));
      EXPECT_EQ(locate_in_ring(ccw, pt(x, y)),
                locate_in_ring(closed, pt(x, y)));
    }
  }
}

TEST(point_in_polygon, ray_through_vertices) {
  // The ray to the right of the tested points passes through vertices: a
  // peak (touching), vertices where the outline continues to the other side
  // (crossing) and horizontal edges.
  auto const r = ring{pt(0, 0),  pt(20, 0),  pt(20, 5), pt(16, 5), pt(16, 9),
                      pt(12, 5), pt(10, 5),  pt(8, 5),  pt(8, 10), pt(0, 10)};
  EXPECT_EQ(kIn, locate_in_ring(r, pt(1, 5)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(9, 5)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(11, 5)));
  EXPECT_EQ(kIn, locate_in_ring(r, pt(14, 5)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(18, 5)));
  EXPECT_EQ(kIn, locate_in_ring(r, pt(15, 7)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(15, 8)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(14, 8)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(21, 5)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(-1, 5)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(-1, 0)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(-1, 10)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(-1, 9)));
  EXPECT_EQ(kIn, locate_in_ring(r, pt(1, 9)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(9, 9)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(15, 9)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(16, 9)));

  // Diamond: the ray from the left passes exactly through the left and the
  // right corner.
  auto const d = ring{pt(0, 5), pt(5, 0), pt(10, 5), pt(5, 10)};
  EXPECT_EQ(kOut, locate_in_ring(d, pt(-3, 5)));
  EXPECT_EQ(kOn, locate_in_ring(d, pt(0, 5)));
  EXPECT_EQ(kIn, locate_in_ring(d, pt(3, 5)));
  EXPECT_EQ(kOn, locate_in_ring(d, pt(10, 5)));
  EXPECT_EQ(kOut, locate_in_ring(d, pt(12, 5)));
  EXPECT_EQ(kOut, locate_in_ring(d, pt(-3, 0)));
  EXPECT_EQ(kOut, locate_in_ring(d, pt(-3, 10)));
}

TEST(point_in_polygon, collinear_spike) {
  // Square with a spike of zero width from (10, 5) to (15, 5), and a spike
  // that points into the square from (0, 5) to (4, 5).
  auto const r = ring{pt(0, 0),  pt(10, 0),  pt(10, 5), pt(15, 5),
                      pt(10, 5), pt(10, 10), pt(0, 10), pt(0, 5),
                      pt(4, 5),  pt(0, 5)};
  EXPECT_EQ(kOn, locate_in_ring(r, pt(12, 5)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(15, 5)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(16, 5)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(12, 4)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(12, 6)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(2, 5)));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(4, 5)));
  EXPECT_EQ(kIn, locate_in_ring(r, pt(5, 5)));
  EXPECT_EQ(kIn, locate_in_ring(r, pt(2, 4)));
  EXPECT_EQ(kIn, locate_in_ring(r, pt(2, 6)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(-1, 5)));
  // Collinear with an edge, but not on it.
  EXPECT_EQ(kOut, locate_in_ring(r, pt(10, 12)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(-4, 0)));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(14, 0)));
}

TEST(point_in_polygon, degenerate_rings) {
  EXPECT_EQ(kOut, locate_in_ring(ring{}, pt(0, 0)));
  EXPECT_EQ(kOn, locate_in_ring(ring{pt(1, 1)}, pt(1, 1)));
  EXPECT_EQ(kOut, locate_in_ring(ring{pt(1, 1)}, pt(1, 2)));
  auto const segment = ring{pt(0, 0), pt(10, 10)};
  EXPECT_EQ(kOn, locate_in_ring(segment, pt(5, 5)));
  EXPECT_EQ(kOut, locate_in_ring(segment, pt(5, 6)));
  EXPECT_EQ(kOut, locate_in_ring(segment, pt(11, 11)));
}

TEST(point_in_polygon, holes) {
  auto const rings = std::vector<ring>{
      {pt(0, 0), pt(20, 0), pt(20, 20), pt(0, 20)},
      {pt(5, 5), pt(5, 10), pt(10, 10), pt(10, 5)},
      // Touches the outer ring in (20, 15).
      {pt(15, 12), pt(20, 15), pt(15, 18)}};
  EXPECT_EQ(kIn, locate_in_rings(rings, pt(2, 2)));
  EXPECT_EQ(kIn, locate_in_rings(rings, pt(12, 7)));
  EXPECT_EQ(kOut, locate_in_rings(rings, pt(7, 7)));
  EXPECT_EQ(kOn, locate_in_rings(rings, pt(5, 7)));
  EXPECT_EQ(kOn, locate_in_rings(rings, pt(10, 10)));
  EXPECT_EQ(kOut, locate_in_rings(rings, pt(16, 15)));
  EXPECT_EQ(kOn, locate_in_rings(rings, pt(15, 15)));
  EXPECT_EQ(kOn, locate_in_rings(rings, pt(20, 15)));
  EXPECT_EQ(kOn, locate_in_rings(rings, pt(20, 3)));
  EXPECT_EQ(kOut, locate_in_rings(rings, pt(21, 7)));
  EXPECT_EQ(kOut, locate_in_rings(std::vector<ring>{}, pt(1, 1)));
}

TEST(point_in_polygon, projection) {
  struct vertex {
    int id_;
    point pos_;
  };
  auto const r = std::vector<vertex>{
      {1, pt(0, 0)}, {2, pt(10, 0)}, {3, pt(10, 10)}, {4, pt(0, 10)}};
  EXPECT_EQ(kIn, locate_in_ring(r, pt(5, 5), &vertex::pos_));
  EXPECT_EQ(kOn, locate_in_ring(r, pt(10, 5), &vertex::pos_));
  EXPECT_EQ(kOut, locate_in_ring(r, pt(15, 5), &vertex::pos_));
}

TEST(point_in_polygon, random_polygons_match_winding_number) {
  // Reference: a point that is not on the boundary is inside iff its winding
  // number is odd. Counts the signed crossings of another ray than the
  // implementation (to the top), with doubles for the intersections.
  auto rng = std::mt19937{815U};
  for (auto i = 0U; i != 200U; ++i) {
    // Star-shaped around (0, 0), small coordinates for many special cases.
    auto const n = 3U + rng() % 12U;
    auto r = ring{};
    for (auto j = 0U; j != n; ++j) {
      auto const angle = 2.0 * 3.14159265358979 * (j + 0.5) / n;
      auto const radius = 3.0 + static_cast<double>(rng() % 8U);
      r.push_back(pt(static_cast<std::int32_t>(std::lround(
                         radius * std::cos(angle))),
                     static_cast<std::int32_t>(std::lround(
                         radius * std::sin(angle)))));
    }

    for (auto x = -12; x != 13; ++x) {
      for (auto y = -12; y != 13; ++y) {
        auto const p = pt(x, y);
        auto on_boundary = false;
        auto winding = 0;
        for (auto j = 0U; j != n; ++j) {
          auto const a = r[j];
          auto const b = r[(j + 1U) % n];
          auto const ax = a.lat_ - p.lat_, ay = a.lng_ - p.lng_;
          auto const bx = b.lat_ - p.lat_, by = b.lng_ - p.lng_;
          auto const turn = std::int64_t{ax} * by - std::int64_t{ay} * bx;
          if (turn == 0 && std::min(ax, bx) <= 0 && std::max(ax, bx) >= 0 &&
              std::min(ay, by) <= 0 && std::max(ay, by) >= 0) {
            on_boundary = true;
          }
          // Upward ray from p, edges half-open in x.
          if ((ax <= 0) != (bx <= 0)) {
            auto const y_at_ray =
                ay + static_cast<double>(by - ay) * (0.0 - ax) / (bx - ax);
            if (y_at_ray > 0.0) {
              winding += bx > ax ? 1 : -1;
            }
          }
        }
        auto const expected =
            on_boundary ? kOn : (winding % 2 != 0 ? kIn : kOut);
        EXPECT_EQ(expected, locate_in_ring(r, p))
            << "polygon " << i << ", point " << x << " " << y;
      }
    }
  }
}
