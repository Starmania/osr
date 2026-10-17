#pragma once

#include <cstdint>
#include <algorithm>
#include <functional>
#include <iterator>

#include "osr/point.h"

namespace osr {

enum class point_location : std::uint8_t { kOutside, kBoundary, kInside };

// Where is p relative to a ring (even-odd rule)?
//
// Exact: only integer orientation tests. The ring may be open or closed
// (first vertex repeated at the end) and may be given in any winding order.
// `proj` maps a ring element to its position (osr::point).
//
// All coordinate differences between p and the ring have to be below 2^31
// (the caller checks the bounding box first), so that the cross products fit
// into 64 bits.
template <typename Ring, typename Proj = std::identity>
point_location locate_in_ring(Ring const& ring, point const p, Proj proj = {}) {
  using coord_t = std::int64_t;
  auto const px = coord_t{p.lat_};
  auto const py = coord_t{p.lng_};

  auto inside = false;
  auto on_boundary = false;
  auto const edge = [&](point const a, point const b) {
    auto const ax = coord_t{a.lat_} - px;
    auto const ay = coord_t{a.lng_} - py;
    auto const bx = coord_t{b.lat_} - px;
    auto const by = coord_t{b.lng_} - py;
    // > 0: p is left of a->b, 0: collinear.
    auto const turn = ax * by - ay * bx;
    if (turn == 0 && std::min(ax, bx) <= 0 && std::max(ax, bx) >= 0 &&
        std::min(ay, by) <= 0 && std::max(ay, by) >= 0) {
      on_boundary = true;
    } else if ((ay > 0) != (by > 0) && ((turn > 0) == (by > ay))) {
      // The edge crosses the horizontal ray that starts at p on the right
      // of p. Edges are half-open: a vertex on the ray belongs to the edge
      // above it, so that every crossing is counted exactly once.
      inside = !inside;
    }
  };

  auto it = std::begin(ring);
  auto const last = std::end(ring);
  if (it == last) {
    return point_location::kOutside;
  }
  auto const first = point{std::invoke(proj, *it)};
  auto prev = first;
  for (++it; it != last; ++it) {
    auto const curr = point{std::invoke(proj, *it)};
    edge(prev, curr);
    prev = curr;
  }
  edge(prev, first);

  return on_boundary ? point_location::kBoundary
                     : (inside ? point_location::kInside
                               : point_location::kOutside);
}

// Where is p relative to a polygon with holes? rings[0] is the outer ring,
// all others are holes. kInside = strictly inside the outer ring and strictly
// outside of all holes. Same preconditions as locate_in_ring.
template <typename Rings, typename Proj = std::identity>
point_location locate_in_rings(Rings const& rings,
                               point const p,
                               Proj proj = {}) {
  auto it = std::begin(rings);
  auto const last = std::end(rings);
  if (it == last) {
    return point_location::kOutside;
  }
  if (auto const outer = locate_in_ring(*it, p, proj);
      outer != point_location::kInside) {
    return outer;
  }
  for (++it; it != last; ++it) {
    switch (locate_in_ring(*it, p, proj)) {
      case point_location::kBoundary: return point_location::kBoundary;
      case point_location::kInside: return point_location::kOutside;
      case point_location::kOutside: break;
    }
  }
  return point_location::kInside;
}

}  // namespace osr
