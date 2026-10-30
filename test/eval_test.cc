// Evaluation harness (copied into osr/test/ temporarily, never committed).
//
// OSR_DIR=<graph>
//   OSR_STATS=1          graph size: nodes, ways, edges, polyline points
//   OSR_GEN=<n>          print n random foot queries (seeded, OSR_SEED)
//   OSR_QUERIES=<file>   route each query, print distance, cost, time
//   OSR_PATH=1           with OSR_QUERIES: also print the segments (Dijkstra)
#include "gtest/gtest.h"

#include <chrono>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <random>
#include <string>

#include "osr/lookup.h"
#include "osr/routing/parameters.h"
#include "osr/routing/route.h"
#include "osr/ways.h"

#if __has_include("osr/extract/area_ways.h")
#define OSR_EVAL_HAS_AREAS 1
#endif

namespace {

bool is_generated(osr::ways const& w, osr::way_idx_t const way) {
#ifdef OSR_EVAL_HAS_AREAS
  return osr::is_area_way(w.way_osm_idx_[way]);
#else
  (void)w;
  (void)way;
  return false;
#endif
}

}  // namespace

TEST(tmp, eval) {
  auto const* dir = std::getenv("OSR_DIR");
  if (dir == nullptr) {
    return;
  }
  auto w = osr::ways{dir, cista::mmap::protection::READ};
  auto l = osr::lookup{w, dir, cista::mmap::protection::READ};
  std::cout.precision(10);

  if (std::getenv("OSR_STATS") != nullptr) {
    auto n_edges = std::uint64_t{0U}, n_points = std::uint64_t{0U},
         n_gen_ways = std::uint64_t{0U}, n_gen_edges = std::uint64_t{0U},
         n_gen_points = std::uint64_t{0U};
    for (auto i = osr::way_idx_t{0U}; i != w.n_ways(); ++i) {
      auto const e = w.r_->way_nodes_[i].size() > 1U
                         ? w.r_->way_nodes_[i].size() - 1U
                         : 0U;
      auto const p = w.way_polylines_[i].size();
      n_edges += e;
      n_points += p;
      if (is_generated(w, i)) {
        ++n_gen_ways;
        n_gen_edges += e;
        n_gen_points += p;
      }
    }
    std::cout << "STATS nodes=" << w.n_nodes() << " ways=" << w.n_ways()
              << " edges=" << n_edges << " points=" << n_points
              << " gen_ways=" << n_gen_ways << " gen_edges=" << n_gen_edges
              << " gen_points=" << n_gen_points << '\n';
  }

  if (auto const* gen = std::getenv("OSR_GEN"); gen != nullptr) {
    auto const n = std::stoul(gen);
    auto const* seed = std::getenv("OSR_SEED");
    auto rng = std::mt19937_64{seed == nullptr ? 42U : std::stoul(seed)};
    auto way_dist =
        std::uniform_int_distribution<std::uint32_t>{0U, w.n_ways() - 1U};
    auto km = std::uniform_real_distribution<double>{0.2, 4.0};
    auto bearing = std::uniform_real_distribution<double>{0.0, 2.0 * M_PI};
    for (auto i = 0U; i != n;) {
      auto const way = osr::way_idx_t{way_dist(rng)};
      if (is_generated(w, way) ||
          !w.r_->way_properties_[way].is_foot_accessible()) {
        continue;
      }
      auto const& poly = w.way_polylines_[way];
      auto const p =
          poly[std::uniform_int_distribution<std::size_t>{0U, poly.size() - 1U}(
              rng)];
      auto const d = km(rng) * 1000.0, b = bearing(rng);
      auto const lat = p.lat() + (d * std::cos(b)) / 111320.0;
      auto const lng = p.lng() + (d * std::sin(b)) /
                                     (111320.0 * std::cos(p.lat() * M_PI / 180.0));
      std::cout << "Q q" << i << ' ' << p.lat() << ' ' << p.lng() << ' ' << lat
                << ' ' << lng << '\n';
      ++i;
    }
  }

  if (auto const* qf = std::getenv("OSR_QUERIES"); qf != nullptr) {
    auto in = std::ifstream{qf};
    auto tag = std::string{}, id = std::string{};
    double a, b, c, d;
    while (in >> tag >> id >> a >> b >> c >> d) {
      for (auto const algo : {osr::routing_algorithm::kDijkstra,
                              osr::routing_algorithm::kAStarBi}) {
        auto const start = std::chrono::steady_clock::now();
        auto const p = osr::route(
            osr::get_parameters(osr::search_profile::kFoot), w, l,
            osr::search_profile::kFoot, osr::location{{a, b}, osr::kNoLevel},
            osr::location{{c, d}, osr::kNoLevel}, std::chrono::seconds{3600},
            osr::direction::kForward, 100.0, nullptr, nullptr, nullptr, algo);
        auto const us = std::chrono::duration_cast<std::chrono::microseconds>(
                            std::chrono::steady_clock::now() - start)
                            .count();
        auto gen = false;
        if (p.has_value()) {
          for (auto const& s : p->segments_) {
            gen |= s.way_ != osr::way_idx_t::invalid() && s.way_ < w.n_ways() &&
                   is_generated(w, s.way_);
          }
        }
        std::cout << "R " << id << ' '
                  << (algo == osr::routing_algorithm::kDijkstra ? "dij"
                                                                : "astar")
                  << ' ' << (p.has_value() ? p->dist_ : -1.0) << ' '
                  << (p.has_value() ? p->cost_ : 0) << ' ' << gen << ' ' << us
                  << '\n';
        if (std::getenv("OSR_PATH") != nullptr && p.has_value() &&
            algo == osr::routing_algorithm::kDijkstra) {
          for (auto const& s : p->segments_) {
            auto const valid =
                s.way_ != osr::way_idx_t::invalid() && s.way_ < w.n_ways();
            std::cout << "S " << id << ' '
                      << (valid ? cista::to_idx(w.way_osm_idx_[s.way_]) : 0U)
                      << ' ' << (valid && is_generated(w, s.way_)) << ' '
                      << s.dist_ << ' ' << s.cost_ << ' '
                      << s.polyline_.size();
            for (auto const& x : s.polyline_) {
              std::cout << ' ' << x.lat_ << ',' << x.lng_;
            }
            std::cout << '\n';
          }
        }
      }
    }
  }
}
