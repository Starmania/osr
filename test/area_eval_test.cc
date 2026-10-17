// Offline evaluation of compute_area_ways on area dumps (OSR_AREA_DUMP of
// osr-extract). Copied into osr/test/ temporarily, never committed.
//
// OSR_AREAS=<dump>  [OSR_STRETCH=1.02] [OSR_CAP=12] [OSR_LIST=<min ratio>]
//                   [OSR_ID=<area id>] [OSR_MODE=approach|junction]
//
// Per area: ways of the approach and of the optimum (stretch 1, no cap),
// realized entry-to-entry distance of the approach vs the optimum (graph =
// outline + links + generated ways), and the size of the union of the
// optimum's paths when corners are shared junctions.
#include "gtest/gtest.h"

#include <cmath>
#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <fstream>
#include <ranges>
#include <iostream>
#include <map>
#include <queue>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>

#include "geo/latlng.h"

#include "osr/extract/area_routing.h"

namespace {

struct graph {
  std::uint32_t node(std::uint64_t const id) {
    auto const [it, added] =
        idx_.emplace(id, static_cast<std::uint32_t>(out_.size()));
    if (added) {
      out_.emplace_back();
    }
    return it->second;
  }
  void add(std::uint64_t const a, std::uint64_t const b, double const d) {
    auto const x = node(a), y = node(b);
    out_[x].emplace_back(y, d);
    out_[y].emplace_back(x, d);
  }
  std::vector<double> run(std::uint64_t const start) {
    auto const s = node(start);
    auto dist = std::vector<double>(out_.size(), 1e18);
    using l = std::pair<double, std::uint32_t>;
    auto pq = std::priority_queue<l, std::vector<l>, std::greater<>>{};
    dist[s] = 0.0;
    pq.push({0.0, s});
    while (!pq.empty()) {
      auto const [d, u] = pq.top();
      pq.pop();
      if (d > dist[u]) {
        continue;
      }
      for (auto const& [v, w] : out_[u]) {
        if (d + w < dist[v]) {
          dist[v] = d + w;
          pq.push({d + w, v});
        }
      }
    }
    return dist;
  }
  std::unordered_map<std::uint64_t, std::uint32_t> idx_;
  std::vector<std::vector<std::pair<std::uint32_t, double>>> out_;
};

double len(osr::point const a, osr::point const b) {
  return geo::distance(a.as_latlng(), b.as_latlng());
}

// Outline (ring_routable_) + links + generated ways (one edge per way).
graph routing_graph(osr::area const& a, osr::area_routing_result const& r) {
  auto g = graph{};
  for (auto i = 0U; i != a.rings_.size(); ++i) {
    auto const& ring = a.rings_[i];
    if (ring.size() < 2U) {
      continue;
    }
    for (auto j = 0U; j != ring.size(); ++j) {
      auto const& x = ring[j];
      auto const& y = ring[(j + 1U) % ring.size()];
      if (x.id_ != y.id_ && !x.blocked_ && !y.blocked_) {
        g.add(x.id_, y.id_, len(x.pos_, y.pos_));
      }
    }
  }
  for (auto const& l : a.links_) {
    g.add(l.from_, l.to_, l.length_);
  }
  for (auto const& w : r.ways_) {
    auto d = 0.0;
    for (auto i = 1U; i < w.size(); ++i) {
      d += len(get_vertex(a, w[i - 1U]).pos_, get_vertex(a, w[i]).pos_);
    }
    g.add(get_vertex(a, w.front()).id_, get_vertex(a, w.back()).id_, d);
  }
  return g;
}

std::vector<std::uint64_t> entry_ids(osr::area const& a) {
  auto e = std::vector<std::uint64_t>{};
  for (auto const& r : a.rings_) {
    for (auto const& v : r) {
      if (v.entry_ && !v.blocked_) {
        e.push_back(v.id_);
      }
    }
  }
  for (auto const& v : a.interior_) {
    e.push_back(v.id_);
  }
  std::ranges::sort(e);
  e.erase(std::unique(begin(e), end(e)), end(e));
  return e;
}

bool read_area(std::istream& in, std::uint64_t& id, osr::area& a) {
  auto tag = std::string{};
  auto n_rings = 0U;
  if (!(in >> tag >> id >> n_rings)) {
    return false;
  }
  a = {};
  auto const vertex = [&]() {
    auto v = osr::area::vertex{};
    auto e = 0;
    in >> v.id_ >> v.pos_.lat_ >> v.pos_.lng_ >> e;
    v.entry_ = (e & 1) != 0;
    v.blocked_ = (e & 2) != 0;
    return v;
  };
  for (auto i = 0U; i != n_rings; ++i) {
    auto n = 0U;
    in >> tag >> n;
    auto& r = a.rings_.emplace_back();
    for (auto j = 0U; j != n; ++j) {
      r.push_back(vertex());
    }
  }
  auto n = 0U;
  in >> tag >> n;
  for (auto j = 0U; j != n; ++j) {
    a.interior_.push_back(vertex());
  }
  in >> tag >> n;
  for (auto j = 0U; j != n; ++j) {
    auto l = osr::area::link{};
    in >> l.from_ >> l.to_ >> l.length_;
    a.links_.push_back(l);
  }
  return true;
}

double env(char const* name, double const def) {
  auto const* x = std::getenv(name);
  return x == nullptr ? def : std::stod(x);
}

}  // namespace

TEST(tmp, areas) {
  auto const* file = std::getenv("OSR_AREAS");
  if (file == nullptr) {
    return;
  }
  auto opt = osr::area_routing_options{};
  opt.max_stretch_ = env("OSR_STRETCH", 1.02);
  opt.max_ways_per_node_ = static_cast<std::uint32_t>(env("OSR_CAP", 24));
  auto exact = osr::area_routing_options{};
  exact.max_stretch_ = 1.0;
  exact.max_ways_per_node_ = 1000000U;
  auto const list = env("OSR_LIST", 1e9);
  auto const* only_str = std::getenv("OSR_ID");
  auto const only = only_str == nullptr ? -1 : 0;
  auto const only_id = only_str == nullptr ? 0U : std::stoull(only_str);

  auto in = std::ifstream{file};
  auto id = std::uint64_t{};
  auto a = osr::area{};
  auto seen = std::set<std::uint64_t>{};

  auto n_ok = 0U, n_capped_areas = 0U, n_bad_areas = 0U,
       n_bad_capped_areas = 0U;
  auto n_pairs = std::uint64_t{0U}, n_bad_pairs = std::uint64_t{0U},
       n_ways = std::uint64_t{0U}, n_exact_ways = std::uint64_t{0U},
       n_union_segments = std::uint64_t{0U}, n_union_chains = std::uint64_t{0U},
       n_junction_corners = std::uint64_t{0U}, n_capped = std::uint64_t{0U},
       n_junctions = std::uint64_t{0U}, n_way_points = std::uint64_t{0U};
  auto const simulate = env("OSR_SIM", 0) != 0.0;
  auto ratios = std::vector<double>{};
  auto n_j_segments = std::uint64_t{0U}, n_j_chains = std::uint64_t{0U},
       n_j_corners = std::uint64_t{0U}, n_j_over = std::uint64_t{0U};
  auto max_j_degree = 0U;
  auto j_worst = 1.0;
  auto t_approach = 0.0, t_exact = 0.0;
  std::cout.precision(10);

  while (read_area(in, id, a)) {
    if (only >= 0 && id != only_id) {
      continue;
    }
    // Groups that are too large are dumped again as single areas.
    auto const t0 = std::chrono::steady_clock::now();
    auto const r = osr::compute_area_ways(a, opt);
    auto const t1 = std::chrono::steady_clock::now();
    if (r.status_ != osr::area_routing_status::kOk) {
      continue;
    }
    if (std::getenv("OSR_TIMES") != nullptr) {
      auto n_interior = 0U;
      for (auto const& v : a.interior_) {
        n_interior += v.entry_ ? 1U : 0U;
      }
      std::cout << "T " << id << ' ' << r.n_entries_ << ' '
                << a.interior_.size() << ' ' << r.ways_.size() << ' '
                << r.n_visibility_edges_ << ' '
                << std::chrono::duration<double>(t1 - t0).count() << '\n';
      continue;
    }
    auto const x = osr::compute_area_ways(a, exact);
    auto const t2 = std::chrono::steady_clock::now();
    t_approach += std::chrono::duration<double>(t1 - t0).count();
    t_exact += std::chrono::duration<double>(t2 - t1).count();
    ++n_ok;
    n_ways += r.ways_.size();
    n_junctions += r.junctions_.size();
    for (auto const& w : r.ways_) {
      n_way_points += w.size();
    }
    n_exact_ways += x.ways_.size();
    n_capped += r.n_capped_;
    n_capped_areas += r.n_capped_ != 0U ? 1U : 0U;

    if (auto const* gj = std::getenv("OSR_GEOJSON"); gj != nullptr && only >= 0) {
      auto out = std::ofstream{gj};
      out.precision(10);
      auto first = true;
      auto const feature = [&](std::string const& props,
                               std::string const& geom) {
        out << (first ? "" : ",\n") << R"({"type":"Feature","properties":{)"
            << props << R"(},"geometry":)" << geom << '}';
        first = false;
      };
      auto const coord = [](osr::point const p) {
        auto const l = p.as_latlng();
        return "[" + std::to_string(l.lng_) + "," + std::to_string(l.lat_) +
               "]";
      };
      out << R"({"type":"FeatureCollection","features":[)";
      for (auto i = 0U; i != a.rings_.size(); ++i) {
        auto const& ring = a.rings_[i];
        if (ring.size() < 2U) {
          continue;
        }
        auto g = std::string{R"({"type":"LineString","coordinates":[)"};
        for (auto const& v : ring) {
          g += coord(v.pos_) + ",";
        }
        g += coord(ring.front().pos_) + "]}";
        feature(R"("kind":"outline","ring":)" + std::to_string(i) +
                    R"(,"stroke":"#555","stroke-width":2)",
                g);
      }
      auto const ways = [&](osr::area_routing_result const& res,
                            std::string const& kind, std::string const& color) {
        for (auto const& w : res.ways_) {
          auto g = std::string{R"({"type":"LineString","coordinates":[)"};
          for (auto i = 0U; i != w.size(); ++i) {
            g += (i == 0U ? "" : ",") + coord(get_vertex(a, w[i]).pos_);
          }
          g += "]}";
          feature(R"("kind":")" + kind + R"(","from":")" +
                      std::to_string(get_vertex(a, w.front()).id_) +
                      R"(","to":")" +
                      std::to_string(get_vertex(a, w.back()).id_) +
                      R"(","stroke":")" + color + R"(","stroke-width":3)",
                  g);
        }
      };
      ways(x, "exact", "#9c9");
      ways(r, "generated", "#d00");
      auto const point = [&](osr::area::vertex const& v,
                             std::string const& kind,
                             std::string const& color) {
        feature(R"("kind":")" + kind + R"(","node":"n)" +
                    std::to_string(v.id_) + R"(","marker-color":")" + color +
                    '"',
                R"({"type":"Point","coordinates":)" + coord(v.pos_) + "}");
      };
      for (auto const& ring : a.rings_) {
        for (auto const& v : ring) {
          if (v.blocked_) {
            point(v, "blocked", "#000");
          } else if (v.entry_) {
            point(v, "entry", "#0a0");
          } else if (std::ranges::binary_search(r.junctions_, v.id_)) {
            point(v, "junction", "#f80");
          }
        }
      }
      for (auto const& v : a.interior_) {
        point(v, "interior entry", "#0a0");
      }
      out << "]}\n";
    }

    // Union of the exact paths with shared corners.
    auto segs = std::set<std::pair<std::uint64_t, std::uint64_t>>{};
    auto deg = std::map<std::uint64_t, std::uint32_t>{};
    for (auto const& w : x.ways_) {
      for (auto i = 1U; i < w.size(); ++i) {
        auto const p = get_vertex(a, w[i - 1U]).id_;
        auto const q = get_vertex(a, w[i]).id_;
        if (segs.emplace(std::min(p, q), std::max(p, q)).second) {
          ++deg[p];
          ++deg[q];
        }
      }
    }
    auto const entries = entry_ids(a);
    auto chains_ends = 0U;
    for (auto const& [v, d] : deg) {
      auto const is_entry = std::ranges::binary_search(entries, v);
      if (!is_entry && d >= 3U) {
        ++n_junction_corners;
      }
      if (is_entry || d != 2U) {
        chains_ends += d;
      }
    }
    n_union_segments += segs.size();
    n_union_chains += chains_ends / 2U;  // ignores pure cycles

    if (only >= 0) {
      auto per_entry = std::map<std::uint64_t, std::pair<int, int>>{};
      for (auto const& w : r.ways_) {
        ++per_entry[get_vertex(a, w.front()).id_].first;
        ++per_entry[get_vertex(a, w.back()).id_].first;
      }
      for (auto const& w : x.ways_) {
        ++per_entry[get_vertex(a, w.front()).id_].second;
        ++per_entry[get_vertex(a, w.back()).id_].second;
      }
      for (auto const& [e, n] : per_entry) {
        std::cout << "ENTRY " << e << " ways=" << n.first
                  << " exact=" << n.second << '\n';
      }
    }

    // Realized distances.
    auto ga = routing_graph(a, r);
    auto gx = routing_graph(a, x);
    auto worst = 1.0;
    auto bad = 0U;
    for (auto i = 0U; i != entries.size(); ++i) {
      auto const da = ga.run(entries[i]);
      auto const dx = gx.run(entries[i]);
      for (auto j = i + 1U; j != entries.size(); ++j) {
        auto const ia = ga.idx_.find(entries[j]);
        auto const ix = gx.idx_.find(entries[j]);
        if (ia == end(ga.idx_) || ix == end(gx.idx_)) {
          continue;
        }
        auto const opt_d = dx[ix->second];
        auto const app_d = da[ia->second];
        if (opt_d >= 1e17 || opt_d < 0.5) {
          continue;
        }
        auto const ratio = app_d / opt_d;
        ++n_pairs;
        ratios.push_back(ratio);
        worst = std::max(worst, ratio);
        if (ratio > opt.max_stretch_ * 1.0001 + 1e-3 / opt_d) {
          ++bad;
          if (only >= 0) {
            std::cout << "BAD " << entries[i] << ' ' << entries[j] << ' '
                      << app_d << ' ' << opt_d << '\n';
          }
        }
      }
    }
    // Junction spanner: the exact paths go through shared corners. Going
    // through them from short to long, a path is only added if the graph
    // so far (outline + links + added segments) is more than the stretch
    // longer.
    if (simulate) {
      auto j = graph{};
      for (auto const& ring : a.rings_) {
        for (auto k = 0U; ring.size() > 1U && k != ring.size(); ++k) {
          auto const& p = ring[k];
          auto const& q = ring[(k + 1U) % ring.size()];
          if (p.id_ != q.id_ && !p.blocked_ && !q.blocked_) {
            j.add(p.id_, q.id_, len(p.pos_, q.pos_));
          }
        }
      }
      for (auto const& l : a.links_) {
        j.add(l.from_, l.to_, l.length_);
      }
      auto paths = std::vector<std::pair<double, std::size_t>>{};
      for (auto k = 0U; k != x.ways_.size(); ++k) {
        auto d = 0.0;
        auto const& w = x.ways_[k];
        for (auto i = 1U; i < w.size(); ++i) {
          d += len(get_vertex(a, w[i - 1U]).pos_, get_vertex(a, w[i]).pos_);
        }
        paths.emplace_back(d, k);
      }
      std::ranges::sort(paths);
      auto added = std::set<std::pair<std::uint64_t, std::uint64_t>>{};
      for (auto const& [d, k] : paths) {
        auto const& w = x.ways_[k];
        auto const from = get_vertex(a, w.front()).id_;
        auto const to = get_vertex(a, w.back()).id_;
        auto const dist = j.run(from);
        auto const it = j.idx_.find(to);
        if (it != end(j.idx_) &&
            dist[it->second] <= d * opt.max_stretch_ * (1.0 + 1e-5) + 1e-3) {
          continue;
        }
        for (auto i = 1U; i < w.size(); ++i) {
          auto const& p = get_vertex(a, w[i - 1U]);
          auto const& q = get_vertex(a, w[i]);
          if (added.emplace(std::min(p.id_, q.id_), std::max(p.id_, q.id_))
                  .second) {
            j.add(p.id_, q.id_, len(p.pos_, q.pos_));
          }
        }
      }
      auto jdeg = std::map<std::uint64_t, std::uint32_t>{};
      for (auto const& [p, q] : added) {
        ++jdeg[p];
        ++jdeg[q];
      }
      auto ends = 0U;
      for (auto const& [v, dg] : jdeg) {
        auto const is_entry = std::ranges::binary_search(entries, v);
        if (is_entry || dg != 2U) {
          ends += dg;
        }
        if (!is_entry && dg >= 3U) {
          ++n_j_corners;
        }
        max_j_degree = std::max(max_j_degree, dg);
        n_j_over += dg > 12U ? 1U : 0U;
      }
      n_j_segments += added.size();
      n_j_chains += ends / 2U;
      for (auto i = 0U; i != entries.size(); ++i) {
        auto const dj = j.run(entries[i]);
        auto const dx = gx.run(entries[i]);
        for (auto k = i + 1U; k != entries.size(); ++k) {
          auto const ij = j.idx_.find(entries[k]);
          auto const ix = gx.idx_.find(entries[k]);
          if (ix == end(gx.idx_) || dx[ix->second] >= 1e17 ||
              dx[ix->second] < 0.5) {
            continue;
          }
          auto const ratio =
              ij == end(j.idx_) ? 1e9 : dj[ij->second] / dx[ix->second];
          j_worst = std::max(j_worst, ratio);
        }
      }
    }

    n_bad_pairs += bad;
    if (bad != 0U) {
      ++n_bad_areas;
      n_bad_capped_areas += r.n_capped_ != 0U ? 1U : 0U;
    }
    if (worst >= list) {
      auto const& p = a.rings_.front().front().pos_.as_latlng();
      std::cout << "AREA " << id << " entries=" << entries.size()
                << " ways=" << r.ways_.size() << " exact=" << x.ways_.size()
                << " capped=" << r.n_capped_ << " bad=" << bad
                << " worst=" << worst << " at " << p.lat_ << ',' << p.lng_
                << '\n';
    }
  }
  std::ranges::sort(ratios);
  auto const q = [&](double const p) {
    return ratios.empty()
               ? 0.0
               : ratios[std::min(ratios.size() - 1U,
                                 static_cast<std::size_t>(p * ratios.size()))];
  };
  std::cout << "SUMMARY areas=" << n_ok << " ways=" << n_ways
            << " junctions=" << n_junctions << " way_points=" << n_way_points
            << " exact_ways=" << n_exact_ways << " union_segments="
            << n_union_segments << " union_chains=" << n_union_chains
            << " junction_corners=" << n_junction_corners
            << " capped_pairs=" << n_capped
            << " capped_areas=" << n_capped_areas << " pairs=" << n_pairs
            << " bad_pairs=" << n_bad_pairs << " bad_areas=" << n_bad_areas
            << " bad_capped_areas=" << n_bad_capped_areas
            << " p99=" << q(0.99) << " p999=" << q(0.999)
            << " max=" << (ratios.empty() ? 0.0 : ratios.back())
            << " junction: segments=" << n_j_segments
            << " chains=" << n_j_chains << " corners=" << n_j_corners
            << " max_degree=" << max_j_degree << " nodes_over_12=" << n_j_over
            << " worst=" << j_worst << " t_approach=" << t_approach << " t_exact=" << t_exact << '\n';
}
