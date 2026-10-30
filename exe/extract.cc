#include <iostream>

#include "fmt/core.h"
#include "fmt/std.h"

#include "conf/options_parser.h"

#include "utl/progress_tracker.h"

#include "osr/extract/extract.h"

using namespace osr;
using namespace boost::program_options;
namespace fs = std::filesystem;

struct config : public conf::configuration {
  config(std::filesystem::path in, std::filesystem::path out)
      : configuration{"Options"}, in_{std::move(in)}, out_{std::move(out)} {
    param(in_, "in,i", "OpenStreetMap .osm.pbf input path");
    param(out_, "out,o", "output directory");
    param(elevation_data_, "elevation_data,e", "directory with elevation data");
    param(with_platforms_, "with_platforms,p", "extract platform info");
    param(areas_, "areas,a",
          "ways across walkable areas: none, spanner (shortest paths between "
          "the access points) or visibility (full visibility graph)");
    param(area_stretch_, "area_stretch",
          "spanner: allowed detour factor between two access points");
    param(area_max_ways_, "area_max_ways",
          "spanner: maximum number of generated ways at one node");
  }

  std::filesystem::path in_, out_, elevation_data_;
  bool with_platforms_{false};
  std::string areas_{"none"};
  double area_stretch_{area_routing_config{}.max_stretch_};
  unsigned area_max_ways_{area_routing_config{}.max_ways_per_node_};
};

int main(int ac, char const** av) {
  auto c = config{"./planet-latest.osm.pbf", "./osr"};

  conf::options_parser parser({&c});
  parser.read_command_line_args(ac, av);

  parser.read_configuration_file();

  parser.print_unrecognized(std::cout);
  parser.print_used(std::cout);

  if (!fs::is_regular_file(c.in_)) {
    fmt::println("input file {} not found", c.in_);
    return 1;
  }

  auto areas = area_routing_config{.max_stretch_ = c.area_stretch_,
                                   .max_ways_per_node_ = c.area_max_ways_};
  if (c.areas_ == "spanner") {
    areas.approach_ = area_routing_config::approach::kSpanner;
  } else if (c.areas_ == "visibility") {
    areas.approach_ = area_routing_config::approach::kVisibilityGraph;
  } else if (c.areas_ != "none") {
    fmt::println("unknown area routing approach {}", c.areas_);
    return 1;
  }

  utl::activate_progress_tracker("osr");
  auto const silencer = utl::global_progress_bars{false};

  extract(c.with_platforms_, c.in_, c.out_, c.elevation_data_, areas);
}