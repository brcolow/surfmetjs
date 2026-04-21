#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <variant>
#include <vector>

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/IO/WKT.h>
#include <CGAL/IO/io.h>
#include <CGAL/Kernel/global_functions_2.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/convex_hull_2.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using FT = K::FT;
using Point = K::Point_2;
using Segment = K::Segment_2;
using Polygon = CGAL::Polygon_2<K>;

using GT = CGAL::Segment_Delaunay_graph_traits_2<K>;
using SDG = CGAL::Segment_Delaunay_graph_2<GT>;
using AT = CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG>;
using VD = CGAL::Voronoi_diagram_2<SDG, AT>;

static std::string decimal_str(const FT &x, int digits = 17) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss.precision(digits);
  oss << CGAL::to_double(x);
  return oss.str();
}

static std::string decimal_double(double x, int digits = 17) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss.precision(digits);
  oss << x;
  return oss.str();
}

static std::string exact_str(const FT &x) {
  std::ostringstream oss;
  CGAL::IO::set_pretty_mode(oss);
  oss << std::setprecision(17);
  oss << CGAL::exact(x);
  return oss.str();
}

static double radius_dec(const FT &r2) {
  return std::sqrt(CGAL::to_double(r2));
}

static std::vector<Point> polygon_vertices(const Polygon &poly) {
  std::vector<Point> points;
  points.reserve(poly.size());
  for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    points.push_back(*v);

  if (points.size() > 1 && points.front() == points.back())
    points.pop_back();

  return points;
}

static std::vector<Point> unique_points(std::vector<Point> points) {
  std::sort(points.begin(), points.end());
  points.erase(std::unique(points.begin(), points.end()), points.end());
  return points;
}

static bool all_collinear(const std::vector<Point> &points) {
  if (points.size() < 3)
    return true;

  std::size_t second = 1;
  while (second < points.size() && points[second] == points[0])
    ++second;
  if (second == points.size())
    return true;

  for (std::size_t i = second + 1; i < points.size(); ++i) {
    if (!CGAL::collinear(points[0], points[second], points[i]))
      return false;
  }

  return true;
}

static FT min_distance2(const std::vector<Point> &points, const Point &center) {
  FT best = CGAL::squared_distance(center, points.front());
  for (std::size_t i = 1; i < points.size(); ++i) {
    const FT d2 = CGAL::squared_distance(center, points[i]);
    if (d2 < best)
      best = d2;
  }
  return best;
}

static FT max_distance2(const std::vector<Point> &points, const Point &center) {
  FT best = CGAL::squared_distance(center, points.front());
  for (std::size_t i = 1; i < points.size(); ++i) {
    const FT d2 = CGAL::squared_distance(center, points[i]);
    if (best < d2)
      best = d2;
  }
  return best;
}

static FT clearance2_edges(const Polygon &poly, const Point &center) {
  const std::size_t n = poly.size();
  FT best =
      CGAL::squared_distance(center, Segment(poly[0], poly[(0 + 1) % n]));
  for (std::size_t i = 1; i < n; ++i) {
    const FT d2 =
        CGAL::squared_distance(center, Segment(poly[i], poly[(i + 1) % n]));
    if (d2 < best)
      best = d2;
  }
  return best;
}

static void print_points_exact(const std::string &label,
                               const std::vector<Point> &points) {
  if (points.empty())
    return;

  std::cout << label << "=";
  for (std::size_t i = 0; i < points.size(); ++i) {
    if (i)
      std::cout << "; ";
    std::cout << "(" << exact_str(points[i].x()) << ", "
              << exact_str(points[i].y()) << ")";
  }
  std::cout << "\n";
}

struct CircleResult {
  Point center;
  FT r2 = FT(0);
  bool valid = false;
  std::vector<Point> support;
};

static CircleResult maximum_inscribed_circle(const Polygon &poly_in) {
  CircleResult out;
  if (poly_in.size() < 3 || !poly_in.is_simple())
    return out;

  Polygon poly = poly_in;
  if (!poly.is_counterclockwise_oriented())
    poly.reverse_orientation();

  SDG sdg;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Point &a = poly[i];
    const Point &b = poly[(i + 1) % poly.size()];
    sdg.insert(a, b);
    sdg.insert(a);
  }

  VD vd(sdg);
  for (auto v = vd.vertices_begin(); v != vd.vertices_end(); ++v) {
    const Point &center = v->point();
    const CGAL::Bounded_side side = poly.bounded_side(center);
    if (side != CGAL::ON_BOUNDED_SIDE && side != CGAL::ON_BOUNDARY)
      continue;

    const FT r2 = clearance2_edges(poly, center);
    if (!out.valid || out.r2 < r2) {
      out.center = center;
      out.r2 = r2;
      out.valid = true;
    }
  }

  if (out.valid) {
    for (std::size_t i = 0; i < poly.size(); ++i) {
      Segment edge(poly[i], poly[(i + 1) % poly.size()]);
      if (CGAL::compare(CGAL::squared_distance(out.center, edge), out.r2) ==
          CGAL::EQUAL) {
        out.support.push_back(poly[i]);
      }
    }
  }

  return out;
}

static CircleResult minimum_circumscribed_circle(const std::vector<Point> &in) {
  CircleResult out;
  const std::vector<Point> points = unique_points(in);
  if (points.empty())
    return out;

  std::vector<Point> hull;
  hull.reserve(points.size());
  CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(hull));
  if (hull.empty())
    return out;

  using Traits = CGAL::Min_circle_2_traits_2<K>;
  using MinCircle = CGAL::Min_circle_2<Traits>;

  MinCircle mc(hull.begin(), hull.end(), true);
  if (mc.is_empty())
    return out;

  const auto circle = mc.circle();
  out.center = circle.center();
  out.r2 = circle.squared_radius();
  out.valid = true;

  for (auto it = mc.support_points_begin(); it != mc.support_points_end(); ++it)
    out.support.push_back(*it);

  return out;
}

struct MZCResult {
  Point center;
  FT r2 = FT(0);
  FT R2 = FT(0);
  FT width = FT(0);
  bool valid = false;
  std::size_t candidates_evaluated = 0;
  std::vector<Point> inner_support;
  std::vector<Point> outer_support;
};

static bool on_nearest_edge(const Point &center, const Point &a,
                            const Point &b,
                            const std::vector<Point> &points) {
  const FT d2a = CGAL::squared_distance(center, a);
  if (CGAL::compare(d2a, CGAL::squared_distance(center, b)) != CGAL::EQUAL)
    return false;

  for (const Point &p : points) {
    if (CGAL::compare(CGAL::squared_distance(center, p), d2a) ==
        CGAL::SMALLER)
      return false;
  }

  return true;
}

static bool on_farthest_edge(const Point &center, const Point &a,
                             const Point &b,
                             const std::vector<Point> &hull) {
  const FT d2a = CGAL::squared_distance(center, a);
  if (CGAL::compare(d2a, CGAL::squared_distance(center, b)) != CGAL::EQUAL)
    return false;

  for (const Point &p : hull) {
    if (CGAL::compare(CGAL::squared_distance(center, p), d2a) == CGAL::LARGER)
      return false;
  }

  return true;
}

static std::optional<Point> bisector_intersection(const Point &a,
                                                  const Point &b,
                                                  const Point &u,
                                                  const Point &v) {
  const FT a1 = FT(2) * (b.x() - a.x());
  const FT b1 = FT(2) * (b.y() - a.y());
  const FT c1 = b.x() * b.x() + b.y() * b.y() -
                a.x() * a.x() - a.y() * a.y();

  const FT a2 = FT(2) * (v.x() - u.x());
  const FT b2 = FT(2) * (v.y() - u.y());
  const FT c2 = v.x() * v.x() + v.y() * v.y() -
                u.x() * u.x() - u.y() * u.y();

  const FT det = a1 * b2 - a2 * b1;
  if (CGAL::compare(det, FT(0)) == CGAL::EQUAL)
    return std::nullopt;

  return Point((c1 * b2 - c2 * b1) / det,
               (a1 * c2 - a2 * c1) / det);
}


static MZCResult minimum_zone_circle_width(const std::vector<Point> &in) {
  const std::vector<Point> points = unique_points(in);
  MZCResult best;

  if (points.empty())
    return best;

  if (points.size() == 1) {
    best.center = points.front();
    best.valid = true;
    best.candidates_evaluated = 1;
    best.inner_support = points;
    best.outer_support = points;
    return best;
  }

  if (points.size() == 2) {
    best.center = CGAL::midpoint(points[0], points[1]);
    best.r2 = CGAL::squared_distance(best.center, points[0]);
    best.R2 = best.r2;
    best.width = FT(0);
    best.valid = true;
    best.candidates_evaluated = 1;
    best.inner_support = points;
    best.outer_support = points;
    return best;
  }

  if (all_collinear(points))
    return best;

  std::vector<Point> hull;
  hull.reserve(points.size());
  CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(hull));
  if (hull.size() < 3)
    return best;

  auto consider = [&](const Point &center) {
    ++best.candidates_evaluated;

    const FT r2 = min_distance2(points, center);
    const FT R2 = max_distance2(hull, center);
    const FT width = CGAL::sqrt(R2) - CGAL::sqrt(r2);

    if (best.valid && !(width < best.width))
      return;

    best.center = center;
    best.r2 = r2;
    best.R2 = R2;
    best.width = width;
    best.valid = true;
    best.inner_support.clear();
    best.outer_support.clear();

    for (const Point &p : points) {
      if (CGAL::compare(CGAL::squared_distance(center, p), r2) ==
          CGAL::EQUAL)
        best.inner_support.push_back(p);
    }
    for (const Point &p : hull) {
      if (CGAL::compare(CGAL::squared_distance(center, p), R2) ==
          CGAL::EQUAL)
        best.outer_support.push_back(p);
    }
  };

  using DT = CGAL::Delaunay_triangulation_2<K>;
  DT dt(points.begin(), points.end());

  std::vector<std::pair<Point, Point>> nearest_pairs;
  nearest_pairs.reserve(static_cast<std::size_t>(3 * points.size()));
  for (auto e = dt.finite_edges_begin(); e != dt.finite_edges_end(); ++e) {
    auto face = e->first;
    const int i = e->second;
    Point a = face->vertex(dt.cw(i))->point();
    Point b = face->vertex(dt.ccw(i))->point();
    if (b < a)
      std::swap(a, b);
    nearest_pairs.emplace_back(a, b);
  }
  std::sort(nearest_pairs.begin(), nearest_pairs.end());
  nearest_pairs.erase(std::unique(nearest_pairs.begin(), nearest_pairs.end()),
                      nearest_pairs.end());

  // 3 inner support points: nearest Voronoi vertices.
  for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end();
       ++face) {
    const Point &a = face->vertex(0)->point();
    const Point &b = face->vertex(1)->point();
    const Point &c = face->vertex(2)->point();
    if (!CGAL::collinear(a, b, c))
      consider(CGAL::circumcenter(a, b, c));
  }

  // 1 inner + 3 outer support points: farthest Voronoi vertices. CGAL does
  // not ship a farthest Voronoi adapter here, so enumerate hull triples and
  // keep only triples whose circumcircle contains the hull.
  for (std::size_t i = 0; i < hull.size(); ++i) {
    for (std::size_t j = i + 1; j < hull.size(); ++j) {
      for (std::size_t k = j + 1; k < hull.size(); ++k) {
        if (CGAL::collinear(hull[i], hull[j], hull[k]))
          continue;

        const Point center = CGAL::circumcenter(hull[i], hull[j], hull[k]);
        const FT R2 = CGAL::squared_distance(center, hull[i]);
        bool is_farthest_vertex = true;
        for (const Point &p : hull) {
          if (CGAL::compare(CGAL::squared_distance(center, p), R2) ==
              CGAL::LARGER) {
            is_farthest_vertex = false;
            break;
          }
        }
        if (is_farthest_vertex)
          consider(center);
      }
    }
  }

  // 2 inner + 2 outer support points: overlay vertices of the nearest and
  // farthest Voronoi diagrams.
  std::vector<std::pair<Point, Point>> farthest_pairs;
  for (std::size_t i = 0; i < hull.size(); ++i) {
    for (std::size_t j = i + 1; j < hull.size(); ++j) {
      Point a = hull[i];
      Point b = hull[j];
      if (b < a)
        std::swap(a, b);
      farthest_pairs.emplace_back(a, b);
    }
  }
  std::sort(farthest_pairs.begin(), farthest_pairs.end());
  farthest_pairs.erase(std::unique(farthest_pairs.begin(), farthest_pairs.end()),
                       farthest_pairs.end());

  for (const auto &near_pair : nearest_pairs) {
    for (const auto &far_pair : farthest_pairs) {
      const auto center = bisector_intersection(
          near_pair.first, near_pair.second, far_pair.first, far_pair.second);
      if (!center)
        continue;

      if (!on_nearest_edge(*center, near_pair.first, near_pair.second, points))
        continue;
      if (!on_farthest_edge(*center, far_pair.first, far_pair.second, hull))
          continue;
      consider(*center);
    }
  }

  return best;
}

static bool verify_mic(const Polygon &poly, const CircleResult &mic) {
  if (!mic.valid || poly.size() < 3 || !poly.is_simple())
    return false;

  if (CGAL::compare(clearance2_edges(poly, mic.center), mic.r2) !=
      CGAL::EQUAL)
    return false;

  SDG sdg;
  for (std::size_t i = 0; i < poly.size(); ++i) {
    const Point &a = poly[i];
    const Point &b = poly[(i + 1) % poly.size()];
    sdg.insert(a, b);
    sdg.insert(a);
  }

  VD vd(sdg);
  for (auto v = vd.vertices_begin(); v != vd.vertices_end(); ++v) {
    const Point &center = v->point();
    const CGAL::Bounded_side side = poly.bounded_side(center);
    if (side != CGAL::ON_BOUNDED_SIDE && side != CGAL::ON_BOUNDARY)
      continue;

    if (CGAL::compare(clearance2_edges(poly, center), mic.r2) ==
        CGAL::LARGER)
      return false;
  }

  return true;
}

static bool verify_mcc(const std::vector<Point> &points,
                       const CircleResult &mcc) {
  if (!mcc.valid || points.empty())
    return false;

  for (const Point &p : points) {
    if (CGAL::compare(CGAL::squared_distance(mcc.center, p), mcc.r2) ==
        CGAL::LARGER)
      return false;
  }

  return !mcc.support.empty() && mcc.support.size() <= 3;
}

static bool verify_mzc(const std::vector<Point> &points,
                       const MZCResult &mzc) {
  if (!mzc.valid || points.empty())
    return false;

  if (CGAL::compare(min_distance2(points, mzc.center), mzc.r2) != CGAL::EQUAL)
    return false;
  if (CGAL::compare(max_distance2(points, mzc.center), mzc.R2) != CGAL::EQUAL)
    return false;

  const FT width = CGAL::sqrt(mzc.R2) - CGAL::sqrt(mzc.r2);
  if (CGAL::compare(width, mzc.width) != CGAL::EQUAL)
    return false;

  return !mzc.inner_support.empty() && !mzc.outer_support.empty();
}

static void print_circle_result(const std::string &label,
                                const CircleResult &circle,
                                const std::string &domain) {
  if (!circle.valid) {
    std::cout << label << "_VALID=0\n";
    return;
  }

  std::cout << label << "_VALID=1\n";
  std::cout << label << "_DOMAIN=" << domain << "\n";
  std::cout << label << "_CENTER_DEC=(" << decimal_str(circle.center.x())
            << ", " << decimal_str(circle.center.y()) << ")\n";
  std::cout << label << "_RADIUS_DEC=" << decimal_double(radius_dec(circle.r2))
            << "\n";
  std::cout << label << "_CENTER_EXACT_X=" << exact_str(circle.center.x())
            << "\n";
  std::cout << label << "_CENTER_EXACT_Y=" << exact_str(circle.center.y())
            << "\n";
  std::cout << label << "_RADIUS2_EXACT=" << exact_str(circle.r2) << "\n";
  std::cout << label << "_RADIUS_EXACT=sqrt(" << exact_str(circle.r2) << ")\n";
  print_points_exact(label + "_SUPPORT_POINTS", circle.support);
}

static void print_mzc_result(const MZCResult &mzc) {
  if (!mzc.valid) {
    std::cout << "MZC_VALID=0\n";
    return;
  }

  std::cout << "MZC_VALID=1\n";
  std::cout << "MZC_DOMAIN=POINT_SET\n";
  std::cout << "MZC_OBJECTIVE=RADIAL_WIDTH\n";
  std::cout << "MZC_SOLVER_EXACT_FOR=FINITE_POINT_SET_RADIAL_WIDTH\n";
  std::cout << "MZC_CENTER_DEC=(" << decimal_str(mzc.center.x()) << ", "
            << decimal_str(mzc.center.y()) << ")\n";
  std::cout << "MZC_INNER_RADIUS_DEC=" << decimal_double(radius_dec(mzc.r2))
            << "\n";
  std::cout << "MZC_OUTER_RADIUS_DEC=" << decimal_double(radius_dec(mzc.R2))
            << "\n";
  std::cout << "MZC_WIDTH_DEC=" << decimal_double(CGAL::to_double(mzc.width))
            << "\n";
  std::cout << "MZC_CENTER_EXACT_X=" << exact_str(mzc.center.x()) << "\n";
  std::cout << "MZC_CENTER_EXACT_Y=" << exact_str(mzc.center.y()) << "\n";
  std::cout << "MZC_R_IN2_EXACT=" << exact_str(mzc.r2) << "\n";
  std::cout << "MZC_R_OUT2_EXACT=" << exact_str(mzc.R2) << "\n";
  std::cout << "MZC_WIDTH_EXACT=sqrt(" << exact_str(mzc.R2) << ")-sqrt("
            << exact_str(mzc.r2) << ")\n";
  std::cout << "MZC_CANDIDATES_EVALUATED=" << mzc.candidates_evaluated
            << "\n";
  print_points_exact("MZC_INNER_SUPPORT", mzc.inner_support);
  print_points_exact("MZC_OUTER_SUPPORT", mzc.outer_support);
}

static void print_mzc_polyline_diagnostic(const Polygon &poly,
                                          const MZCResult &mzc) {
  if (!mzc.valid || poly.size() < 3)
    return;

  const FT r2 = clearance2_edges(poly, mzc.center);
  const FT width = CGAL::sqrt(mzc.R2) - CGAL::sqrt(r2);

  std::cout << "MZC_POLYLINE_DIAGNOSTIC=EVAL_AT_POINT_SET_MZC_CENTER\n";
  std::cout << "MZC_POLYLINE_INNER_RADIUS_DEC="
            << decimal_double(radius_dec(r2)) << "\n";
  std::cout << "MZC_POLYLINE_OUTER_RADIUS_DEC="
            << decimal_double(radius_dec(mzc.R2)) << "\n";
  std::cout << "MZC_POLYLINE_WIDTH_DEC="
            << decimal_double(CGAL::to_double(width)) << "\n";
  std::cout << "MZC_POLYLINE_R_IN2_EXACT=" << exact_str(r2) << "\n";
  std::cout << "MZC_POLYLINE_R_OUT2_EXACT=" << exact_str(mzc.R2) << "\n";
  std::cout << "MZC_POLYLINE_WIDTH_EXACT=sqrt(" << exact_str(mzc.R2)
            << ")-sqrt(" << exact_str(r2) << ")\n";
}

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  CGAL::IO::set_pretty_mode(std::cout);

  std::ifstream fin;
  std::istream *input = &std::cin;
  if (argc > 1) {
    fin.open(argv[1]);
    if (!fin) {
      std::cerr << "Cannot open: " << argv[1] << "\n";
      return 1;
    }
    input = &fin;
  }

  Polygon poly;
  if (!CGAL::IO::read_polygon_WKT(*input, poly)) {
    std::cerr << "Failed to read WKT POLYGON.\n";
    return 1;
  }

  const std::vector<Point> vertices = polygon_vertices(poly);
  const std::vector<Point> points = unique_points(vertices);

  std::cout << "INPUT_VERTEX_COUNT=" << vertices.size() << "\n";
  std::cout << "INPUT_UNIQUE_POINT_COUNT=" << points.size() << "\n";
  std::cout << "INPUT_SIMPLE=" << (poly.is_simple() ? 1 : 0) << "\n";
  std::cout << "INPUT_POINT_SET_SOURCE=POLYGON_VERTICES\n";

  const CircleResult mic = maximum_inscribed_circle(poly);
  const CircleResult mcc = minimum_circumscribed_circle(points);
  const MZCResult mzc = minimum_zone_circle_width(points);

  print_circle_result("MIC", mic, "POLYGON_INTERIOR");
  print_circle_result("MCC", mcc, "POINT_SET_EQUIVALENT_TO_POLYGON");

  if (mic.valid && mcc.valid) {
    const FT dx = mcc.center.x() - mic.center.x();
    const FT dy = mcc.center.y() - mic.center.y();
    const double dr = radius_dec(mcc.r2) - radius_dec(mic.r2);
    const double outer = radius_dec(mcc.r2);
    std::cout << "CENTER_DELTA_DEC=(" << decimal_str(dx) << ", "
              << decimal_str(dy) << ")\n";
    std::cout << "RADIUS_DELTA_DEC=" << decimal_double(dr) << "\n";
    std::cout << "RADIUS_REL_GAP="
              << decimal_double(outer > 0.0 ? dr / outer : 0.0) << "\n";
  }

  print_mzc_result(mzc);
  print_mzc_polyline_diagnostic(poly, mzc);

  std::cout << "MIC_CERT=" << (verify_mic(poly, mic) ? "OK" : "FAIL")
            << "\n";
  std::cout << "MCC_CERT=" << (verify_mcc(points, mcc) ? "OK" : "FAIL")
            << "\n";
  std::cout << "MZC_CERT=" << (verify_mzc(points, mzc) ? "OK" : "FAIL")
            << "\n";

  if (!mzc.valid && points.size() > 2 && all_collinear(points)) {
    std::cerr << "MZC has no finite minimizer for a collinear point set; "
                 "the width infimum is approached at infinity.\n";
  }

  return 0;
}
