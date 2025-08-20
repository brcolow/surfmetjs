#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <string>
#include <vector>

#include <CGAL/Default.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel_with_sqrt.h>
#include <CGAL/IO/io.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/Segment_Delaunay_graph_2.h>
#include <CGAL/Segment_Delaunay_graph_adaptation_traits_2.h>
#include <CGAL/Segment_Delaunay_graph_traits_2.h>
#include <CGAL/Voronoi_diagram_2.h>
#include <CGAL/convex_hull_2.h>

#include <CGAL/IO/WKT.h>

using K = CGAL::Exact_predicates_exact_constructions_kernel_with_sqrt;
using FT = K::FT;
using Point = K::Point_2;
using Segment = K::Segment_2;
using Polygon = CGAL::Polygon_2<K>;

using GT = CGAL::Segment_Delaunay_graph_traits_2<K>;
using SDG = CGAL::Segment_Delaunay_graph_2<GT>;
using AT = CGAL::Segment_Delaunay_graph_adaptation_traits_2<SDG>;
using VD = CGAL::Voronoi_diagram_2<SDG, AT>;

// Pretty-print an exact number (FT) as rational/surd.
static std::string exact_str(const FT &x) {
  std::ostringstream oss;
  oss << x; // With this Kernel, operator<< prints exact forms (e.g., (a +
            // b*sqrt(c))/d)
  return oss.str();
}

// Decimal for WKT output (WKT is numeric-decimal only).
static std::string decimal_str(const FT &x, int digits = 17) {
  std::ostringstream oss;
  oss.setf(std::ios::fixed);
  oss.precision(digits);
  oss << CGAL::to_double(
      x); // for interoperability; exact form is printed separately
  return oss.str();
}

struct MIC_Result {
  Point center;
  FT r2; // exact squared radius
  bool valid = false;
};

static MIC_Result maximum_inscribed_circle(const Polygon &poly_in) {
  MIC_Result out;
  if (!poly_in.is_simple())
    return out;

  // Ensure CCW orientation for consistency
  Polygon poly = poly_in;
  if (!poly.is_counterclockwise_oriented())
    poly.reverse_orientation();

  const std::size_t n = poly.size();
  if (n < 3)
    return out;

  // Build SDG from polygon segments (+ optional point sites at vertices)
  SDG sdg;
  for (std::size_t i = 0; i < n; ++i) {
    const Point &a = poly[i];
    const Point &b = poly[(i + 1) % n];
    sdg.insert(a, b); // segment site
    sdg.insert(a);    // point site (helps in corner cases)
  }

  VD vd(sdg);

  // Precompute segments for exact distance queries
  std::vector<Segment> edges;
  edges.reserve(n);
  for (std::size_t i = 0; i < n; ++i)
    edges.emplace_back(poly[i], poly[(i + 1) % n]);

  auto clearance2 = [&](const Point &c) -> FT {
    FT best = CGAL::squared_distance(c, edges[0]);
    for (std::size_t i = 1; i < n; ++i) {
      const FT d2 = CGAL::squared_distance(c, edges[i]);
      if (d2 < best)
        best = d2;
    }
    return best;
  };

  // Scan Voronoi vertices inside the polygon
  for (VD::Vertex_iterator v = vd.vertices_begin(); v != vd.vertices_end();
       ++v) {
    const Point &c = v->point();
    auto side = poly.bounded_side(c);
    if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY) {
      FT r2 = clearance2(c);
      if (!out.valid || r2 > out.r2) {
        out.center = c;
        out.r2 = r2;
        out.valid = true;
      }
    }
  }

  return out;
}

struct MCC_Result {
  Point center;
  FT r2; // exact squared radius
  bool valid = false;
  std::vector<Point> support;
};

MCC_Result minimum_circumscribed_circle(const Polygon &poly_in) {
  MCC_Result out;
  if (poly_in.size() == 0)
    return out;

  // 1) Collect vertices (Polygon_2 stores them without duplicate closure)
  std::vector<Point> pts;
  pts.reserve(poly_in.size());
  for (auto v = poly_in.vertices_begin(); v != poly_in.vertices_end(); ++v)
    pts.push_back(*v);

  // 2) Convex hull (not strictly required, but faster for large n)
  std::vector<Point> hull;
  hull.reserve(pts.size());
  CGAL::convex_hull_2(pts.begin(), pts.end(), std::back_inserter(hull));
  if (hull.empty())
    return out;

  // 3) Minimum enclosing circle on the hull
  using Traits = CGAL::Min_circle_2_traits_2<K>;
  using Min_circle = CGAL::Min_circle_2<Traits>;
  Min_circle mc(hull.begin(), hull.end(), true); // randomize input

  if (mc.is_empty())
    return out;
  auto c = mc.circle(); // K::Circle_2
  out.center = c.center();
  out.r2 = c.squared_radius();
  out.valid = true;

  // Record support set: 2 or 3 points define the circle
  for (auto it = mc.support_points_begin(); it != mc.support_points_end(); ++it)
    out.support.push_back(*it);

  return out;
}

struct MZC_Result {
  Point center;
  FT r2, R2;
  bool valid = false;
  std::vector<Point> inner_support, outer_support;
};

static MZC_Result minimum_zone_circle_exact(const Polygon &poly_in) {
  std::vector<Point> P;
  P.reserve(poly_in.size());
  for (auto v = poly_in.vertices_begin(); v != poly_in.vertices_end(); ++v)
    P.push_back(*v);
  MZC_Result best;
  if (P.size() < 2)
    return best;

  // Convex hull vertices (for farthest)
  std::vector<Point> H;
  CGAL::convex_hull_2(P.begin(), P.end(), std::back_inserter(H));
  if (H.size() < 1)
    return best;

  // Delaunay for nearest-VD edges
  using DT = CGAL::Delaunay_triangulation_2<K>;
  DT dt(P.begin(), P.end());

  // Collect unique Delaunay edges (unordered pairs)
  struct PairLess {
    bool operator()(const std::pair<Point, Point> &x,
                    const std::pair<Point, Point> &y) const {
      if (x.first < y.first)
        return true;
      if (y.first < x.first)
        return false;
      return x.second < y.second;
    }
  };
  std::set<std::pair<Point, Point>, PairLess> nearest_pairs;
  for (auto e = dt.finite_edges_begin(); e != dt.finite_edges_end(); ++e) {
    auto fh = e->first;
    int i = e->second;
    Point a = fh->vertex(dt.cw(i))->point();
    Point b = fh->vertex(dt.ccw(i))->point();
    if (b < a)
      std::swap(a, b);
    nearest_pairs.insert({a, b});
  }

  // Farthest-VD edges correspond to adjacent hull pairs
  std::vector<std::pair<Point, Point>> farthest_pairs;
  for (std::size_t i = 0, m = H.size(); i < m; ++i) {
    Point u = H[i], v = H[(i + 1) % m];
    if (v < u)
      std::swap(u, v);
    farthest_pairs.push_back({u, v});
  }

  auto r2_of = [&](const Point &c) {
    FT bestd2 = CGAL::squared_distance(c, P[0]);
    for (std::size_t i = 1; i < P.size(); ++i) {
      FT d2 = CGAL::squared_distance(c, P[i]);
      if (d2 < bestd2)
        bestd2 = d2;
    }
    return bestd2;
  };
  auto R2_of = [&](const Point &c) {
    FT worst = CGAL::squared_distance(c, H[0]);
    for (std::size_t i = 1; i < H.size(); ++i) {
      FT d2 = CGAL::squared_distance(c, H[i]);
      if (worst < d2)
        worst = d2;
    }
    return worst;
  };

  auto on_nearest_edge = [&](const Point &c, const Point &a, const Point &b) {
    FT d2a = CGAL::squared_distance(c, a);
    FT d2b = CGAL::squared_distance(c, b);
    if (CGAL::compare(d2a, d2b) != CGAL::EQUAL)
      return false;
    for (const Point &p : P)
      if (CGAL::compare(CGAL::squared_distance(c, p), d2a) == CGAL::SMALLER)
        return false;
    return true;
  };
  auto on_farthest_edge = [&](const Point &c, const Point &u, const Point &v) {
    FT D2u = CGAL::squared_distance(c, u);
    FT D2v = CGAL::squared_distance(c, v);
    if (CGAL::compare(D2u, D2v) != CGAL::EQUAL)
      return false;
    for (const Point &h : H)
      if (CGAL::compare(CGAL::squared_distance(c, h), D2u) == CGAL::LARGER)
        return false;
    return true;
  };

  auto width_dec = [&](const FT &R2, const FT &r2) {
    // Comparator only; final stored R2/r2 remain exact.
    return std::sqrt(CGAL::to_double(R2)) - std::sqrt(CGAL::to_double(r2));
  };

  auto maybe_update = [&](const Point &c) {
    FT r2 = r2_of(c);
    FT R2 = R2_of(c);
    double w = width_dec(R2, r2);
    if (!best.valid || w < width_dec(best.R2, best.r2)) {
      best.valid = true;
      best.center = c;
      best.r2 = r2;
      best.R2 = R2;
      best.inner_support.clear();
      best.outer_support.clear();
      for (const Point &p : P)
        if (CGAL::compare(CGAL::squared_distance(c, p), r2) == CGAL::EQUAL)
          best.inner_support.push_back(p);
      for (const Point &h : H)
        if (CGAL::compare(CGAL::squared_distance(c, h), R2) == CGAL::EQUAL)
          best.outer_support.push_back(h);
    }
  };

  // 1) Overlay vertices: intersections of nearest-edge bisector with
  // farthest-edge bisector
  for (const auto &ne : nearest_pairs) {
    auto Ln = CGAL::bisector(ne.first, ne.second); // line
    for (const auto &fe : farthest_pairs) {
      auto Lf = CGAL::bisector(fe.first, fe.second);
      auto inter = CGAL::intersection(Ln, Lf);
      if (!inter)
        continue;
      if (const Point *pc = std::get_if<Point>(&*inter)) {
        const Point &c = *pc;
        if (!on_nearest_edge(c, ne.first, ne.second))
          continue;
        if (!on_farthest_edge(c, fe.first, fe.second))
          continue;
        maybe_update(c);
      }
      // (skip coincident lines — degenerate/measure-zero)
    }
  }

  // 2) Nearest-VD vertices (circumcenters of finite Delaunay faces)
  for (auto f = dt.finite_faces_begin(); f != dt.finite_faces_end(); ++f) {
    Point a = f->vertex(0)->point();
    Point b = f->vertex(1)->point();
    Point c = f->vertex(2)->point();
    Point cc = CGAL::circumcenter(a, b, c);
    maybe_update(cc);
  }

  return best;
}

// --- helper: exact squared clearance to polygon boundary (min dist to any
// edge)
static FT clearance2_edges(const Polygon &poly, const Point &c) {
  const std::size_t n = poly.size();
  FT best = CGAL::squared_distance(c, Segment(poly[0], poly[1 % n]));
  for (std::size_t i = 1; i < n; ++i) {
    FT d2 = CGAL::squared_distance(c, Segment(poly[i], poly[(i + 1) % n]));
    if (d2 < best)
      best = d2;
  }
  return best;
}

// --- main method: MZC over the polygonal region (exact r^2 via edges; exact
// R^2 via hull)
static MZC_Result minimum_zone_circle_exact_region(const Polygon &poly) {
  MZC_Result best;
  if (poly.size() < 3 || !poly.is_simple())
    return best;

  // Gather polygon vertices and convex hull (only hull vertices can be
  // farthest)
  std::vector<Point> P;
  P.reserve(poly.size());
  for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    P.push_back(*v);

  std::vector<Point> H;
  CGAL::convex_hull_2(P.begin(), P.end(), std::back_inserter(H));
  if (H.empty())
    return best;

  auto R2_of = [&](const Point &c) {
    FT worst = CGAL::squared_distance(c, H[0]);
    for (std::size_t i = 1; i < H.size(); ++i) {
      FT d2 = CGAL::squared_distance(c, H[i]);
      if (worst < d2)
        worst = d2;
    }
    return worst;
  };

  // Candidate centers:
  // (A) SDG Voronoi vertices (medial-axis vertices of the polygon boundary)
  std::vector<Point> cands;
  {
    SDG sdg;
    const std::size_t n = poly.size();
    for (std::size_t i = 0; i < n; ++i) {
      const Point &a = poly[i];
      const Point &b = poly[(i + 1) % n];
      sdg.insert(a, b); // boundary segment
      sdg.insert(a);    // also add vertex as a point site
    }
    VD vd(sdg);
    for (auto v = vd.vertices_begin(); v != vd.vertices_end(); ++v) {
      const Point &c = v->point();
      auto side = poly.bounded_side(c);
      if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY)
        cands.push_back(c);
    }
  }

  // (B) Farthest-VD vertices ≈ circumcenters of hull triples (filter by
  // farthest check)
  {
    const std::size_t h = H.size();
    if (h >= 3) {
      for (std::size_t i = 0; i < h; ++i)
        for (std::size_t j = i + 1; j < h; ++j)
          for (std::size_t k = j + 1; k < h; ++k) {
            if (CGAL::collinear(H[i], H[j], H[k]))
              continue;
            Point cc = CGAL::circumcenter(H[i], H[j], H[k]);
            // keep it only if those three are co-farthest among H
            FT D2 = CGAL::squared_distance(cc, H[i]);
            if (CGAL::compare(D2, CGAL::squared_distance(cc, H[j])) !=
                CGAL::EQUAL)
              continue;
            if (CGAL::compare(D2, CGAL::squared_distance(cc, H[k])) !=
                CGAL::EQUAL)
              continue;
            bool ok = true;
            for (const Point &hpt : H) {
              if (CGAL::compare(CGAL::squared_distance(cc, hpt), D2) ==
                  CGAL::LARGER) {
                ok = false;
                break;
              }
            }
            if (ok) {
              auto side = poly.bounded_side(cc);
              if (side == CGAL::ON_BOUNDED_SIDE || side == CGAL::ON_BOUNDARY)
                cands.push_back(cc);
            }
          }
    }
  }

  // Evaluate all candidates exactly; choose the one minimizing (sqrt(R2) -
  // sqrt(r2)). We store exact R2/r2; for ordering we use high-precision doubles
  // (robust in practice).
  auto consider = [&](const Point &c) {
    FT r2 = clearance2_edges(poly, c); // min distance to boundary (edges)
    FT R2 = R2_of(c);                  // max distance to hull vertices
    double width =
        std::sqrt(CGAL::to_double(R2)) - std::sqrt(CGAL::to_double(r2));
    if (!best.valid || width < (std::sqrt(CGAL::to_double(best.R2)) -
                                std::sqrt(CGAL::to_double(best.r2)))) {
      best.valid = true;
      best.center = c;
      best.r2 = r2;
      best.R2 = R2;
      best.inner_support.clear();
      best.outer_support.clear();
      const std::size_t n = poly.size();
      for (std::size_t i = 0; i < n; ++i) {
        Segment e(poly[i], poly[(i + 1) % n]);
        if (CGAL::compare(CGAL::squared_distance(c, e), r2) == CGAL::EQUAL)
          best.inner_support.push_back(
              poly[i]); // store an endpoint for reference
      }
      for (const Point &hv : H)
        if (CGAL::compare(CGAL::squared_distance(c, hv), R2) == CGAL::EQUAL)
          best.outer_support.push_back(hv);
    }
  };

  for (const Point &c : cands)
    consider(c);

  return best;
}

// Verifies MCC against all polygon vertices (exact).
// Returns true if: every vertex is inside/on, and there are 2 or 3 exact
// supporters.
static bool verify_mcc(const Polygon &poly, const Point &C, const FT &R2) {
  if (poly.size() == 0)
    return false;

  // All vertices must be <= R2 from center
  for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) {
    FT d2 = CGAL::squared_distance(C, *v);
    if (CGAL::compare(d2, R2) == CGAL::LARGER)
      return false;
  }

  // Count exact boundary supporters (should be 2 or 3 for a tight MEC)
  int on_cnt = 0;
  for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v) {
    FT d2 = CGAL::squared_distance(C, *v);
    if (CGAL::compare(d2, R2) == CGAL::EQUAL)
      ++on_cnt;
  }
  return (on_cnt >= 2 && on_cnt <= 3);
}

// Returns true iff r2_star is the maximum clearance among all SDG Voronoi
// vertices. (This is an exact, sufficiency-based certificate for polygons with
// straight edges.)
static bool verify_mic(const Polygon &poly, const Point &C, const FT &r2_star) {
  if (!poly.is_simple() || poly.size() < 3)
    return false;

  // Build SDG of the boundary (same as your MIC computation)
  SDG sdg;
  const std::size_t n = poly.size();
  for (std::size_t i = 0; i < n; ++i) {
    const Point &a = poly[i];
    const Point &b = poly[(i + 1) % n];
    sdg.insert(a, b);
    sdg.insert(a);
  }
  VD vd(sdg);

  // 1) r2_star must equal clearance at the chosen center
  if (CGAL::compare(clearance2_edges(poly, C), r2_star) != CGAL::EQUAL)
    return false;

  // 2) No Voronoi vertex has strictly larger clearance
  for (auto v = vd.vertices_begin(); v != vd.vertices_end(); ++v) {
    FT r2_here = clearance2_edges(poly, v->point());
    if (CGAL::compare(r2_here, r2_star) == CGAL::LARGER)
      return false;
  }

  return true; // MIC certified
}

struct WidthVal {
  FT r2, R2;
};

static FT R2_from_hull(const std::vector<Point> &H, const Point &c) {
  FT worst = CGAL::squared_distance(c, H[0]);
  for (std::size_t i = 1; i < H.size(); ++i) {
    FT d2 = CGAL::squared_distance(c, H[i]);
    if (worst < d2)
      worst = d2;
  }
  return worst;
}

static bool verify_mzc_region(const Polygon &poly, const Point &Cstar,
                              const FT &r2_star, const FT &R2_star) {
  if (poly.size() < 3 || !poly.is_simple())
    return false;

  // Gather vertices and hull (outer distance only needs hull)
  std::vector<Point> P;
  P.reserve(poly.size());
  for (auto v = poly.vertices_begin(); v != poly.vertices_end(); ++v)
    P.push_back(*v);

  std::vector<Point> H;
  CGAL::convex_hull_2(P.begin(), P.end(), std::back_inserter(H));
  if (H.empty())
    return false;

  // Target width (decimal only for ordering)
  auto width_dec = [](const FT &R2, const FT &r2) {
    return std::sqrt(CGAL::to_double(R2)) - std::sqrt(CGAL::to_double(r2));
  };
  const double w_star = width_dec(R2_star, r2_star);

  // A) SDG vertices of the boundary
  std::vector<Point> cands;
  {
    SDG sdg;
    const std::size_t n = poly.size();
    for (std::size_t i = 0; i < n; ++i) {
      const Point &a = poly[i];
      const Point &b = poly[(i + 1) % n];
      sdg.insert(a, b); // edge site
      sdg.insert(a);    // vertex site (harmless)
    }
    VD vd(sdg);
    for (auto v = vd.vertices_begin(); v != vd.vertices_end(); ++v)
      cands.push_back(
          v->point()); // note: we DO include centers outside; that’s allowed
  }

  // B) Farthest-VD vertices (circumcenters of hull triples that are
  // co-farthest)
  {
    const std::size_t h = H.size();
    if (h >= 3) {
      for (std::size_t i = 0; i < h; ++i)
        for (std::size_t j = i + 1; j < h; ++j)
          for (std::size_t k = j + 1; k < h; ++k) {
            if (CGAL::collinear(H[i], H[j], H[k]))
              continue;
            Point cc = CGAL::circumcenter(H[i], H[j], H[k]);
            FT D2 = CGAL::squared_distance(cc, H[i]);
            if (CGAL::compare(D2, CGAL::squared_distance(cc, H[j])) !=
                CGAL::EQUAL)
              continue;
            if (CGAL::compare(D2, CGAL::squared_distance(cc, H[k])) !=
                CGAL::EQUAL)
              continue;
            bool max_ok = true;
            for (const Point &hpt : H) {
              if (CGAL::compare(CGAL::squared_distance(cc, hpt), D2) ==
                  CGAL::LARGER) {
                max_ok = false;
                break;
              }
            }
            if (max_ok)
              cands.push_back(cc);
          }
    }
  }

  // C) Evaluate every candidate; none may beat (R2*, r2*)
  for (const Point &c : cands) {
    FT r2 = clearance2_edges(poly, c); // inner: edge distance (exact)
    FT R2 = R2_from_hull(H, c);        // outer: hull vertex distance (exact)
    if (width_dec(R2, r2) < w_star)
      return false; // someone beats us → not optimal
  }

  // D) Our own center must realize the stored values exactly
  if (CGAL::compare(clearance2_edges(poly, Cstar), r2_star) != CGAL::EQUAL)
    return false;
  if (CGAL::compare(R2_from_hull(H, Cstar), R2_star) != CGAL::EQUAL)
    return false;

  return true; // MZC certified against the complete region candidate set
               // (A)+(B)
}

int main(int argc, char **argv) {
  std::ios::sync_with_stdio(false);
  CGAL::IO::set_pretty_mode(
      std::cout); // exact fractions/surds when streaming FT

  // High-precision decimal for WKT/diagnostics
  std::cout.setf(std::ios::fixed);
  std::cout << std::setprecision(17);

  // Read WKT POLYGON from file or stdin
  std::ifstream fin;
  std::istream *in = &std::cin;
  if (argc > 1) {
    fin.open(argv[1]);
    if (!fin) {
      std::cerr << "Cannot open: " << argv[1] << "\n";
      return 1;
    }
    in = &fin;
  }

  Polygon poly;
  if (!CGAL::IO::read_polygon_WKT(*in, poly)) {
    std::cerr << "Failed to read WKT POLYGON.\n";
    return 1;
  }

  // MIC
  auto mic = maximum_inscribed_circle(poly);
  if (!mic.valid) {
    std::cerr << "MIC not found (check that the polygon is simple and "
                 "non-degenerate).\n";
  } else {
    std::cout << "MIC_CENTER_DEC=(" << decimal_str(mic.center.x()) << ", "
              << decimal_str(mic.center.y()) << ")\n";
    std::cout << "MIC_RADIUS_DEC=" << std::sqrt(CGAL::to_double(mic.r2))
              << "\n";

    std::cout << "MIC_CENTER_EXACT_X=" << exact_str(mic.center.x()) << "\n";
    std::cout << "MIC_CENTER_EXACT_Y=" << exact_str(mic.center.y()) << "\n";
    std::cout << "MIC_RADIUS2_EXACT=" << exact_str(mic.r2) << "\n";
    std::cout << "MIC_RADIUS_EXACT=sqrt(" << exact_str(mic.r2) << ")\n";
  }

  // MCC
  auto mcc = minimum_circumscribed_circle(poly);
  if (!mcc.valid) {
    std::cerr << "MCC not found (empty/degenerate input?).\n";
  } else {
    std::cout << "MCC_CENTER_DEC=(" << decimal_str(mcc.center.x()) << ", "
              << decimal_str(mcc.center.y()) << ")\n";
    std::cout << "MCC_RADIUS_DEC=" << std::sqrt(CGAL::to_double(mcc.r2))
              << "\n";

    std::cout << "MCC_CENTER_EXACT_X=" << exact_str(mcc.center.x()) << "\n";
    std::cout << "MCC_CENTER_EXACT_Y=" << exact_str(mcc.center.y()) << "\n";
    std::cout << "MCC_RADIUS2_EXACT=" << exact_str(mcc.r2) << "\n";
    std::cout << "MCC_RADIUS_EXACT=sqrt(" << exact_str(mcc.r2) << ")\n";

    if (!mcc.support.empty()) {
      std::cout << "MCC_SUPPORT_POINTS=";
      for (std::size_t i = 0; i < mcc.support.size(); ++i) {
        if (i)
          std::cout << "; ";
        std::cout << "(" << exact_str(mcc.support[i].x()) << ", "
                  << exact_str(mcc.support[i].y()) << ")";
      }
      std::cout << "\n";
    }
  }

  // If both valid, print deltas (decimal) to show they’re close but not
  // identical
  if (mic.valid && mcc.valid) {
    auto dx = decimal_str(mcc.center.x() - mic.center.x());
    auto dy = decimal_str(mcc.center.y() - mic.center.y());
    auto dr =
        std::sqrt(CGAL::to_double(mcc.r2)) - std::sqrt(CGAL::to_double(mic.r2));
    std::cout << "CENTER_DELTA_DEC=(" << dx << "," << dy << ")\n";
    std::cout << "RADIUS_DELTA_DEC=" << dr << "\n";

    // Relative gap (how “annulus-like” the polygon is)
    double r_out = std::sqrt(CGAL::to_double(mcc.r2));
    double rel = (r_out > 0.0) ? (dr / r_out) : 0.0;
    std::cout << "RADIUS_REL_GAP=" << rel << "  # fraction of outer radius\n";
  }

  // MZC
  auto mzc = minimum_zone_circle_exact(poly);
  if (mzc.valid) {
    std::cout << "MZC_CENTER_DEC=(" << decimal_str(mzc.center.x()) << ", "
              << decimal_str(mzc.center.y()) << ")\n";
    std::cout << "MZC_WIDTH_DEC="
              << (std::sqrt(CGAL::to_double(mzc.R2)) -
                  std::sqrt(CGAL::to_double(mzc.r2)))
              << "\n";
    std::cout << "MZC_CENTER_EXACT_X=" << exact_str(mzc.center.x()) << "\n";
    std::cout << "MZC_CENTER_EXACT_Y=" << exact_str(mzc.center.y()) << "\n";
    std::cout << "MZC_R_IN2_EXACT=" << exact_str(mzc.r2) << "\n";
    std::cout << "MZC_R_OUT2_EXACT=" << exact_str(mzc.R2) << "\n";
    if (!mzc.inner_support.empty()) {
      std::cout << "MZC_INNER_SUPPORT=";
      for (std::size_t i = 0; i < mzc.inner_support.size(); ++i) {
        if (i)
          std::cout << "; ";
        std::cout << "(" << exact_str(mzc.inner_support[i].x()) << ", "
                  << exact_str(mzc.inner_support[i].y()) << ")";
      }
      std::cout << "\n";
    }
    if (!mzc.outer_support.empty()) {
      std::cout << "MZC_OUTER_SUPPORT=";
      for (std::size_t i = 0; i < mzc.outer_support.size(); ++i) {
        if (i)
          std::cout << "; ";
        std::cout << "(" << exact_str(mzc.outer_support[i].x()) << ", "
                  << exact_str(mzc.outer_support[i].y()) << ")";
      }
      std::cout << "\n";
    }
  } else {
    std::cerr << "MZC not found.\n";
  }

  /*
  bool mic_ok = verify_mic(poly, mic.center, mic.r2);
  bool mcc_ok = verify_mcc(poly, mcc.center, mcc.r2);
  bool mzc_ok = verify_mzc_region(poly, mzc.center, mzc.r2, mzc.R2);

  std::cerr << "MIC_CERT=" << (mic_ok ? "OK" : "FAIL") << "\n";
  std::cerr << "MCC_CERT=" << (mcc_ok ? "OK" : "FAIL") << "\n";
  std::cerr << "MZC_CERT=" << (mzc_ok ? "OK" : "FAIL") << "\n";
  */
  return 0;
}