#ifndef _DELAUNAYHEADER_
#define _DELAUNAYHEADER_
#endif

#define CGAL_EIGEN3_ENABLED 1

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/Triangulation_face_base_with_info_2.h>

#include <CGAL/Delaunay_triangulation_2.h>

#include <CGAL/Projection_traits_xy_3.h>
#include <CGAL/utility.h>

#include <CGAL/Delaunay_triangulation_3.h>
#include <CGAL/Tetrahedron_3.h>
#include <CGAL/Triangulation_vertex_base_with_info_3.h>

#include <RcppEigen.h>

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

struct FaceInfo2 {
  FaceInfo2() {}
  int nesting_level;
  bool in_domain() { return nesting_level % 2 == 1; }
};
typedef CGAL::Triangulation_vertex_base_with_info_2<int, K> Vbi2;
typedef CGAL::Triangulation_face_base_with_info_2<FaceInfo2, K> Tfbi2;
typedef CGAL::Constrained_triangulation_face_base_2<K, Tfbi2> Ctfb2;
typedef CGAL::Triangulation_data_structure_2<Vbi2, Ctfb2> Tds2;
typedef CGAL::Exact_predicates_tag Itag;
typedef CGAL::Constrained_Delaunay_triangulation_2<K, Tds2, Itag> CDT;

typedef CGAL::Delaunay_triangulation_2<K, Tds2> DT2D;
typedef K::Point_2 Point2;
typedef std::pair<Point2, int> IPoint2;

typedef CGAL::Projection_traits_xy_3<K> Pxy;
typedef CGAL::Triangulation_vertex_base_with_info_2<unsigned, Pxy> Vb_xy;
typedef CGAL::Triangulation_data_structure_2<Vb_xy> Tds_xy;
typedef CGAL::Delaunay_triangulation_2<Pxy, Tds_xy> DTXY;
typedef K::Point_3 Point3;

typedef std::pair<Point3, int> IPoint3;
typedef CGAL::Triangulation_vertex_base_with_info_3<int, K> Vb3;
typedef CGAL::Triangulation_data_structure_3<Vb3> Tds3;
typedef CGAL::Delaunay_triangulation_3<K, Tds3> DT3D;
typedef K::Tetrahedron_3 Tetrahedron3;
typedef Rcpp::NumericVector Dvector;