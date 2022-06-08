#ifndef _DELAUNAYHEADER_
#include "delaunay.h"
#endif

double volume_under_triangle(Dvector v0, Dvector v1, Dvector v2) {
  double x0 = v0(0);
  double y0 = v0(1);
  double z0 = v0(2);
  double x1 = v1(0);
  double y1 = v1(1);
  double z1 = v1(2);
  double x2 = v2(0);
  double y2 = v2(1);
  double z2 = v2(2);
  return (z0 + z1 + z2) *
         (x0 * y1 - x1 * y0 + x1 * y2 - x2 * y1 + x2 * y0 - x0 * y2) / 6.0;
}

// [[Rcpp::export]]
Rcpp::List delXY_cpp(Rcpp::NumericMatrix pts) {
  const int npoints = pts.ncol();
  // compute Delaunay mesh
  DTXY mesh;
  DTXY::Vertex_handle vh;
  for(int i = 0; i < npoints; i++) {
    Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    vh = mesh.insert(Point3(pt_i(0), pt_i(1), pt_i(2)));
    vh->info() = i + 1;
  }
  //
  const size_t nfaces = mesh.number_of_faces();
  const size_t h =
      2 * npoints - 2 - nfaces;  // number of vertices of convex hull
  const size_t nedges = 3 * npoints - 3 - h;
  //
  Rcpp::IntegerMatrix Faces(3, nfaces);
  Rcpp::NumericVector Volumes(nfaces);
  double totalVolume = 0.0;
  {
    size_t i = 0;
    for(DTXY::Finite_faces_iterator fit = mesh.finite_faces_begin();
        fit != mesh.finite_faces_end(); fit++) {
      const int i0 = fit->vertex(0)->info();
      const int i1 = fit->vertex(1)->info();
      const int i2 = fit->vertex(2)->info();
      Faces(Rcpp::_, i) = Rcpp::IntegerVector::create(i0, i1, i2);
      double volume = volume_under_triangle(
          pts(Rcpp::_, i0 - 1), pts(Rcpp::_, i1 - 1), pts(Rcpp::_, i2 - 1));
      Volumes(i) = volume;
      totalVolume += volume;
      i++;
    }
    Faces.attr("volumes") = Volumes;
  }
  //
  const DTXY::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix Edges(2, nedges);
  {
    size_t i = 0;
    for(DTXY::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
        eit++) {
      const std::pair<DTXY::Face_handle, int> edge = *eit;
      const int i0 = edge.first->vertex((edge.second + 1) % 3)->info();
      const int i1 = edge.first->vertex((edge.second + 2) % 3)->info();
      Edges(Rcpp::_, i) = Rcpp::IntegerVector::create(i0, i1);
      i++;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("faces") = Faces,
                            Rcpp::Named("edges") = Edges,
                            Rcpp::Named("volume") = totalVolume);
}