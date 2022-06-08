#ifndef _DELAUNAYHEADER_
#include "delaunay.h"
#endif

// [[Rcpp::export]]
Rcpp::List del2D_cpp(Rcpp::NumericMatrix pts) {
  const int npoints = pts.ncol();
  std::vector<IPoint2> points(npoints);
  for(int i = 0; i < npoints; i++) {
    Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    points[i] = std::make_pair(Point2(pt_i(0), pt_i(1)), i + 1);
  }
  // compute Delaunay mesh
  const DT2D mesh(points.begin(), points.end());
  //
  const size_t nfaces = mesh.number_of_faces();
  const size_t h =
    2 * npoints - 2 - nfaces;  // number of vertices of convex hull
  const size_t nedges = 3 * npoints - 3 - h;
  //
  Rcpp::IntegerMatrix Faces(3, nfaces);
  {
    int i = 0;
    for(DT2D::Finite_faces_iterator fit = mesh.finite_faces_begin();
        fit != mesh.finite_faces_end(); fit++) {
      const int id0 = fit->vertex(0)->info();
      const int id1 = fit->vertex(1)->info();
      const int id2 = fit->vertex(2)->info();
      Faces(Rcpp::_, i) = Rcpp::IntegerVector::create(id0, id1, id2);
      i++;
    }
  }
  //
  const DT2D::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix Edges(2, nedges);
  {
    int i = 0;
    for(DT2D::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
    eit++) {
      const std::pair<DT2D::Face_handle, int> edge = *eit;
      const int i0 = edge.first->vertex((edge.second + 1) % 3)->info();
      const int i1 = edge.first->vertex((edge.second + 2) % 3)->info();
      Edges(Rcpp::_, i) = Rcpp::IntegerVector::create(i0, i1);
      i++;
    }
  }
  //
  return Rcpp::List::create(Rcpp::Named("faces") = Faces,
                            Rcpp::Named("edges") = Edges);
}

