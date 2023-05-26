#ifndef _DELAUNAYHEADER_
#include "delaunay.h"
#endif

// [[Rcpp::export]]
Rcpp::List del2D_cpp(Rcpp::NumericMatrix pts) {
  const int npoints = pts.ncol();
  std::vector<IPoint2> points(npoints);
  Mesh mesh;
  for(int i = 0; i < npoints; i++) {
    Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    Point2 pt(pt_i(0), pt_i(1));
    points[i] = std::make_pair(pt, i);
    mesh.add_vertex(pt);
  }
  // compute Delaunay mesh
  const DT2D dmesh(points.begin(), points.end());
  //
  const size_t nfaces = dmesh.number_of_faces();
  const size_t h =
      2 * npoints - 2 - nfaces;  // number of vertices of convex hull
  const size_t nedges = 3 * npoints - 3 - h;
  //
  Rcpp::IntegerMatrix Faces(3, nfaces);
  {
    int i = 0;
    for(DT2D::Finite_faces_iterator fit = dmesh.finite_faces_begin();
                                    fit != dmesh.finite_faces_end(); fit++) {
      const int id0 = fit->vertex(0)->info();
      const int id1 = fit->vertex(1)->info();
      const int id2 = fit->vertex(2)->info();
      Faces(Rcpp::_, i++) = Rcpp::IntegerVector::create(id0+1, id1+1, id2+1);
      const face_descriptor fd = mesh.add_face(
        vertex_descriptor(id0), 
        vertex_descriptor(id1), 
        vertex_descriptor(id2)
      );
      if(fd == Mesh::null_face()) { // that should not happen
        Rcpp::stop("The face could not be added.");
      }
    }
  }
  //
  const DT2D::Finite_edges itedges = dmesh.finite_edges();
  Rcpp::IntegerMatrix Edges(2, nedges);
  {
    int i = 0;
    for(DT2D::Finite_edges_iterator eit = itedges.begin(); 
                                    eit != itedges.end(); eit++) {
      const std::pair<DT2D::Face_handle, int> edge = *eit;
      const int i0 = edge.first->vertex((edge.second + 1) % 3)->info();
      const int i1 = edge.first->vertex((edge.second + 2) % 3)->info();
      Edges(Rcpp::_, i++) = Rcpp::IntegerVector::create(i0 + 1, i1 + 1);
    }
  }
  //
  return Rcpp::List::create(
    Rcpp::Named("faces") = Faces,
    Rcpp::Named("edges") = Edges,
    Rcpp::Named("mesh")  = Rcpp::List::create(
      Rcpp::Named("vertices") = getVertices(mesh),
      Rcpp::Named("edges")    = getEdges(mesh),
      Rcpp::Named("faces")    = getFacesInfo(mesh) 
    )
  );
}