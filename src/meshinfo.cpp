#ifndef _DELAUNAYHEADER_
#include "delaunay.h"
#endif


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
Rcpp::DataFrame getEdges(Mesh& mesh) {
  const size_t nedges = mesh.number_of_edges();
  Rcpp::IntegerVector I1(nedges);
  Rcpp::IntegerVector I2(nedges);
  Rcpp::NumericVector Length(nedges);
  Rcpp::LogicalVector Border(nedges);
  Rcpp::IntegerVector F1(nedges);
  Rcpp::IntegerVector F2(nedges);
  {
    size_t i = 0;
    for(edge_descriptor ed : mesh.edges()) {
      const vertex_descriptor s = source(ed, mesh);
      const vertex_descriptor t = target(ed, mesh);
      I1(i) = (int)s + 1;
      I2(i) = (int)t + 1;
      const halfedge_descriptor h0 = mesh.halfedge(ed, 0);
      Length(i) = PMP::edge_length(h0, mesh);
      const bool isBorder = mesh.is_border(ed);
      Border(i) = isBorder;
      F1(i) = int(mesh.face(h0)) + 1;
      if(isBorder) {
        F2(i) = Rcpp::IntegerVector::get_na();
      } else {
        const halfedge_descriptor h1 = mesh.halfedge(ed, 1);
        F2(i) = int(mesh.face(h1)) + 1;
      }
      i++;
    }
  }
  Rcpp::DataFrame Edges = Rcpp::DataFrame::create(
    Rcpp::Named("i1")       = I1,
    Rcpp::Named("i2")       = I2,
    Rcpp::Named("length")   = Length,
    Rcpp::Named("border")   = Border,
    Rcpp::Named("f1")       = F1,
    Rcpp::Named("f2")       = F2
  );
  return Edges;
}

// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
Rcpp::NumericMatrix getFacesInfo(Mesh& mesh) {
  const Rcpp::CharacterVector rownames = {"ccx", "ccy", "area"};
  Rcpp::NumericMatrix FacesInfo(3, mesh.number_of_faces());
  int i = 0;
  for(face_descriptor fd : mesh.faces()) {
    auto vs = vertices_around_face(mesh.halfedge(fd), mesh).begin();
    const Point2 p1 = mesh.point(*(vs++));
    const Point2 p2 = mesh.point(*(vs++));
    const Point2 p3 = mesh.point(*vs);
    const Point2 circumcenter = CGAL::circumcenter(p1, p2, p3);
    const double area         = CGAL::area(p1, p2, p3);
    Rcpp::NumericVector col_i = {
      circumcenter.x(),
      circumcenter.y(),
      area
    };
    FacesInfo(Rcpp::_, i++) = col_i;
  }
  Rcpp::rownames(FacesInfo) = rownames;
  return Rcpp::transpose(FacesInfo);
}


// -------------------------------------------------------------------------- //
// -------------------------------------------------------------------------- //
Rcpp::NumericMatrix getVertices(Mesh& mesh) {
  const size_t nvertices = mesh.number_of_vertices();
  Rcpp::NumericMatrix Vertices(2, nvertices);
  {
    size_t i = 0;
    for(vertex_descriptor vd : mesh.vertices()) {
      Rcpp::NumericVector col_i(2);
      const Point2 vertex = mesh.point(vd);
      col_i(0) = vertex.x();
      col_i(1) = vertex.y();
      Vertices(Rcpp::_, i++) = col_i;
    }
  }
  return Rcpp::transpose(Vertices);
}
