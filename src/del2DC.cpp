#ifndef _DELAUNAYHEADER_
#include "delaunay.h"
#endif

void mark_domains0(CDT& ct,
                   CDT::Face_handle start,
                   int index,
                   std::list<CDT::Edge>& border) {
  if(start->info().nesting_level != -1) {
    return;
  }
  std::list<CDT::Face_handle> queue;
  queue.push_back(start);
  while(!queue.empty()) {
    CDT::Face_handle fh = queue.front();
    queue.pop_front();
    if(fh->info().nesting_level == -1) {
      fh->info().nesting_level = index;
      for(int i = 0; i < 3; i++) {
        CDT::Edge e(fh, i);
        CDT::Face_handle n = fh->neighbor(i);
        if(n->info().nesting_level == -1) {
          if(ct.is_constrained(e))
            border.push_back(e);
          else
            queue.push_back(n);
        }
      }
    }
  }
}

void mark_domains(CDT& cdt) {
  for(CDT::Face_handle f : cdt.all_face_handles()) {
    f->info().nesting_level = -1;
  }
  std::list<CDT::Edge> border;
  mark_domains0(cdt, cdt.infinite_face(), 0, border);
  while(!border.empty()) {
    CDT::Edge e = border.front();
    border.pop_front();
    CDT::Face_handle n = e.first->neighbor(e.second);
    if(n->info().nesting_level == -1) {
      mark_domains0(cdt, n, e.first->info().nesting_level + 1, border);
    }
  }
}

// [[Rcpp::export]]
Rcpp::List del2DC_cpp(Rcpp::NumericMatrix pts,
                               Rcpp::IntegerMatrix edges) {
  const int npoints = pts.ncol();
  std::vector<std::pair<CDT::Point, int>> points(npoints);
  Mesh mesh;
  for(int i = 0; i < npoints; ++i) {
    const Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    points[i] = std::make_pair(CDT::Point(pt_i(0), pt_i(1)), i);
    mesh.add_vertex(Point2(pt_i(0), pt_i(1)));
  }
  CDT cdt;
  {
    const size_t nedges = edges.ncol();
    for(size_t k = 0; k < nedges; ++k) {
      const Rcpp::IntegerVector edge_k = edges(Rcpp::_, k);
      cdt.insert_constraint(points[edge_k(0) - 1].first,
                            points[edge_k(1) - 1].first);
    }
  }
  cdt.insert(points.begin(), points.end());
  const size_t nfaces = cdt.number_of_faces();
  Rcpp::IntegerMatrix faces(3, nfaces);
  mark_domains(cdt);
  size_t nfaces_out;
  {
    size_t i = 0;
    for(CDT::Face_handle f : cdt.finite_face_handles()){
      if(f->info().in_domain()){
        const int id0 = f->vertex(0)->info();
        const int id1 = f->vertex(1)->info();
        const int id2 = f->vertex(2)->info();
        faces(Rcpp::_, i++) = Rcpp::IntegerVector::create(id0+1, id1+1, id2+1);
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
    nfaces_out = i;
  }
  //
  return Rcpp::List::create(
    Rcpp::Named("faces") = faces(Rcpp::_, Rcpp::Range(0, nfaces_out-1)), 
    Rcpp::Named("mesh")  = Rcpp::List::create(
      Rcpp::Named("vertices") = getVertices(mesh),
      Rcpp::Named("edges")    = getEdges(mesh),
      Rcpp::Named("faces")    = getFacesInfo(mesh) 
    )
  );
}
