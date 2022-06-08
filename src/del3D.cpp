#ifndef _DELAUNAYHEADER_
#include "delaunay.h"
#endif

Rcpp::String stringPair(const int i, const int j) {
  const int i0 = std::min(i, j);
  const int i1 = std::max(i, j);
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
      std::to_string(i0), "-", std::to_string(i1));
  return Rcpp::collapse(stringids);
}

std::array<Rcpp::String, 3> triangleEdges(const int i0,
                                          const int i1,
                                          const int i2) {
  return {stringPair(i0, i1), stringPair(i0, i2), stringPair(i1, i2)};
}

Rcpp::String stringTriple(const int i, const int j, const int k) {
  const Rcpp::CharacterVector stringids = Rcpp::CharacterVector::create(
      std::to_string(i), "-", std::to_string(j), "-", std::to_string(k));
  return Rcpp::collapse(stringids);
}

// [[Rcpp::export]]
Rcpp::List del3D_cpp(Rcpp::NumericMatrix pts) {
  const int npoints = pts.ncol();
  std::vector<IPoint3> points(npoints);
  for(int i = 0; i < npoints; i++) {
    Rcpp::NumericVector pt_i = pts(Rcpp::_, i);
    points[i] = std::make_pair(Point3(pt_i(0), pt_i(1), pt_i(2)), i + 1);
  }
  // compute Delaunay mesh
  const DT3D mesh(points.begin(), points.end());
  //
  const size_t nfacets = mesh.number_of_finite_facets();
  const size_t ncells = mesh.number_of_finite_cells();
  const size_t nedges = mesh.number_of_finite_edges();
  //
  std::map<Rcpp::String, size_t> facetsMap = {};
  Rcpp::IntegerMatrix Facets(4, nfacets);
  std::vector<Rcpp::String> edgesOnHull(0);
  {
    size_t i = 0;
    for(DT3D::Finite_facets_iterator fit = mesh.finite_facets_begin();
        fit != mesh.finite_facets_end(); fit++) {
      std::pair<DT3D::Cell_handle, int> facet = *fit;
      DT3D::Vertex_handle v0 = facet.first->vertex((facet.second + 1) % 4);
      DT3D::Vertex_handle v1 = facet.first->vertex((facet.second + 2) % 4);
      DT3D::Vertex_handle v2 = facet.first->vertex((facet.second + 3) % 4);
      const int id0 = v0->info();
      const int id1 = v1->info();
      const int id2 = v2->info();
      Facets(Rcpp::_, i) =
          Rcpp::IntegerVector::create(id0, id1, id2, 0);
      bool onhull = mesh.is_infinite(facet.first) ||
                    mesh.is_infinite(mesh.mirror_facet(facet).first);
      if(onhull) {
        Facets(3, i) = onhull;
        std::array<Rcpp::String, 3> triangle = triangleEdges(id0, id1, id2);
        edgesOnHull.insert(edgesOnHull.end(), triangle.begin(), triangle.end());
      }
      std::array<int, 3> ids = {id0, id1, id2};
      std::sort(ids.begin(), ids.end());
      const Rcpp::String facetAsString = stringTriple(ids[0], ids[1], ids[2]);
      i++;
      facetsMap[facetAsString] = i;
    }
    Rcpp::CharacterVector rowNames =
        Rcpp::CharacterVector::create("i1", "i2", "i3", "onhull");
    Rcpp::rownames(Facets) = rowNames;
    Facets.attr("info") =
        "The `onhull` column indicates whether the face is on the convex hull.";
  }

  const DT3D::Finite_edges itedges = mesh.finite_edges();
  Rcpp::IntegerMatrix Edges(3, nedges);
  {
    size_t i = 0;
    for(DT3D::Finite_edges_iterator eit = itedges.begin(); eit != itedges.end();
        eit++) {
      const CGAL::Triple<DT3D::Cell_handle, int, int> edge = *eit;
      const int i0 = edge.first->vertex(edge.second)->info();
      const int i1 = edge.first->vertex(edge.third)->info();
      Edges(Rcpp::_, i) =
          Rcpp::IntegerVector::create(i0, i1, 0);
      const Rcpp::String i0i1 = stringPair(i0, i1);
      if(std::find(edgesOnHull.begin(), edgesOnHull.end(), i0i1) !=
         edgesOnHull.end()) {
        Edges(2, i) = 1;
      }
      i++;
    }
    Rcpp::CharacterVector rowNames =
        Rcpp::CharacterVector::create("i1", "i2", "onhull");
    Rcpp::rownames(Edges) = rowNames;
    Edges.attr("info") =
        "The `onhull` column indicates whether the edge is on the convex hull.";
  }

  Rcpp::List Cells(ncells);
  double totalVolume = 0.0;
  {
    size_t i = 0;
    for(DT3D::Finite_cells_iterator cit = mesh.finite_cells_begin();
        cit != mesh.finite_cells_end(); cit++) {
      const int id0 = cit->vertex(0)->info();
      const int id1 = cit->vertex(1)->info();
      const int id2 = cit->vertex(2)->info();
      const int id3 = cit->vertex(3)->info();
      Rcpp::IntegerVector Cell =
          Rcpp::IntegerVector::create(id0, id1, id2, id3);
      std::array<int, 4> ids = {id0, id1, id2, id3};
      std::sort(ids.begin(), ids.end());
      const Rcpp::IntegerVector Faces = Rcpp::IntegerVector::create(
          facetsMap[stringTriple(ids[0], ids[1], ids[2])],
          facetsMap[stringTriple(ids[0], ids[1], ids[3])],
          facetsMap[stringTriple(ids[0], ids[2], ids[3])],
          facetsMap[stringTriple(ids[1], ids[2], ids[3])]);
      Tetrahedron3 th = CGAL::Tetrahedron_3<K>(
          cit->vertex(0)->point(), cit->vertex(1)->point(),
          cit->vertex(2)->point(), cit->vertex(3)->point());
      const double volume = th.volume();
      totalVolume += volume;
      Cells[i] = Rcpp::List::create(Rcpp::Named("cell") = Cell,
                                    Rcpp::Named("faces") = Faces,
                                    Rcpp::Named("volume") = volume);
      i++;
    }
  }
  //
  return Rcpp::List::create(
      Rcpp::Named("cells") = Cells, Rcpp::Named("facets") = Facets,
      Rcpp::Named("edges") = Edges, Rcpp::Named("volume") = totalVolume);
}