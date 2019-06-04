/*
  Copyright (c) 2010-12 Andreas Paffenholz
 
  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the
  Free Software Foundation; either version 2, or (at your option) any
  later version: http://www.gnu.org/licenses/gpl.txt.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.
*/

#include <polymake/polytope/subdivision_tools.h>
#include <polymake/polytope/face_lattice_tools.h>

namespace polymake { namespace polytope { namespace subdivision {
      
      template <typename Iterator>
      void fill_dual_graph(Graph<>& G, Iterator rif_it, int lowest_index) {
	for (; !rif_it.at_end(); ++rif_it) 
	  G.edge(rif_it->front() - lowest_index, rif_it->back() - lowest_index);
      }
      
      IncidenceMatrix<> construct_vif (const Matrix<Integer>& F, const Matrix<Rational>& V) {
	return IncidenceMatrix<> (F.rows(), V.rows(),
				  attach_operation(product(rows(F), rows(V), operations::mul()),
						   operations::is_zero()).begin());
      }
      
      Graph<> construct_dual_graph (const IncidenceMatrix<> VIF) {
	
	perl::Object q("Polytope<Rational>");
	q.take("VERTICES_IN_FACETS") << VIF;
	graph::Lattice<graph::lattice::BasicDecoration, graph::lattice::Sequential> HD = q.give("HASSE_DIAGRAM");
	const graph::Lattice<graph::lattice::BasicDecoration, graph::lattice::Sequential>::nodes_of_rank_type facet_nodes=HD.nodes_of_rank(HD.rank()-1);
	Graph<> G(facet_nodes.size());
	
	fill_dual_graph(G, 
			entire(select(rows(adjacency_matrix(HD.graph())),HD.nodes_of_rank(HD.rank()-2))), 
			facet_nodes.front() );
	
	return G;
      }
      
      // checks whether the affine hyperplane h has points of a point set V on both sides
      bool hyperplane_intersects ( const Vector<Integer> & h, const Matrix<Rational> & V ) {
	
	int i = 0;
	while ( i < V.rows() && h * V.row(i) == 0 ) 
	  ++i;
	
	if ( i == V.rows() )
	  return false;
	
	if ( h * V.row(i) < 0 )
	  while ( i < V.rows() && h * V.row(i) <= 0 ) 
	    ++i;
	else
	  while ( i < V.rows() && h * V.row(i) >= 0 ) 
	    ++i;
	
	if ( i == V.rows() )
	  return false;
	
	return true;
      }

      // remove duplicate rows from a matrix
      // FIXME not very efficient
      Matrix<Integer> make_equations_unique(const Matrix<Integer> & M ) {
	Matrix<Integer> N(0,M.cols());
	for ( auto mit = entire(rows(M)); !mit.at_end(); ++mit ) {
	  bool known = false;
	  for ( auto nit = entire(rows(N)); !nit.at_end() && !known; ++nit )
	    if ( (*nit) == (*mit) || (*nit) == -(*mit) )
	      known = true;
	  if ( !known )
	    N /= *mit;
	}
	
	return N;
      }

    
  } } }

