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


#ifndef _POLYMAKE_SUBDIVIDEDPOLYTOPE_H
#define _POLYMAKE_SUBDIVIDEDPOLYTOPE_H

#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/IncidenceMatrix.h>
#include <polymake/Rational.h>
#include <polymake/Integer.h>
#include <polymake/Set.h>
#include <polymake/Graph.h>
#include <polymake/Array.h>
#include <polymake/polytope/cdd_interface.h>
#include <polymake/graph/HasseDiagram.h>
#include <polymake/polytope/face_lattice_tools.h>
#include <polymake/linalg.h>
#include <polymake/polytope/subdivision_tools.h>

namespace polymake { namespace polytope { namespace subdivision {

    class SubdividedPolytope {

      Matrix<Integer> _facets;
      Matrix<Rational> _vertices;
      Matrix<Integer> _subdivision_hyperplanes;
      IncidenceMatrix<> _vif;
      Graph<> _dg;
      int _dim;
      int _number;
      
    public: 

      const Matrix<Rational> vertices() const { return _vertices; }
      const Matrix<Integer> facets() const { return _facets; }
      const Matrix<Integer> subdivision_hyperplanes() const { return _subdivision_hyperplanes; }
      const IncidenceMatrix<> vertices_in_facets() const { return _vif; }
      const Graph<> dual_graph() const { return _dg; }
      const int dim() const { return _dim; }
      void set_number(int n) { _number = n; }
      const int number() const { return _number; }
      
      Set<int> valid_projection_coordinates() {
       
	Matrix<Integer> M = _facets/_subdivision_hyperplanes;
	Set<int> VPC;
	for ( int i = 1; i < M.cols(); ++i ) {
	  bool valid = true;
	  for ( int j = 0; j < M.rows() && valid; ++j )
	    if ( abs(M(j,i) ) > 1 )
	      valid = false;
	  if ( valid ) 
	    VPC += i;
	}
	
	return VPC;
      }


      const perl::Object polytope() const {
	perl::Object p("LatticePolytope");

	p.take("FACETS") << _facets;
	p.take("VERTICES") << _vertices;
	p.take("VERTICES_IN_FACETS") << _vif;
	p.take("CONE_DIM") << _dim+1;
	p.take("DUAL_GRAPH.ADJACENCY") << _dg;
	p.take("SUBDIVISION.SUBDIVISION_HYPERPLANES") << _subdivision_hyperplanes;
	return p;
      }

      SubdividedPolytope ( const SubdividedPolytope & SP ) {
	_facets = SP.facets();
	_vertices = SP.vertices();
	_subdivision_hyperplanes = SP.subdivision_hyperplanes();
	_number = SP.number();
	_vif = SP.vertices_in_facets();
	_dg = SP.dual_graph();
	_dim = SP.dim();
      }

      SubdividedPolytope(const Matrix<Integer> & F, const Matrix<Integer> & S) {
	
	_facets = F;
	_subdivision_hyperplanes = S;
	Matrix<Rational> AH(0,F.cols());
	cdd_interface::solver<Rational> solver;
	_vertices = solver.enumerate_vertices(Matrix<Rational>(_facets),AH).first;
	_dim = F.cols()-1;
	_vif = construct_vif(_facets, _vertices);
	_dg = construct_dual_graph(_vif);
      }
      
      SubdividedPolytope(const Matrix<Rational> & F, const Matrix<Rational> & S) {
	
	_facets = primitive(multiply_by_common_denominator(F));
	_subdivision_hyperplanes = primitive(multiply_by_common_denominator(S));
	Matrix<Rational> AH(0,F.cols());
	cdd_interface::solver<Rational> solver;
	_vertices = solver.enumerate_vertices(Matrix<Rational>(_facets),AH).first;
	_dim = F.cols()-1;
	_vif = construct_vif(_facets, _vertices);
	_dg = construct_dual_graph(_vif);
      }

      SubdividedPolytope(const Matrix<Integer> & F, const Matrix<Rational>& V, const Matrix<Integer> & S) {
	_facets = F;
	_subdivision_hyperplanes = S;
	_vertices = V;
	_dim = F.cols()-1;
	IncidenceMatrix<> vifcheck = construct_vif(_facets, _vertices);
	perl::Object q("Polytope<Rational>");
	q.take("VERTICES") << V;
	q.take("FACETS") << F;
	IncidenceMatrix<> VIF2 = q.give("VERTICES_IN_FACETS");
	_vif = VIF2;
	Graph<> DG = q.give("DUAL_GRAPH.ADJACENCY");
	_dg = DG;
      }

      SubdividedPolytope(const Matrix<Rational> & F, const Matrix<Rational>& V, const Matrix<Rational> & S, 
			 const IncidenceMatrix<>& VIF, const Graph<> DG) {
	
	_facets = primitive(multiply_by_common_denominator(F));
	_subdivision_hyperplanes = primitive(multiply_by_common_denominator(Matrix<Rational>(S)));
	_vertices = V;
	_dim = F.cols()-1;
	_vif = VIF;
	_dg = DG;
      }

    }; // end of class SubdividedPolytope

    } } }

#endif
