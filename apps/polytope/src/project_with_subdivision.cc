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


#include <cstdlib>
#include <polymake/Rational.h>
#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Vector.h>
#include <polymake/IncidenceMatrix.h>
#include <polymake/Graph.h>
#include <polymake/Array.h>
#include <polymake/Set.h>
#include <polymake/polytope/subdivisions.h>
#include <polymake/common/lattice_tools.h>

namespace polymake { namespace polytope {

    perl::Object project_with_subdivision(perl::Object p, int c) {

      int adim = p.give("CONE_AMBIENT_DIM");
      int dim = p.give("CONE_DIM");
      if ( dim != adim ) {
	throw std::runtime_error("project_with_subdivision: polytope must be full dimensional\n");       
      }

      Matrix<Rational> RSH(0,adim-1);	
      Matrix<Integer> SH(0,adim-1);	
      if ( p.lookup("SUBDIVISION") ) 
	p.lookup("SUBDIVISION.SUBDIVISION_HYPERPLANES") >> RSH;
      
      if ( RSH.rows() > 0 ) {
	SH.resize(RSH.rows(),RSH.cols());
	SH = common::primitive(RSH);
      }

      Matrix<Rational> RF = p.give("FACETS");
      Matrix<Integer> F = common::primitive(RF);
      Matrix<Rational> V = p.give("VERTICES");
      IncidenceMatrix<> VIF = p.give("VERTICES_IN_FACETS");
      Graph<> DG = p.give("DUAL_GRAPH.ADJACENCY");

      subdivision::SubdividedPolytope SP(F, V, SH, VIF, DG);
      subdivision::SubdividedPolytope SQ = subdivision::project_with_subdivision_step(SP,c);
      return SQ.polytope();
    }
    
    UserFunction4perl("# Category: Extension: projection with subdivision\n"
		      "# produces the push-forward subdivision of a polytope"
		      "# @param Polytope<Rational> the polytope"
		      "# @param c the projection coordinate (not that this is in the range 1..[[CONE_DIM]-1])",
		      &project_with_subdivision, 
		      "project_with_subdivision(Polytope,$)");
    
  } }

