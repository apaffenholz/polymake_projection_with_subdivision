
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


#ifndef _POLYMAKE_SUBDIVISION_TOOLS_H
#define _POLYMAKE_SUBDIVISION_TOOLS_H

#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Vector.h>
#include <polymake/IncidenceMatrix.h>
#include <polymake/Rational.h>
#include <polymake/Integer.h>
#include <polymake/Set.h>
#include <polymake/Graph.h>
#include <polymake/graph/HasseDiagram.h>
#include <polymake/polytope/face_lattice_tools.h>
#include <polymake/linalg.h>

namespace polymake { namespace polytope { namespace subdivision {
      
      IncidenceMatrix<> construct_vif (const Matrix<Integer>& F, const Matrix<Rational>& V);
	
	template <typename Iterator>
	void fill_dual_graph(Graph<>& G, Iterator rif_it, int lowest_index);
      
      Graph<> construct_dual_graph (const IncidenceMatrix<> VIF); 

      bool hyperplane_intersects ( const Vector<Integer> & h, const Matrix<Rational> & V );

      Matrix<Integer> make_equations_unique(const Matrix<Integer> & M );
      
    } } }

#endif
