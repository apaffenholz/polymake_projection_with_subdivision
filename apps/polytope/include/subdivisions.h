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


#ifndef _POLYMAKE_SUBDIVISIONS_H
#define _POLYMAKE_SUBDIVISIONS_H

#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Vector.h>
#include <polymake/IncidenceMatrix.h>
#include <polymake/Rational.h>
#include <polymake/Integer.h>
#include <polymake/Set.h>
#include <polymake/Graph.h>
#include <polymake/Array.h>
#include <polymake/linalg.h>
#include <polymake/polytope/SubdividedPolytope.h>
#include <polymake/polytope/SubdivisionCell.h>
#include <polymake/polytope/subdivision_tools.h>

namespace polymake { namespace polytope { namespace subdivision {
    
      // checks whether a subdivision hyperplane has points of a point set on both sides
      bool sub_intersects ( const Vector<Integer> & H, const Matrix<Rational> & V );

      // computes the cells of a subdivision
      Vector<SubdivisionCell> compute_cells ( const Matrix<Rational> & V, 
					      const Matrix<Integer> & F, 
					      const Matrix<Integer> & S );
     

      // do a projection step of a SubdividedPolytope in direction c
      SubdividedPolytope 
	project_with_subdivision_step(SubdividedPolytope & SP, int c); 
      
    }
    
    // computes the cells of a subdivision and matches the vertices with the lattice points of the polytope
    Array< Set<int> > subdivision_cells (  const Matrix<Rational> & V, 
					   const Matrix<Integer> & F, 
					   const Matrix<Integer> & S,
					   const Matrix<Integer> & L );
    
    // computes the cells of a subdivision 
    // and checks whether all cells have integral vertices
    bool lattice_subdivision(const Matrix<Rational> & V, 
			     const Matrix<Integer> & F, 
			     const Matrix<Integer> & S ); 
    
    
    
  } }

#endif // _POLYMAKE_SUBDIVISIONS_H
