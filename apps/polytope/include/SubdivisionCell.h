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



#ifndef _POLYMAKE_SUBDIVISIONCELL_H
#define _POLYMAKE_SUBDIVISIONCELL_H

#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Rational.h>
#include <polymake/Integer.h>
#include <polymake/Set.h>
#include <polymake/Graph.h>
#include <polymake/polytope/cdd_interface.h>
#include <polymake/linalg.h>
#include <polymake/polytope/subdivision_tools.h>
#include <polymake/common/lattice_tools.h>

namespace polymake { namespace polytope { namespace subdivision {

    class SubdivisionCell {
      
      Matrix<Rational> _vertices;
      Matrix<Integer> _facets;

    public: 

      const Matrix<Rational> vertices() const { return _vertices; }
      const Matrix<Integer> facets() const { return _facets; }
      

      SubdivisionCell ( ) { }

      SubdivisionCell ( const Matrix<Integer> & I ) { 
	cdd_interface::ConvexHullSolver<Rational> solver;
	Matrix<Rational> AH(0,I.cols());
	Matrix<Rational> Lin(0,I.cols());
	_vertices = solver.enumerate_vertices(Matrix<Rational>(I),AH,0).first;
	auto F=solver.enumerate_facets(_vertices,Lin,0);
	_facets = common::primitive(F.first); 
      }
      
      SubdivisionCell ( const SubdivisionCell & SC ) {
	_vertices = SC.vertices();
	_facets = SC.facets();
      }

      ~SubdivisionCell () {}

      const bool integral () const {
	for (auto x=entire(concat_rows(_vertices)); !x.at_end(); ++x)
	  if ( denominator(*x)!=1 ) return false;
	return true;
      }

    }; // end of class SubdivisionCell    


    } } }
#endif
