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
#include <polymake/Set.h>
#include <polymake/linalg.h>



namespace polymake { namespace polytope {

  perl::Object scale_to_lattice(perl::Object p) {

    Matrix<Rational> M = p.give("VERTICES");

    Integer d = 1;

    for ( int i = 0; i < M.rows(); ++i )
      for ( int j = 1; j < M.cols(); ++j ) 
	d = lcm ( d, denominator(M.row(i)[j]) );

    M = M*d;
    M.col(0) = ones_vector<Rational>(M.rows());

    perl::Object po("Polytope<Rational>");

    po.take("VERTICES") << M;
    return po;
  }

  UserFunction4perl("# Category: Other\n"
		    "# ",
		    &scale_to_lattice, 
		    "scale_to_lattice(Polytope<Rational>)");

} }

