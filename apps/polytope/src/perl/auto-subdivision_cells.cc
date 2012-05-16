/* Copyright (c) 1997-2010
   Ewgenij Gawrilow, Michael Joswig (Technische Universitaet Darmstadt, Germany)
   http://www.polymake.de

   This program is free software; you can redistribute it and/or modify it
   under the terms of the GNU General Public License as published by the
   Free Software Foundation; either version 2, or (at your option) any
   later version: http://www.gnu.org/licenses/gpl.txt.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.
*/

///==== this line controls the automatic file splitting: max.instances=40

#include "polymake/client.h"
#include "polymake/polytope/subdivisions.h"
#include "polymake/Matrix.h"
#include "polymake/Rational.h"
#include "polymake/IncidenceMatrix.h"
#include "polymake/Integer.h"
namespace polymake { namespace polytope {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0, typename T1, typename T2, typename T3>
   FunctionInterface4perl( subdivision_cells_X_X_X_X, T0,T1,T2,T3 ) {
      perl::Value arg0(stack[0]), arg1(stack[1]), arg2(stack[2]), arg3(stack[3]);
      WrapperReturn( subdivision_cells(arg0.get<T0>(), arg1.get<T1>(), arg2.get<T2>(), arg3.get<T3>()) );
   };

   FunctionInstance4perl(subdivision_cells_X_X_X_X, perl::Canned< const Matrix< Rational > >, perl::Canned< const Matrix< Integer > >, perl::Canned< const Matrix< Integer > >, perl::Canned< const Matrix< Integer > >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
