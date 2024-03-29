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
#include "polymake/Array.h"
#include "polymake/Integer.h"
namespace polymake { namespace common {
///==== Automatically generated contents follow.    Please do not delete this line. ====
   template <typename T0>
   FunctionInterface4perl( new, T0 ) {
      WrapperReturnNew(T0, () );
   };

   template <typename T0, typename T1>
   FunctionInterface4perl( new_X, T0,T1 ) {
      perl::Value arg0(stack[1]);
      WrapperReturnNew(T0, (arg0.get<T1>()) );
   };

   Class4perl("Polymake::common::Array__Array__Integer", Array< Array< Integer > >);
   OperatorInstance4perl(assign, Array< int >, perl::Canned< const Array< Array< Integer > > >);
   FunctionInstance4perl(new_X, Array< Integer >, perl::Canned< const Array< Integer > >);
   FunctionInstance4perl(new, Array< Array< Integer > >);
   FunctionInstance4perl(new_X, Array< Array< Integer > >, perl::Canned< const Array< Array< Integer > > >);
   Class4perl("Polymake::common::Array__Integer", Array< Integer >);
///==== Automatically generated contents end here.  Please do not delete this line. ====
} }
