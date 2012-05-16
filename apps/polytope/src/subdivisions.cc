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


#include <polymake/polytope/subdivisions.h>

namespace polymake { namespace polytope { namespace subdivision {
    
      // checks whether a subdivision hyperplane has points of a point set on both sides
      bool sub_intersects ( const Vector<Integer> & H, const Matrix<Rational> & V ) {
	
	int i = 0;
	while ( H * V.row(i) == 0 && i < V.rows() ) i++;
	if ( i < V.rows() ) {
	  int a = sign( H * V.row(i++) );
	  for ( int j = i; j < V.rows(); ++j )
	    if ( sign( H *V.row(j) ) == -a ) return true;
	}
	
	return false;
      } // end of sub_intersects
      
      // computes the cells of a subdivision
      Vector<SubdivisionCell> compute_cells ( const Matrix<Rational> & V, 
					      const Matrix<Integer> & F, 
					      const Matrix<Integer> & S ) {
	
	Vector<SubdivisionCell> cells;
	
	// we start with the whole polytope
	SubdivisionCell initial(F);
	cells |= initial;
	
	// and successively intersect all cells with a subdivision hyperplane
	for ( Entire<Rows<Matrix<Integer> > >::const_iterator sub = entire(rows(S)); !sub.at_end(); ++sub ) {
	  
	  Vector<SubdivisionCell> cells2;
	  for ( Entire<Vector<SubdivisionCell> >::const_iterator c = entire(cells); !c.at_end(); ++c ) {
	    
	    if ( sub_intersects(*sub,(*c).vertices()) ) {
	      SubdivisionCell c1((*c).facets()/(*sub));
	      SubdivisionCell c2((*c).facets()/(-(*sub)));
	      cells2 |= c1;
	      cells2 |= c2;
	    } else {
	      cells2 |= *c;
	    }
	  } // end of loop over cells
	  cells.resize(cells2.dim());
	  cells = cells2;
	} // end of loop over subdivision hyperplanes
	
	
	return cells;
      }  // end of compute_cells

      // do a projection step of a SubdividedPolytope in direction c
    SubdividedPolytope 
    project_with_subdivision_step(SubdividedPolytope & SP, int c) {
	
	// we need these quite often
	Matrix<Integer> PF = SP.facets();
	Matrix<Integer> PSubdivision = SP.subdivision_hyperplanes();

	// check that we are allowed to project in direction c 
	for ( Entire<Rows<Matrix<Integer> > >::const_iterator fit = entire(rows(PF)); !fit.at_end(); ++fit ) 
	  if ( abs((*fit)[c]) > 1 )
	    throw std::runtime_error("project_with_subdivision: not a valid projection coordinate\n");       
	
	for ( Entire<Rows<Matrix<Integer> > >::const_iterator sdit = entire(rows(PSubdivision)); !sdit.at_end(); ++sdit ) 
	  if ( abs((*sdit)[c]) > 1 )
	    throw std::runtime_error("project_with_subdivision: not a valid projection coordinate\n");       

	// compute vertices and facets of the projection
	// FIXME should be done without the detour to perl
	perl::Object q("LatticePolytope");
	q.take("POINTS") << SP.vertices().minor(All,~range(c,c));      

	Matrix<Rational> QV = q.give("VERTICES");
	Matrix<Integer> QF = primitive(multiply_by_common_denominator(q.give("FACETS")));

	Matrix<Integer>  QSubdivision(0,QF.cols());  // the new subdivision hyperplanes

	// if there is a subdivision, then subdivision hyperplanes orthogonal 
	// to the projection direction can immediately be projected.
	for ( int i = 0; i < PSubdivision.rows(); ++i ) 
	  if ( PSubdivision(i,c) == 0 ) {
	    Vector<Integer> h = PSubdivision.row(i).slice(~range(c,c));
	    if ( hyperplane_intersects(h,QV) ) // check whether the projected hyperplane intersects the projection
	      QSubdivision /= h;
	  }

	// combine subdivision with subdivision and facets
	for ( int i = 0; i < PSubdivision.rows(); ++i )        // loop over subdivision hyperplanes
	  if ( PSubdivision(i,c) != 0 ) { // otherwise it's already projected
	    
	    // FIXME 
	    // determine whether intersection fo subdivision hyperplanes intersecs inside polytope
	    // should be done via a 1D linear program
	    perl::Object spp("Polytope<Rational>");
	    perl::Object spm("Polytope<Rational>");
	    spp.take("INEQUALITIES") << PF/PSubdivision.row(i);
	    spm.take("INEQUALITIES") << PF/(-PSubdivision.row(i));
	    Matrix<Rational> VSP = spp.give("VERTICES");
	    Matrix<Rational> VSM = spm.give("VERTICES");
	    
	    
	    for ( int j = 0; j < PSubdivision.rows(); ++j )  		// combine with subdivision hyperplanes
	      if ( PSubdivision(j,c) != 0 ) // otherwise it's already projected
		if ( hyperplane_intersects(PSubdivision.row(j),VSP) &&  hyperplane_intersects(PSubdivision.row(j),VSM) ) {
		  Vector<Integer> w(QF.cols());
		  if ( PSubdivision(i,c) == PSubdivision(j,c) ) 
		    w = (PSubdivision.row(i) - PSubdivision.row(j)).slice(~range(c,c));
		  else
		    w = (PSubdivision.row(i) + PSubdivision.row(j)).slice(~range(c,c));
		  if ( hyperplane_intersects(w,QV) )
		    QSubdivision /= w;
		}
	    
	  
	    for ( int j = 0; j < PF.rows(); ++j ) // combine with subdivision hyperplanes
	      if ( PF(j,c) != 0 ) // otherwise it's already projected
		if ( hyperplane_intersects(  PSubdivision.row(i), SP.vertices().minor(SP.vertices_in_facets().row(j),All) ) ) {
		  Vector<Integer> w(QF.cols());
		  if ( PSubdivision(i,c) == PF(j,c) ) 
		    w = (PSubdivision.row(i) - PF.row(j)).slice(~range(c,c));
		  else
		    w = (PSubdivision.row(i) + PF.row(j)).slice(~range(c,c));
		  if ( hyperplane_intersects(w,QV) )
		    QSubdivision /= w;
		}
	  }
	
	// combine facets with facets, but only adjacent ones!
	for ( int i = 0; i < PF.rows(); ++i ) 
	  if ( PF(i,c) != 0 )  // otherwise it's already projected
	    for ( int j = 0; j < PF.rows(); ++j )  // combine with subdivision hyperplanes
	      if ( SP.dual_graph().edge_exists(i,j) && PF(j,c) != 0 ) {
		Vector<Integer> w(QF.cols());
		if ( PF(i,c) == PF(j,c) ) 
		  w = (PF.row(i) - PF.row(j)).slice(~range(c,c));
		else
		  w = (PF.row(i) + PF.row(j)).slice(~range(c,c));
		if ( hyperplane_intersects(w,QV) ) 
		  QSubdivision /= w;
	      }

	return SubdividedPolytope(QF,QV,make_equations_unique(primitive(QSubdivision)));
	
      }

    }
    
    // computes the cells of a subdivision and matches the vertices with the lattice points of the polytope
    Array< Set<int> > subdivision_cells (  const Matrix<Rational> & V, 
					   const Matrix<Integer> & F, 
					   const Matrix<Integer> & S,
					   const Matrix<Integer> & L ) {
      
      Vector<subdivision::SubdivisionCell> cells = subdivision::compute_cells ( V, F, S );
      Array<Set<int> > lp_in_cells(cells.dim());
      int j = 0;
      for ( Entire<Vector<subdivision::SubdivisionCell> >::const_iterator c = entire( cells );  !c.at_end(); ++c ) {
	Set<int> S;
	for ( Entire<Rows<Matrix<Rational> > >:: const_iterator v = entire(rows((*c).vertices()));
	      !v.at_end(); ++v ) {
	  int i = 0;
	  while ( i < L.rows() && *v != L.row(i) ) ++i;
	  if ( i == L.rows() ) {
	    exit(1);
	  } else {
	    S += i;
	  }
	}
	lp_in_cells[j++] = S;
      }

      return lp_in_cells;
    }

    
    // computes the cells of a subdivision 
    // and checks whether all cells have integral vertices
    bool lattice_subdivision(const Matrix<Rational> & V, 
			     const Matrix<Integer> & F, 
			     const Matrix<Integer> & S ){
      
      const Vector<subdivision::SubdivisionCell> cells = subdivision::compute_cells ( V, F, S );
      
      for ( Entire<Vector<subdivision::SubdivisionCell> >::const_iterator c = entire( cells ); 
	    !c.at_end(); ++c ) 
	if ( !(*c).integral() ) {
	  return false;
	}
      
      return true;
    } //  end of lattice_subdivision
    
    
  } }
