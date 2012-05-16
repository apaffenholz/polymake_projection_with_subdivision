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
#include <polymake/Integer.h>
#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Vector.h>
#include <polymake/IncidenceMatrix.h>
#include <polymake/Graph.h>
#include <polymake/Array.h>
#include <polymake/Set.h>
#include <polymake/polytope/subdivisions.h>
#include <polymake/polytope/SubdividedPolytope.h>
#include <polymake/linalg.h>

using std::list;

namespace polymake { namespace polytope { 
    
    typedef list<list<int> > PSequence;
    
    namespace subdivision {
      
      bool projection_step ( SubdividedPolytope SP, const int c, 
			     const int startdim, PSequence & pseq, list<int> seq, 
			     int & lowest, const bool & only_lowest ) {      
	
	// do a projection along c
	SubdividedPolytope SQ = project_with_subdivision_step( SP, c);
	
	Matrix<Integer> F = SQ.facets();
	Matrix<Integer> S = SQ.subdivision_hyperplanes();
	bool ls = lattice_subdivision(SQ.vertices(), F, S);
	
	if ( ls ) {
	  list<int> seq2(seq);
	  seq2.push_back(c);
	  if ( SQ.dim() < lowest ) {
	    lowest = SQ.dim();
	    pseq.clear();
	  }
	  
	  Set<int> VPC = SQ.valid_projection_coordinates();
	  
	  if ( SQ.dim() > 1 && VPC.size() > 0 ) {
	    for(Entire<Set<int> >::const_iterator vpc = entire(VPC); !vpc.at_end(); ++vpc ) {
	      if ( projection_step(SQ,*vpc,startdim,pseq,seq2,lowest, only_lowest) )
		break;
	    }
	  } else {
	    if ( SQ.dim() == lowest )
	      pseq.push_back(seq2);
	    if ( only_lowest && lowest == 1 )
	      return true;
	  }
	} else 
	  if ( SP.dim() == lowest ) 
	    pseq.push_back(seq);
	
	return ( only_lowest && lowest == 1 );
      }
    }
    
    std::pair<int , std::pair<int, Array<Array<Integer> > > > 
    check_all_projections(perl::Object p, bool only_lowest = false, bool all_facets = false) {
      
      int adim = p.give("CONE_AMBIENT_DIM");
      int dim = p.give("CONE_DIM");
      if ( dim != adim ) {
	throw std::runtime_error("project_with_subdivision: polytope must be full dimensional\n");       
      }
      
      Matrix<Rational> SH(0,adim-1);	
      if ( p.lookup("SUBDIVISION") ) 
	p.lookup("SUBDIVISION.SUBDIVISION_HYPERPLANES") >> SH;
      
      Matrix<Rational> F = p.give("FACETS");
      Matrix<Rational> V = p.give("VERTICES");
      IncidenceMatrix<> VIF = p.give("VERTICES_IN_FACETS");
      Graph<> DG = p.give("DUAL_GRAPH.ADJACENCY");

      std::pair<int, std::pair<int,Array<Array<int> > > > pair_pseq;

      
      if ( !all_facets ) {
	subdivision::SubdividedPolytope SP(F, V, SH, VIF, DG);
	
	PSequence pseq;
	
	Set<int> VPC =  SP.valid_projection_coordinates();
	int lowest = SP.dim();
	
	for(Entire<Set<int> >::const_iterator vpc = entire(VPC); 
	    !vpc.at_end(); ++vpc ) {
	  list<int> seq;
	  if ( subdivision::projection_step(SP,*vpc,dim-1,pseq,seq,lowest, only_lowest) )
	    break;
	}

	pair_pseq.first = lowest;
	pair_pseq.second.first = -1;
	pair_pseq.second.second = pseq;

      } else {
	
	IncidenceMatrix<> ftv = p.give("FACETS_THRU_VERTICES");
	
	pair_pseq.first = F.cols();

	for ( Entire<Rows<IncidenceMatrix<> > >::const_iterator ftv_it = entire(rows(ftv)); !ftv_it.at_end(); ++ftv_it ) {
	 
	  Matrix<Rational> MT = vector2row(unit_vector<Rational>(F.cols(),0))/
                	       (zero_vector<Rational>(F.cols()-1)|(-F.minor(*ftv_it,~range(0,0))));
	  Matrix<Rational> IT = inv(MT);

	  Matrix<Rational> NF = F * IT;
	  Matrix<Rational> NSH = SH * IT;
	  Matrix<Rational> NV = V * T(MT);

	  subdivision::SubdividedPolytope SP(NF, NV, NSH, VIF, DG);
	  
	  PSequence pseq;
	  
	  Set<int> VPC =  SP.valid_projection_coordinates();
	  int lowest = SP.dim();
	  
	  for(Entire<Set<int> >::const_iterator vpc = entire(VPC); 
	      !vpc.at_end(); ++vpc ) {
	    list<int> seq;
	    if ( subdivision::projection_step(SP,*vpc,dim-1,pseq,seq,lowest, only_lowest) )
	      break;
	  }

	  if ( lowest < pair_pseq.first ) {
	    pair_pseq.first = lowest;
	    pair_pseq.second.first = ftv_it->index();
	    pair_pseq.second.second = pseq;
	  }
	  
	  if ( lowest == 1 )
	    break;
	}

      }

      return pair_pseq;
    }

    Function4perl( &check_all_projections, "check_all_projections(Polytope; $=0, $=0)");

  } }

