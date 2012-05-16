

#include <polymake/Rational.h>
#include <polymake/Integer.h>
#include <polymake/client.h>
#include <polymake/Matrix.h>
#include <polymake/Array.h>

namespace polymake { namespace polytope { 

    std::pair<int , std::pair<int, Array<Array<Integer> > > >  
    pair_test() {

      std::pair<int , std::pair<int, Array<Array<Integer> > > > p;
      std::pair<int, Array<Array<Integer> >  > q;
      Array<Array<Integer> > r(2);
      q.first=1;
      q.second = r;
      p.first=3;
      p.second=q;
      return p;
    }

    UserFunction4perl( "# ", &pair_test, "pair_test()");

  }}
