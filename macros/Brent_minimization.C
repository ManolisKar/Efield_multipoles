#include "Math/Functor.h"
#include "Math/BrentMinimizer1D.h"
#include "Math/GSLMinimizer1D.h"

double myfunc(double x ) { 
   return 1 + -4*x + 1*x*x; 
}
 
int Brent_minimization(){
   ROOT::Math::Functor1D func(&myfunc);
 
   ROOT::Math::BrentMinimizer1D bm;
   bm.SetFunction(func, -10,10);
   bm.Minimize(100,0,0);
 
   cout << "f(" << bm.XMinimum() << ") = " <<bm.FValMinimum() << endl;
 
   return 0;
}




int GSL_minimization()
{
   ROOT::Math::Functor1D func(&myfunc);
 
   // default (Brent)
   ROOT::Math::GSLMinimizer1D minBrent;
   minBrent.SetFunction(func,1,-10,10);
   minBrent.Minimize(100,0.01,0.01);
   std::cout << "Found minimum: x = " << minBrent.XMinimum() 
             << "  f(x) = " << minBrent.FValMinimum() << std::endl;
 
   // Golden Section
   ROOT::Math::GSLMinimizer1D minGold(ROOT::Math::Minim1D::kGOLDENSECTION);
   minGold.SetFunction(func,1,-10,10);
   minGold.Minimize(100,0.01,0.01);
   std::cout << "Found minimum: x = " << minGold.XMinimum() 
             << "  f(x) = " << minGold.FValMinimum() << std::endl;
 
   return 0;
}