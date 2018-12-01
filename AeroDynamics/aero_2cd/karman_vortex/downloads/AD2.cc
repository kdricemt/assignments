/*******************************************
 * Incompressible Navier-Stokes 2D Flow Solver
 * "Flow around Rectangular Cylinder"
 * Author: K. Shu (160328, U-Tokyo)
 * Language: C++11
 * Compiler: Clang++
 * Date: 
********************************************/

#include <iostream>
#include <cmath>
#include <vector>
#include <fstream>

#include "AD2_setgrd_3.h"
#include "AD2_slvflw_3.h" //CHANGE FILE NAME!!!!

using namespace std;


/***** Global Var *****/
const extern double re;  //Reynolds No.
extern double rei;  //Re inverse

extern int nlast, nlp;  //No. of time steps
/***** End Global Var *****/


/***** Main *****/
int main(){
  
  ios_base::sync_with_stdio(false);

  double cfl=0.2; //CFL
  int mdx=SetGrd::mx+100, mdy=SetGrd::my+100; //size of x[i][j], y[i][j]
  vector<vector<double>> x(mdx, vector<double>(mdy)),
                         y(mdx, vector<double>(mdy)),
                         u(mdx, vector<double>(mdy)),
                         v(mdx, vector<double>(mdy)),
                         p(mdx, vector<double>(mdy)); 

  double dt=SetGrd::setgrd(x, y)*cfl; //call setgrd(return min(dx, dy))

  cout<<">>2D Incompressible Flow Solver\n"
      <<">>Re= "<<re<<"\n"
      <<">>No. of Grid Points: (x, y)= ("<<SetGrd::mx<<", "<<SetGrd::my<<")\n"
      <<">>CFL/dt/Steps= "<<cfl<<"/"<<dt<<"/"<<nlast<<endl;

  SlvFlw::slvflw(dt,x,y,u,v,p); //call slvflw

return 0;
}
/***** End Main *****/
