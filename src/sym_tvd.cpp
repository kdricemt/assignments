/*********************************************
C++ Program for 1D linear-hyperbolic problem
Symmetric-TVD scheme(minmod,superbee)
**********************************************/

#include <cstdio>
#include <cmath>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>
#include <sstream>

namespace std {
  template <typename T>
  std::string to_string(T value)
  {
    //create an output string stream
    std::ostringstream os ;

    //throw the value into the string stream
    os << value ;

    //convert the string stream into a string and return
    return os.str() ;
  }
}

std::vector<double> x(101,0);
std::vector<double> yi(101,0); //initial condition
std::vector<double> ye(101,0); //exact solution
std::vector<double> ys(101,0); //numerical solution
std::vector<double> flux(101,0);
std::vector<double> work(101,0);

double c = 1.0;
int mx = 51; //number of nodes
int nlast = 30; //number of steps
double cfl = 0.5;
double dx = 0.02; // = 1/(mx-1)
double dt = cfl * dx;
double travel = dt * c * double(nlast);
double ecp = 0.1;

double sign(double a, double b){
  int bplus = 0;
  if(b >= 0) bplus = 1;
  else bplus = -1;
  return std::abs(a) * bplus;
}

void exact(int i){
   if(x[i] < (0.5 + travel)) ye[i] = 1.0;
   else ye[i] = 0.0;
}

void minmod1(int i){
  double dm = ys[i] - ys[i-1];
  double d0 = ys[i+1] - ys[i];
  double dp = ys[i+2] - ys[i+1];
  double s = sign(1.0,dm);
  double q = s * std::max(0.0,std::min(std::min(s*dm,s*d0),s*dp));
  double ac = std::abs(c);
  if(ac < ecp) ac = (c*c + ecp*ecp)*0.5/ecp;
  double f = - (dt*c*c/dx*q + ac*(d0-q));
  flux[i] = 0.5*(c*ys[i] + c*ys[i+1] + f);
}

void minmod2(int i){
  double dm = ys[i] - ys[i-1];
  double d0 = ys[i+1] - ys[i];
  double dp = ys[i+2] - ys[i+1];
  double s = sign(1.0,dm);
  double q = s * std::max(0.0,std::min(std::min(std::min(2*s*dm,2*s*d0),2*s*dp),0.5*s*(dm+dp)));
  double ac = std::abs(c);
  if(ac < ecp) ac = (c*c + ecp*ecp)*0.5/ecp;
  double f = - (dt*c*c/dx*q + ac*(d0-q));
  flux[i] = 0.5*(c*ys[i] + c*ys[i+1] + f);
}

void superbee(int i){
  double dm = ys[i] - ys[i-1];
  double d0 = ys[i+1] - ys[i];
  double dp = ys[i+2] - ys[i+1];
  double s = sign(1.0,dm);
  double sb1 = s * std::max(0.0,std::max(std::min(2*std::abs(dm),s*d0),std::min(std::abs(dm),2*s*d0)));
  double sb2 = s * std::max(0.0,std::max(std::min(2*std::abs(d0),s*dp),std::min(std::abs(d0),2*s*dp)));
  double q = sb1 + sb2 - d0;
  double ac = std::abs(c);
  if(ac < ecp) ac = (c*c + ecp*ecp)*0.5/ecp;
  double f = - (dt*c*c/dx*q + ac*(d0-q));
  flux[i] = 0.5*(c*ys[i] + c*ys[i+1] + f);
}


void numerical_sol(int scheme_num,int i){
  switch(scheme_num){
    case 1:
      minmod1(i);
      break;
    case 2:
      minmod2(i);
      break;
    case 3:
      superbee(i);
      break;
    default:
      superbee(i);
  }
}

int main(){
  int scheme_num =0;

  std::string scheme_str;
  std::cout << "Symmetric- TVD scheme" << std::endl;
  std::cout << "Select Scheme Number (1:minmod1,2:minmod2,3:superbee) :" << std::endl;
  std::cin >> scheme_num;
  std::cout << "Selected Scheme: ";
  switch(scheme_num){
    case 1:
      scheme_str = "minmod1";
      break;
    case 2:
      scheme_str = "minmod2";
      break;
    case 3:
      scheme_str = "superbee";
      break;
    default:
      scheme_str = "superbee";
  }
  std::cout << scheme_str << std::endl;
  std::cout << "Input Entropy Correction Parameter :" << std::endl;
  std::cin >> ecp;
  std::cout << "ECP: " << ecp << std::endl;

  // set initial condtition
  for(int i=0;i<mx;i++){
     x[i] = dx * double(i);
     if(x[i] < 0.5) yi[i] = 1.0;
     else yi[i] = 0.0;
   }
   // set exact solution
   for (int i=0;i<mx;i++){
     exact(i);
   }
   // numerical solution
   for(int i=0; i< mx; i++){
      ys[i] = yi[i];
   }
   for(int n=0; n < nlast; n++){
     for(int i =1; i < mx - 1; i++){
        numerical_sol(scheme_num,i);
     }
     for(int i=1; i<mx-1; i++){
        work[i] = ys[i] - cfl*(flux[i]-flux[i-1]);
     }
     for(int i=1; i<mx-1; i++){
        ys[i] = work[i];
     }
  }

   /*** write datafile ***/
   // exact solution
   std::ofstream fout_ex;
   std::string filename_ex = "../datafile/sym_tvd/exact.txt";
   fout_ex.open(filename_ex.c_str(), std::ios::out);
   for (int i=0; i<mx; i++) {
       fout_ex << i*dx << " " << ye[i] << std::endl;
   }
   fout_ex.close();

   // numerical sol
   std::ofstream fout_num;
   std::string filename_num = "../datafile/sym_tvd/" + scheme_str + "_ecp" + std::to_string(ecp) + ".txt";
   fout_num.open(filename_num.c_str(), std::ios::out);
   for (int i=0; i<mx; i++) {
       fout_num << i*dx << " " << ys[i] << std::endl;
   }
   fout_num.close();
}
