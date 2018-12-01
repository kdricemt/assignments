/*******************************************
C++ Program for 1D linear-hyperbolic Problem
********************************************/

#include <cstdio>
#include <vector>
#include <cstdlib>
#include <iostream>
#include <string>
#include <fstream>

std::vector<double> x(100,0);
std::vector<double> yi(100,0); //initial condition
std::vector<double> ye(100,0); //exact solution
std::vector<double> ys(100,0); //numerical solution
std::vector<double> work(100,0);

int mx = 99; // number of nodes
int nstep = 30; // number of steps
double c = 1.0;
double cfl = 0.5; //Courant number
double dx = 0.01;  // =1/mx
double dt = cfl*dx;

void exact(int i){
   double travel = double(nstep) * dt*c;
   if(x[i] < (0.5 + travel)) ye[i] = 1.0;
   else ye[i] = 0.0;
}

void scheme_1(int i){
   work[i] = ys[i] - 0.5*cfl*(ys[i+1]-ys[i-1]);
}

void scheme_2(int i){
   work[i] = ys[i] - cfl*(ys[i]-ys[i-1]);
}

void lax_wendroff(int i){
   work[i] = 0.5*cfl*(1+cfl)*ys[i-1]+(1-cfl*cfl)*ys[i] - 0.5*cfl*(1-cfl)*ys[i+1];
}


void numerical_sol(int scheme_num,int i){
  switch(scheme_num){
    case 1:
      scheme_1(i);
      break;
    case 2:
      scheme_2(i);
      break;
    case 3:
      lax_wendroff(i);
      break;
    default:
      lax_wendroff(i);
  }
}

int main(){
  int scheme_num =0;
  std::string scheme_str;
  std::cout << "Select Scheme Number (1:中心差分,2:風上差分,3:lax_wendroff) :" << std::endl;
  std::cin >> scheme_num;
  std::cout << "Selected Scheme: ";
  switch(scheme_num){
    case 1:
      scheme_str = "center";
      break;
    case 2:
      scheme_str = "upper_stream";
      break;
    case 3:
      scheme_str = "Lax_Wendroff";
      break;
    default:
      scheme_str = "Lax_Wendroff";
  }
  std::cout << scheme_str << std::endl;

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
   for(int n=0; n < nstep; n++){
     for(int i =0; i < mx; i++){
        numerical_sol(scheme_num,i);
     }
     for(int i=0; i<mx; i++){
        ys[i] = work[i];
     }
   }
   // write datafile

   // exact solution
   std::ofstream fout_ex;
   std::string filename_ex = "../datafile/1Dhyperbolic/exact.txt";
   fout_ex.open(filename_ex.c_str(), std::ios::out);
   for (int i=0; i<mx; i++) {
       fout_ex << i*dx << " " << ye[i] << std::endl;
   }
   fout_ex.close();

   // numerical sol
   std::ofstream fout_num;
   std::string filename_num = "../datafile/1Dhyperbolic/" + scheme_str + ".txt";
   fout_num.open(filename_num.c_str(), std::ios::out);
   for (int i=0; i<mx; i++) {
       fout_num << i*dx << " " << ys[i] << std::endl;
   }
   fout_num.close();
}
