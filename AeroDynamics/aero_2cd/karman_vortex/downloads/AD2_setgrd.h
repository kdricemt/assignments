#ifndef AD2_SetGrd_H
#define AD2_SetGrd_H

using namespace std;

class SetGrd{
  public:
    static const int mx=401, i1=96, i2=106
                   , my=201, j1=96, j2=106;
    static constexpr double dx=1/double(i2-i1), dy=1/double(j2-j1);
    static double setgrd(vector<vector<double>>& x, vector<vector<double>>& y);
};

double SetGrd::setgrd(vector<vector<double>>& x, vector<vector<double>>& y){
  double icent=0.5*(i1+i2), jcent=0.5*(j1+j2);
  for(int i=0;i<=mx;i++){
    for(int j=0;j<=my;j++){
      x[i][j]=dx*double(i-icent);
      y[i][j]=dy*double(j-jcent);
    }
  }

  if(dx>=dy){return dx;}else{return dy;}
}
#endif
