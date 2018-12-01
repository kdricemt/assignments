//翼周りの流れのシミュレーション

#include<iostream>
#include<cstdio>
#include<cmath>
#include<complex>
#include<vector>
using namespace std;

double PI = 3.141592;

FILE *gid;
complex<double> cv,ca,cgama,ccent,small_z,large_z,cf_potential,cf_pressure;
double a = 0.6;
double v0 = 1.0;
double c025 = 0.5;
double gama = 5.0;

double x_z[180][100]; //Z平面におけるx
double y_z[180][100]; //Z平面におけるy
double x_zeta[180][100]; //zeta平面におけるx
double y_zeta[180][100]; //zeta平面におけるy
double f_stream[180][100]; //流線関数
double f_cp[180][100]; //圧力Cpの分布

//計算用関数(空力中心以外)
void calculation(double alfa,double beta);
  void cal_z(int i,int j,double alfa,double beta);
  void grid_on_zetaplane(int i,int j);
  void cal_stream(int i,int j);
  void cal_pressure(int i,int j,double alfa);
void cal_force(double alfa,double beta);
//ファイルデータ書き込み用関数
void write_data();
  void write_shape();
  void write_stream();
  void write_pressure();
  void write_surfacepressure();
//空力中心計算
void write_ac(double beta);
  double cal_cm(double x_sample);

//------------------------------------------------main
int main(){
  int wing;
  double rad = (PI/180.0);
  double alfa,beta;
  printf("choose Jekofusuki(0) or symmetric(1):\n");
  scanf("%d",&wing);
  if(wing == 0){
  //ジェコフスキー翼
    alfa = 5.0 * rad;
    beta = 20.0 * rad;
    calculation(alfa,beta);
    write_data();
    cal_force(alfa,beta); //Cl,CD,風圧中心を計算
    write_ac(beta); //空力中心を計算
  }
  else{
  //対称翼
    alfa = 5.0 * rad;
    beta = 0.0 * rad;
    calculation(alfa,beta);
    write_data();
    cal_force(alfa,beta); //Cl,CD,風圧中心を計算
    write_ac(beta); //空力中心を計算
  }
}

//--------------------------------------------------計算用関数(空力中心以外)
/*
考える領域をζではなくZで定義する
Z平面上で同心円上に計算点を配置。j=1がZ面、ζ面共に物体表面に対応している。
複素ポテンシャル(cf)とともに位置ζも求めている
*/
void calculation(double alfa,double beta){
  //Kutta condition
  gama = 4.0 * PI * v0 * a * sin(alfa+beta);

  for(int i=1;i<180;i++){
    for(int j=0;j<100;j++){
      //small_z,large_z
      cal_z(i,j,alfa,beta);
      //zeta-plain
      grid_on_zetaplane(i,j);
      //stream_function
      cal_stream(i,j);
      //pressure
      cal_pressure(i,j,alfa);
   }
  }
}

//large_z,small_zの計算
void cal_z(int i,int j,double alfa,double beta){
  cv = complex<double>(v0,0.0);
  ca = complex<double>(a,0.0);
  cgama = complex<double>(0.0,gama/(2.0*PI));
  //grid on Z-plane
  double dr = 0.05;
  double dtheta = 2.0*(PI/180.0); //rad
  double rrr = a + dr * ((double)j); //大きさ
  double theta0 = -beta;
  double theta = theta0 + dtheta * ((double)i); //角度
  //プリントのzexp(iα)に対応(Z平面上の同心円)
  x_z[i][j] = rrr * cos(theta);
  y_z[i][j] = rrr * sin(theta);
  //プリントのZc
  ccent = complex<double>(c025,0.0) + ca * exp(complex<double>(0.0,PI - beta));
  //プリントのZに対応
  large_z = complex<double>(x_z[i][j],y_z[i][j])+ccent;
  //Step2のz
  small_z = (large_z-ccent) * exp(complex<double>(0.0,-alfa));
}

//ζ平面のx-y関係の計算 j=0が翼表面
void grid_on_zetaplane(int i,int j){
  x_zeta[i][j] = real(large_z + complex<double>(c025*c025,0.0)/large_z);
  y_zeta[i][j] = imag(large_z + complex<double>(c025*c025,0.0)/large_z);
}

//流線関数の計算
void cal_stream(int i,int j){
  //プリントstep2のf
  cf_potential = cv * (small_z + ca*ca/small_z) + cgama * log(small_z);
  f_stream[i][j] = imag(cf_potential);
}

//圧力の計算
void cal_pressure(int i,int j,double alfa){
  cf_pressure = exp(complex<double>(0.0,-alfa))
      *(cv * (complex<double>(1.0,0.0)-ca*ca/(small_z*small_z)) + cgama/small_z)
      /(complex<double>(1.0,0.0)
      - complex<double>(c025,0.0) * complex<double>(c025,0.0) / (large_z * large_z));
  double cp = 1.0
      - (real(cf_pressure) * real(cf_pressure) + imag(cf_pressure)*imag(cf_pressure))/(v0 * v0);
  f_cp[i][j] = cp;
  if(i==0 && j==0){
  printf("real_cp:%f,img_cp:%f\n",real(cf_pressure),real(cf_pressure));
  printf("f_cp[0][0]:%f\n",cp);
  }
}

//Cl,Cd,風圧中心の計算
void cal_force(double alfa,double beta){
  calculation(alfa,beta);
  double cxp = 0.0;
  double cyp = 0.0;
  double sum = 0.0;
  for(int i=1;i<179;i++){
    double dxw = x_zeta[i+1][0] - x_zeta[i][0];
    double dyw = y_zeta[i+1][0] - y_zeta[i][0];
    double dnx = dyw;
    double dny = -dxw;
    double cpm = (f_cp[i+1][0]+f_cp[i][0])/2.0;
    cxp = cxp - cpm * dnx;
    cyp = cyp - cpm * dny;
    double fx = -cpm * dnx;
    double fy = -cpm * dny;
    sum = sum + (fy * x_zeta[i][0] - fx * y_zeta[i][0]);
  }
  double center_of_pressure = sum/cyp; //x座標
  cxp = cxp / (4.0 * c025);
  cyp = cyp / (4.0 * c025);
  double cdp = cxp*cos(alfa) + cyp*sin(alfa);
  double clp = cyp*cos(alfa) - cxp*sin(alfa);

  printf("CL = %f\n",clp);
  printf("CD = %f\n",cdp);
  printf("center_of_pressure=%f(前縁からの距離)\n",center_of_pressure/(4.0 * c025) + 0.5);
}

//---------------------------------------------データのファイル書き込み用関数

void write_data(){
    write_shape();
    write_stream();
    write_pressure();
    write_surfacepressure();
}

//翼型の表示
void write_shape(){
  FILE *fid;
  const char *data={"xy_zeta.txt"};
  fid = fopen(data,"w");
    for(int p=1;p<180;p++){
      fprintf(fid,"%f\t%f\n",x_zeta[p][0],y_zeta[p][0]);
    }
}

//流線関数の表示
void write_stream(){
  FILE *fid2;
  const char *data2={"streamline.txt"};
  fid2 = fopen(data2,"w");
  for(int i=1;i<180;i++){
    for(int j=0;j<100;j++){
      fprintf(fid2,"%f\t%f\t%f\n",x_zeta[i][j],y_zeta[i][j],f_stream[i][j]);
   }
   fprintf(fid2,"\n");
 }
}

//圧力の表示
void write_pressure(){
  FILE *fid3;
  const char *data3={"cp.txt"};
  fid3 = fopen(data3,"w");
  for(int i=1;i<180;i++){
    for(int j=0;j<100;j++){
      fprintf(fid3,"%f\t%f\t%f\n",x_zeta[i][j],y_zeta[i][j],f_cp[i][j]);
   }
   fprintf(fid3,"\n");
 }
}

//翼面上の圧力を表示
void write_surfacepressure(){
  FILE *fid4;
  const char *data4={"cp_surface.txt"};
  fid4 = fopen(data4,"w");
    for(int i=1;i<180;i++){
      fprintf(fid4,"%f\t%f\n",x_zeta[i][0],f_cp[i][0]);
    }
  fprintf(fid4,"\n");
}

//-------------------------------------------空力中心計算用関数

//空力中心計算用のデータの書き込み
void write_ac(double beta){
  FILE *fid5;
  const char *data5={"aerodynamic_center.txt"};
  fid5 = fopen(data5,"w");
  for(int i=0;i<5;i++){
    double alfa = i*2*(PI/180.0);
    calculation(alfa,beta);
    double x_sample = -0.9; //x=-0.9-0.9
    for(int j=0;j<18;j++){
      x_sample = x_sample + 0.1;
      double cm = cal_cm(x_sample);
      fprintf(fid5,"%f\t%f\n",x_sample/(4.0*c025),cm);
    }
    fprintf(fid5,"\n");
  }
}

//参照点(x_sample,0)周りピッチングモーメントの計算
double cal_cm(double x_sample){
  double cm = 0.0;
  for(int i=1;i<179;i++){
    double dxw = x_zeta[i+1][0] - x_zeta[i][0];
    double dyw = y_zeta[i+1][0] - y_zeta[i][0];
    double dnx = dyw;
    double dny = -dxw;
    double cpm = (f_cp[i+1][0]+f_cp[i][0])/2.0;
    double fx = -cpm * dnx;
    double fy = -cpm * dny;
    cm = cm + (fy * (x_sample - x_zeta[i][0]) + fx * y_zeta[i][0]);
  }
  return cm;
}
