//翼周りの流れのシミュレーション
//air.cpp

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

double windpressure_center; //風圧中心のx座標
double cm; //ピッチングモーメント係数
double aerodynamic_center; //空力中心のx座標

typedef struct point{
  double x;
  double y;
}Point;

int open_gnuplot();
int close_gnuplot();
void caluculation(double alfa,double beta);
  void cal_z(int i,int j,double alfa,double beta);
  void grid_on_zetaplane(int i,int j);
  void cal_stream(int i,int j);
  void cal_pressure(int i,int j,double alfa);
void cal_force(double alfa,double beta);
void cal_ac();
  double cal_cm();
void show_data(int mode);
  void shape_data();
  void show_stream();
  void show_pressure();
  void show_surfacepressure();
  void make_contour(double f[][100],double min,double max,double df);

//------------------------------------------------gnu_plot
//gnu_plotにパイプを通す、初期設定
int open_gnuplot(){
  gid = popen("gnuplot", "w");
  if (gid == NULL)
    return -1;
  fprintf(gid,"set angles radians\n"); //角度の単位をラジアンに設定
  return 0;
}

int close_gnuplot(){
  fflush(gid);
  pclose(gid);
  return 0;
}

//--------------------------------------------------計算用関数
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
      cal_pressure(i,j,double alfa);
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
  //プリントのlarge Zに対応
  large_z = complex<double>(x_z[i][j],y_z[i][j])+ccent;
  //Step2のsmall_z
  small_z = (large_z-ccent) * exp(complex<double>(0.0,-alfa));
}

//ζ平面のx-y関係の計算
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
  windpressure_center = sum/cyp; //x座標
  cxp = cxp / (4.0 * c025);
  cyp = cyp / (4.0 * c025);
  double cdp = cxp*cos(alfa) + cyp*sin(alfa);
  double clp = cyp*cos(alfa) - cxp*sin(alfa);

  printf("CL = %f\n",clp);
  printf("CD = %f\n",cdp);
  printf("windpressure_center=%f\n",windpressure_center/(4.0 * c025) - 0.5);
}

//参照点(x_sample,0)周りピッチングモーメントの計算
double cal_cm(double x_sample){
  cm = 0.0;
  for(int i=1;i<179;i++){
    double dxw = x_zeta[i+1][0] - x_zeta[i][0];
    double dyw = y_zeta[i+1][0] - y_zeta[i][0];
    double dnx = dyw;
    double dny = -dxw;
    double cpm = (f_cp[i+1][0]+f_cp[i][0])/2.0;
    double fx = -cpm * dnx;
    double fy = -cpm * dny;
    cm = cm + (fy * (x_sample - x_zeta[i][0]) - fx * y_zeta[i][0]);
    return cm;
  }
}

//空力中心の計算
void cal_ac(double beta){
  FILE *fid7;
  const char *data7={"aerodynamic_center.txt"};
  fid7 = fopen(data7,"w");
  for(int i=0;i<5;i++){
    alfa = i*2*(PI/180.0);
    calculation(alfa,beta);
    double x = -0.9; //x=-0.9-0.9
    for(int j=0;j<18;j++){
      x = x + 0.1;
      double cm = cal_cm(double x);
      fprintf(fid7,"%f\t%f\n",x,cm);
    }
    fprintf(fid7,"\n");
  }
}

//--------------------------------------------------------表示用関数

void show_data(int mode){
  switch(mode){
    case 0 : shape_data(); break;
    case 1 : show_stream();  break;
    case 2 : show_pressure(); break;
    case 3 : show_surfacepressure(); break;
  }
}
//翼型の表示
void shape_data(){
  FILE *fid;
  const char *data={"xy_zeta.txt"};
  fid = fopen(data,"w");
    for(int p=1;p<180;p++){
      fprintf(fid,"%f\t%f\n",x_zeta[p][0],y_zeta[p][0]);
    }
}

//流線関数の表示
void show_stream(){
  FILE *fid2;
  const char *data2={"streamline.txt"};
  fid2 = fopen(data2,"w");
  for(int i=1;i<180;i++){
    for(int j=0;j<100;j++){
      fprintf(fid2,"%f\t%f\t%f\n",x_zeta[i][j],y_zeta[i][j],f_stream[i][j]);
   }
   fprintf(fid2,"\n");
 }
 fprintf(fid2,"\n");
 fprintf(gid, "set xrange [-2:2]\n");
 fprintf(gid, "set yrange [-2:2]\n");
 fprintf(gid, "set zrange [-0.5:1]\n");
 fprintf(gid, "set dgrid3d 150,150 qnorm 3\n");
 fprintf(gid, "set contour\n");
 fprintf(gid, "set cntrparam levels incremental -2.0,0.1,1.0\n");
 //fprintf(gid, "set cntrparam levels discrete -0.5,-0.4,-0.3,-0.2,0,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0,1,1,1.2\n");
 fprintf(gid, "set view 0,0\n");
 fprintf(gid, "set key off\n");
 fprintf(gid, "unset surface\n");
 fprintf(gid, "set table 'streamline.table'\n");
 //fprintf(gid, "splot 'streamline.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7\n");
 fprintf(gid, "splot 'streamline.txt' using 1:2:3 with lines palette\n");
 fprintf(gid, "unset table\n");
}

//圧力の表示
void show_pressure(){
  FILE *fid3;
  const char *data3={"cp.txt"};
  fid3 = fopen(data3,"w");
  for(int i=1;i<180;i++){
    for(int j=0;j<100;j++){
      fprintf(fid3,"%f\t%f\t%f\n",x_zeta[i][j],y_zeta[i][j],f_cp[i][j]);
   }
   fprintf(fid3,"\n");
 }
 fprintf(fid3,"\n");
 fprintf(gid, "set xrange [-2:2]\n");
 fprintf(gid, "set yrange [-2:2]\n");
 fprintf(gid, "set zrange [-2:2]\n");
 //fprintf(gid, "set view 0,0\n");
 fprintf(gid, "set dgrid3d 100,100 qnorm 3\n");
 fprintf(gid, "set contour\n");
 fprintf(gid, "set cntrparam levels incremental -2.0,0.1,1.0\n");
 fprintf(gid, "set view 0,0\n");
 fprintf(gid, "set key off\n");
 fprintf(gid, "unset surface\n");
 //fprintf(gid, "splot 'streamline.txt' using 1:2:3 with points palette pointsize 0.5 pointtype 7\n");
 fprintf(gid, "set table 'pressure.table'\n");
 fprintf(gid, "splot 'cp.txt' using 1:2:3 with lines palette\n");
 fprintf(gid, "unset table\n");
}

//翼面上の圧力を表示
void show_surfacepressure(){
  FILE *fid5;
  const char *data5={"cp_surface.txt"};
  fid5 = fopen(data5,"w");
    for(int i=1;i<180;i++){
      fprintf(fid5,"%f\t%f\n",x_zeta[i][0],f_cp[i][0]);
    }
  fprintf(fid5,"\n");
}

//等高線作成 fは対象とする関数、[min,max]は範囲、dfは刻み幅
/*void make_contour(double f[][100],double min,double max,double df){
  FILE *fid4;
  const char *data4={"contour.dat"};
  fid4 = fopen(data4,"w");
  int num = ((max-min)/df) + 1; //等高線の数
  vector< vector <Point> > contours;
  contours = vector< vector <Point> >(num,vector <Point>(0));
  for (int k=0;k < num;k++){
    double f_num = min + k * df;
    double threshold=0.02;
    for(int i=0;i<181;i++){
      //流線の時の例外措置
      //if(k == 15) break;
      //if(k == 16) threshold = 0.02;
      for(int j=0;j<100;j++){
        if(f[i][j] > (f_num - threshold*df) && f[i][j] < (f_num + threshold*df) ){
          Point point_ij;
          point_ij.x = x_zeta[i][j];
          point_ij.y = y_zeta[i][j];
          contours[k].push_back(point_ij);
        }
      }
    }
  }
  //fprintf(gid,"set multiplot\n");
  for(int n=0;n<num;n++){
    double n_num = min + n * df;
    cout << "number:" << n_num << " contour size:" << contours[n].size() << "\n";
    //plot
    for (int m = 0; m < contours[n].size(); m++) {
        fprintf(fid4, "%f\t%f\n", contours[n][m].x, contours[n][m].y);
      }
    fprintf(fid4, "\n");
    }
  fprintf(fid4,"\n");
    //fprintf(gid,"unset multiplot\n");
}*/

int main(){
  int mode;
  double alfa,beta;
  printf("choose mode:(0=wing_shape,1=stream,2=pressure,3=surface_pressure)\n");
  scanf("%d",&mode);
  printf("enter alfa;\n");
  scanf("%d",&alfa);
  printf("enter beta:\n");
  scanf("%lf",&beta);
  alfa = alfa * (PI/180.0);
  beta = beta * (PI/180.0);

  calculation(alfa,beta);
  show_data(mode);
  cal_force(alfa,beta);
  cal_ac(beta);
  close_gnuplot();
}
