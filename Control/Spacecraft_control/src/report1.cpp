#include "Eigen/Core"
#include "Eigen/Geometry"
#include "non_linear.hpp"
#include <iostream>

int main() {
  Eigen::Vector3d I;
  I.x() = 1.9; // I_x
  I.y() = 1.6; // I_y
  I.z() = 2.0; // I_z

  Eigen::Vector3d omega_init; //ノミナル角速度
  omega_init.x() = 0.1;
  omega_init.y() = 17 * 0.1047 + 0.1;
  omega_init.z() = 0;

  Eigen::Vector4d qt_init; //初期Quartanion
  qt_init << 1.0, 0, 0, 0;

  // Initialize State
  Eigen::VectorXd omega_qt;
  Eigen::Matrix<double, 3, 3> dcm;
  omega_qt = Eigen::VectorXd::Zero(7);

  omega_qt << qt_init[0], qt_init[1], qt_init[2], qt_init[3], omega_init[0],
      omega_init[1], omega_init[2];

  double SD_Q = 0.0; //外乱トルクの標準偏差
  double SD_R = 0.0; //観測ノイズの標準偏差
  double SEED =
      1; //外乱の乱数生成メルセンヌツイスタ用シード,0の時シードもランダムに生成される

  NonLinear nol(I, SD_Q, SD_R, SEED);

  //ルンゲクッタ用パラメータ
  double dt = 0.01;
  double T_END = 100;

  int counter = 0; // print用カウンター

  for (double t = 0; t < T_END; t += dt) {
    if (!(counter % 100)) {
      std::cout << "time:" << t << std::endl;
    }
    nol.update_vW(); // make disterbuance noise
    nol.propagation(dt, omega_qt, dcm, t, t + dt);
    nol.writefile(t, omega_qt, dcm, "R1");
    counter += 1;
  }
}
