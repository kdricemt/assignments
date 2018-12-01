#include "Eigen/Core"
#include "Eigen/Geometry"
#include "iostream"
#include "kalman_filter.hpp"
#include "non_linear.hpp"

int main() {

  // setup for system model

  // settings for satelite attitude
  Eigen::Vector3d I;
  I.x() = 1.9; // I_x
  I.y() = 1.6; // I_y
  I.z() = 2.0; // I_z

  Eigen::Vector3d omega_init; //初期角速度
  Eigen::Vector3d omega_init_error;
  omega_init << 0.1, 17 * 0.1047 + 0.1, 0;
  omega_init_error << 0.01, 0.01, 0.01;

  Eigen::Vector4d qt_init; //初期Quartanion
  Eigen::Vector4d qt_init_errror;
  qt_init << 1.0, 0, 0, 0;
  qt_init_errror << -0.01, 0.08, 0.08, 0.08;

  // initial state
  Eigen::VectorXd vX0_Tru;
  Eigen::VectorXd vX0_Est;
  vX0_Tru = Eigen::VectorXd::Zero(7);
  vX0_Est = Eigen::VectorXd::Zero(7);

  vX0_Tru << qt_init[0], qt_init[1], qt_init[2], qt_init[3], omega_init[0],
      omega_init[1], omega_init[2];
  omega_init += omega_init_error;
  qt_init += qt_init_errror;
  vX0_Est << qt_init[0], qt_init[1], qt_init[2], qt_init[3], omega_init[0],
      omega_init[1], omega_init[2];

  double SD_Q = 0.01; //外乱トルクの標準偏差
  double SD_R = 0.01; //観測ノイズの標準偏差
  double SEED =
      0; //外乱の乱数生成メルセンヌツイスタ用シード,0の時シードもランダムに生成される

  NonLinear nol(I, SD_Q, SD_R, SEED);

  // System Parameters
  int dimI = 3;
  int dimO = 3;
  int dimS = 7;

  // Define Initial Vectors and Matrix
  Eigen::VectorXd vz;
  Eigen::MatrixXd Q;
  Eigen::MatrixXd R;
  Eigen::MatrixXd P0;

  // init Vectors and Matrix
  vz = Eigen::VectorXd::Zero(dimO);
  Q = Eigen::MatrixXd::Identity(dimI, dimI);
  R = Eigen::MatrixXd::Identity(dimO, dimO);
  P0 = Eigen::MatrixXd::Identity(dimS, dimS);

  // Q process noise
  Q(0, 0) = SD_Q * SD_Q;
  Q(1, 1) = SD_Q * SD_Q;
  Q(2, 2) = SD_Q * SD_Q;
  // R measurement noise
  R(0, 0) = SD_R * SD_R;
  R(1, 1) = SD_R * SD_R;
  R(2, 2) = SD_R * SD_R;
  // P0
  for (int i = 0; i < dimS; i++) {
    P0(i, i) = SD_R;
  }

  // simulation用パラメータ
  double dt = 0.01;
  double T_END = 100;
  double obs_interval = 1.0;

  // end model setup

  // make model
  KalmanFilter KFD(dimI, dimO, dimS, vX0_Tru, vX0_Est, Q, R, P0, dt);

  // simulate
  KFD.simulation(nol, T_END, obs_interval);
}
