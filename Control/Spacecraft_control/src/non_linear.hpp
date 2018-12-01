// attitude for spinning satelite
#ifndef NONLINEAR_H
#define NONLINEAR_H

#include "Eigen/Core"
#include "Eigen/Geometry"
#include <cstdio>
#include <string>

class NonLinear {
public:
  NonLinear(Eigen::Vector3d &I, double SD_Q, double SD_R, int SEED) {
    Ix_ = I[0];
    Iy_ = I[1];
    Iz_ = I[2];
    SD_Q_ = SD_Q;
    SD_R_ = SD_R;
    SEED_ = SEED;
  };

  double Ix_, Iy_, Iz_;
  Eigen::Vector3d vV_; //観測ノイズ
  Eigen::Vector3d vW_; //外乱ノイズ

  //  functions for state propagation
  void propagation(double dt, Eigen::VectorXd &omega_qt,
                   Eigen::Matrix<double, 3, 3> &dcm, double t_start,
                   double T_END);
  void cal_omega_qt_dot(Eigen::VectorXd &omega_qt,
                        Eigen::VectorXd &omega_qt_dot);
  void quartanion_to_dcm(Eigen::VectorXd &omega_qt,
                         Eigen::Matrix<double, 3, 3> &dcm_);

  // functions for producing/ deleting noise
  void update_vW(); // disterbuance noise
  void update_vV(); // observation noise
  void del_noise();

  // writing result in files
  void writefile(double t, Eigen::VectorXd &omega_qt,
                 Eigen::Matrix<double, 3, 3> &dcm, std::string fileprefix);

private:
  double SD_Q_, SD_R_; /*ノイズ標準偏差 SD_Q_:外乱トルク SD_R:観測ノイズ*/
  double SEED_; /*メルセンヌツイスタ関数用シード 0-シードも乱数で決定*/
};

#endif
