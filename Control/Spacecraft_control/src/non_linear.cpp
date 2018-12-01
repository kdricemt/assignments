#include "non_linear.hpp"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <string>

/***********************************************************
    state propgation function
      propagation()
      cal_omega_qt_dot()
      quartanion_to_dcm()
*************************************************************/

/* propagation
   @Input
     dt: 時間ステップ
     omega_qt: 伝搬するstate(角速度ω, qt)
     dcm: qtに対応するdcm
     t_start: 積分開始時間
     t_end: 積分終了時間
*/
void NonLinear::propagation(double dt, Eigen::VectorXd &omega_qt,
                            Eigen::Matrix<double, 3, 3> &dcm, double t_start,
                            double T_END) {
  for (double t = t_start; t < T_END; t += dt) { // update
    Eigen::VectorXd k1, k2, k3, k4;
    Eigen::VectorXd omega_qt_2, omega_qt_3, omega_qt_4;
    k1 = Eigen::VectorXd::Zero(7);
    k2 = Eigen::VectorXd::Zero(7);
    k3 = Eigen::VectorXd::Zero(7);
    k4 = Eigen::VectorXd::Zero(7);
    omega_qt_2 = Eigen::VectorXd::Zero(7);
    omega_qt_3 = Eigen::VectorXd::Zero(7);
    omega_qt_4 = Eigen::VectorXd::Zero(7);
    cal_omega_qt_dot(omega_qt, k1);

    omega_qt_2 = omega_qt + dt / 2 * k1;
    cal_omega_qt_dot(omega_qt_2, k2);

    omega_qt_3 = omega_qt + dt / 2 * k2;
    cal_omega_qt_dot(omega_qt_3, k3);

    omega_qt_4 = omega_qt + dt * k3;
    cal_omega_qt_dot(omega_qt_4, k4);

    omega_qt = omega_qt + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    quartanion_to_dcm(omega_qt, dcm);
  }
}

// ω, qtの微分を求める
void NonLinear::cal_omega_qt_dot(Eigen::VectorXd &omega_qt,
                                 Eigen::VectorXd &omega_qt_dot) {
  double q0 = omega_qt(0, 0);
  double q1 = omega_qt(1, 0);
  double q2 = omega_qt(2, 0);
  double q3 = omega_qt(3, 0);
  double ox = omega_qt(4, 0);
  double oy = omega_qt(5, 0);
  double oz = omega_qt(6, 0);

  omega_qt_dot(0, 0) = 0.5 * (-q1 * ox + -q2 * oy + -q3 * oz);
  omega_qt_dot(1, 0) = 0.5 * (q0 * ox + -q3 * oy + q2 * oz);
  omega_qt_dot(2, 0) = 0.5 * (q3 * ox + q0 * oy + -q1 * oz);
  omega_qt_dot(3, 0) = 0.5 * (-q2 * ox + q1 * oy + q0 * oz);
  omega_qt_dot(4, 0) = (Iy_ - Iz_) / Ix_ * oy * oz + vW_[0] / Ix_;
  omega_qt_dot(5, 0) = (Iz_ - Ix_) / Iy_ * oz * ox + vW_[1] / Iy_;
  omega_qt_dot(6, 0) = (Ix_ - Iy_) / Iz_ * ox * oy + vW_[2] / Iz_;
}

// qt->dcm
void NonLinear::quartanion_to_dcm(Eigen::VectorXd &omega_qt,
                                  Eigen::Matrix<double, 3, 3> &dcm) {
  double q0 = omega_qt(0, 0);
  double q1 = omega_qt(1, 0);
  double q2 = omega_qt(2, 0);
  double q3 = omega_qt(3, 0);
  // first column
  dcm(0, 0) = q0 * q0 + q1 * q1 - q2 * q2 - q3 * q3;
  dcm(1, 0) = 2 * (q1 * q3 + q0 * q3);
  dcm(2, 0) = 2 * (q1 * q3 - q0 * q2);
  // second column
  dcm(0, 1) = 2 * (q1 * q2 - q0 * q3);
  dcm(1, 1) = q0 * q0 - q1 * q1 + q2 * q2 - q3 * q3;
  dcm(2, 1) = 2 * (q2 * q3 + q0 * q1);
  // third column
  dcm(0, 2) = 2 * (q1 * q3 + q0 * q2);
  dcm(1, 2) = 2 * (q2 * q3 - q0 * q1);
  dcm(2, 2) = q0 * q0 - q1 * q1 - q2 * q2 + q3 * q3;
}

/********************************************************
    noise producing function
      update_vW() : 外乱ノイズ
      update_vV() : 観測ノイズ
      del_noise() : ノイズ消去
***********************************************************/

// 外乱ノイズ生成関数 seed=0の場合, seedもランダム生成
void NonLinear::update_vW() {
  if (SEED_ == 0) {
    std::random_device seed_gen;
    std::mt19937 randMt(seed_gen());
    static std::normal_distribution<> normW(0.0, SD_Q_);
    vW_(0) = normW(randMt);
    vW_(1) = normW(randMt);
    vW_(2) = normW(randMt);
  } else {
    std::mt19937 randMt(SEED_); //メルセンヌツイスタ乱数
    static std::normal_distribution<> normW(0.0, SD_Q_);
    vW_(0) = normW(randMt);
    vW_(1) = normW(randMt);
    vW_(2) = normW(randMt);
  }
}

// 観測ノイズ生成関数, seed=0の時はseedもランダム生成
void NonLinear::update_vV() {
  if (SEED_ == 0) {
    std::random_device seed_gen;
    std::mt19937 randMt(seed_gen());
    static std::normal_distribution<> normW(0.0, SD_Q_);
    vV_(0) = normW(randMt);
    vV_(1) = normW(randMt);
    vV_(2) = normW(randMt);
  } else {
    std::mt19937 randMt(SEED_); //メルセンヌツイスタ乱数
    static std::normal_distribution<> normW(0.0, SD_Q_);
    vV_(0) = normW(randMt);
    vV_(1) = normW(randMt);
    vV_(2) = normW(randMt);
  }
}

// システムノイズ/観測ノイズ消去
void NonLinear::del_noise() {
  vV_(0) = 0;
  vV_(1) = 0;
  vV_(2) = 0;
  vW_(0) = 0;
  vW_(1) = 0;
  vW_(2) = 0;
}

/********************************************************
    Writing quartanion,omega,dcm data
      writefile()
*********************************************************/

/* @input
     t: 時間
     omega_qt,dcm: state
     fileprefix: ファイル先頭につける文字
*/

void NonLinear::writefile(double t, Eigen::VectorXd &omega_qt,
                          Eigen::Matrix<double, 3, 3> &dcm,
                          std::string fileprefix) {
  // write data of quartanion
  std::fstream fout_q;
  std::string filename_q;
  if (SD_Q_ > 0) {
    filename_q = "../datafile/" + fileprefix + "qt_error.txt";
  } else {
    filename_q = "../datafile/" + fileprefix + "qt.txt";
  }
  if (t == 0)
    fout_q.open(filename_q, std::ios::out);
  else
    fout_q.open(filename_q, std::ios::app);
  fout_q << t << " ";
  for (int n = 0; n < 4; n++) {
    fout_q << omega_qt(n, 0) << " ";
    if (n == 3) {
      fout_q << std::endl;
    }
  }
  // write data of omega
  std::fstream fout_o;
  std::string filename_o;
  if (SD_Q_ > 0) {
    filename_o = "../datafile/" + fileprefix + "omega_error.txt";
  } else {
    filename_o = "../datafile/" + fileprefix + "omega.txt";
  }
  if (t == 0)
    fout_o.open(filename_o, std::ios::out);
  else
    fout_o.open(filename_o, std::ios::app);
  fout_o << t << " ";
  for (int m = 0; m < 3; m++) {
    fout_o << omega_qt(m + 4, 0) << " ";
    if (m == 2) {
      fout_o << std::endl;
    }
  }
  // write data for dcm
  for (int m = 0; m < 3; m++) {
    std::fstream fout_d;
    std::string filename_d;
    if (SD_Q_ > 0) {
      filename_d = "../datafile/" + fileprefix + "dcm_column" +
                   std::to_string(m) + "_error" + ".txt";
    } else {
      filename_d = "../datafile/" + fileprefix + "dcm_column" +
                   std::to_string(m) + ".txt";
    }
    if (t == 0)
      fout_d.open(filename_d, std::ios::out);
    else
      fout_d.open(filename_d, std::ios::app);
    fout_d << t << " " << dcm(0, m) << " " << dcm(1, m) << " " << dcm(2, m)
           << std::endl;
  }
}
