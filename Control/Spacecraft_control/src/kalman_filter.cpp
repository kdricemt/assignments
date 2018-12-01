#include "Eigen/LU"
#include "unsupported/Eigen/MatrixFunctions"
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>

#include "kalman_filter.hpp"

/***********************************************************
    Update Functions for A,B,H
      void UpdateMatrixA(NonLinear &nol)
      void UpdateMatrixB(NonLinear &nol)
      void UpdateMatrixH(int col)
**********************************************************/

void KalmanFilter::UpdateMatrixA(NonLinear &nol) {
  double Ix = nol.Ix_;
  double Iy = nol.Iy_;
  double Iz = nol.Iz_;
  double q0 = vSysEst(0);
  double q1 = vSysEst(1);
  double q2 = vSysEst(2);
  double q3 = vSysEst(3);
  double wx = vSysEst(4);
  double wy = vSysEst(5);
  double wz = vSysEst(6);
  mA(0, 0) = 0.0;
  mA(0, 1) = -0.5 * wx;
  mA(0, 2) = -0.5 * wy;
  mA(0, 3) = -0.5 * wz;

  mA(1, 0) = +0.5 * wx;
  mA(1, 1) = 0.0;
  mA(1, 2) = +0.5 * wz;
  mA(1, 3) = -0.5 * wy;

  mA(2, 0) = +0.5 * wy;
  mA(2, 1) = -0.5 * wz;
  mA(2, 2) = 0.0;
  mA(2, 3) = +0.5 * wx;

  mA(3, 0) = +0.5 * wz;
  mA(3, 1) = +0.5 * wy;
  mA(3, 2) = -0.5 * wx;
  mA(3, 3) = 0.0;

  mA(4, 0) = 0.0;
  mA(4, 1) = 0.0;
  mA(4, 2) = 0.0;
  mA(4, 3) = 0.0;

  mA(5, 0) = 0.0;
  mA(5, 1) = 0.0;
  mA(5, 2) = 0.0;
  mA(5, 3) = 0.0;

  mA(6, 0) = 0.0;
  mA(6, 1) = 0.0;
  mA(6, 2) = 0.0;
  mA(6, 3) = 0.0;

  mA(0, 4) = -0.5 * q1;
  mA(0, 5) = -0.5 * q2;
  mA(0, 6) = -0.5 * q3;

  mA(1, 4) = +0.5 * q0;
  mA(1, 5) = -0.5 * q3;
  mA(1, 6) = +0.5 * q2;

  mA(2, 4) = +0.5 * q3;
  mA(2, 5) = +0.5 * q0;
  mA(2, 6) = -0.5 * q1;

  mA(3, 4) = -0.5 * q2;
  mA(3, 5) = +0.5 * q1;
  mA(3, 6) = +0.5 * q0;

  mA(4, 4) = 0.0;
  mA(4, 5) = wz * (Iy - Iz) / Ix;
  mA(4, 6) = wy * (Iy - Iz) / Ix;

  mA(5, 4) = wz * (Iz - Ix) / Iy;
  mA(5, 5) = 0.0;
  mA(5, 6) = wx * (Iz - Ix) / Iy;

  mA(6, 4) = wy * (Ix - Iy) / Iz;
  mA(6, 5) = wx * (Ix - Iy) / Iz;
  mA(6, 6) = 0.0;
}

void KalmanFilter::UpdateMatrixB(NonLinear &nol) {
  mB(4, 0) = 1.0 / nol.Ix_;
  mB(5, 1) = 1.0 / nol.Iy_;
  mB(6, 2) = 1.0 / nol.Iz_;
}

/*
    @Input: col: column of DCM to observe
*/

void KalmanFilter::UpdateMatrixH(int col) {
  double q0 = vSysEst(0) * 2.0;
  double q1 = vSysEst(1) * 2.0;
  double q2 = vSysEst(2) * 2.0;
  double q3 = vSysEst(3) * 2.0;

  for (int i = 0; i < 3; i++) {
    mH(i, 4) = 0.0;
    mH(i, 5) = 0.0;
    mH(i, 6) = 0.0;
  }
  if (col == 0) {
    mH(0, 0) = +q0;
    mH(0, 1) = +q1;
    mH(0, 2) = -q2;
    mH(0, 3) = -q3;
    mH(1, 0) = +q3;
    mH(1, 1) = +q2;
    mH(1, 2) = +q1;
    mH(1, 3) = +q0;
    mH(2, 0) = -q2;
    mH(2, 1) = +q3;
    mH(2, 2) = -q0;
    mH(2, 3) = +q1;
  } else if (col == 1) {
    mH(0, 0) = -q3;
    mH(0, 1) = +q2;
    mH(0, 2) = +q1;
    mH(0, 3) = -q0;
    mH(1, 0) = +q0;
    mH(1, 1) = -q1;
    mH(1, 2) = +q2;
    mH(1, 3) = -q3;
    mH(2, 0) = +q1;
    mH(2, 1) = +q0;
    mH(2, 2) = +q3;
    mH(2, 3) = +q2;
  } else if (col == 2) {
    mH(0, 0) = +q2;
    mH(0, 1) = +q3;
    mH(0, 2) = +q0;
    mH(0, 3) = +q1;
    mH(1, 0) = -q1;
    mH(1, 1) = -q0;
    mH(1, 2) = +q3;
    mH(1, 3) = +q2;
    mH(2, 0) = +q0;
    mH(2, 1) = -q1;
    mH(2, 2) = -q2;
    mH(2, 3) = +q3;
  } else {
    std::cout << "DCM COL INDEX ERROR" << std::endl;
  }
}

/******************************************************************************
    Update Kalman-filter Matrices
      void UpdateMatrixPhi()
      void UpdateMatrixGmm()
*****************************************************************************/

void KalmanFilter::UpdateMatrixPhi() {
  Eigen::MatrixXd mTemp = dt * mA;
  mPhi = mTemp.exp();
}

void KalmanFilter::UpdateMatrixGmm() {
  static Eigen::MatrixXd mI = Eigen::MatrixXd::Identity(dimS, dimS);
  Eigen::MatrixXd mTemp = dt * mA;
  mGmm = mA.inverse() * (mTemp.exp() - mI) * mB;
}

/**********************************************************
    Fuctions for simulation
      void KalmanFilter::simulation(NonLinear &nol, double T_END,
                                  double obs_interval)
      double KalmanFilter::CalcQuaternionNorm(Eigen::VectorXd &vX)
      void KalmanFilter::NormalizeQuaternion(Eigen::VectorXd &vX)
      double KalmanFilter::CalcErrNorm(Eigen::VectorXd &vSysTru,
                                       Eigen::VectorXd &vSysEst)
      int KalmanFilter::select_column(int seed)
*********************************************************/

/*
    @input
      nol: NonLinear Class for non-linear model propagation
      T_END: simulation time
      obs_interval: observation interval
*/

void KalmanFilter::simulation(NonLinear &nol, double T_END,
                              double obs_interval) {
  Eigen::Matrix<double, 3, 3> dcm_Tru;
  Eigen::Matrix<double, 3, 3> dcm_Est;
  Eigen::VectorXd vY_Tru; // one column of dcm_Tru + error
  Eigen::VectorXd vY_Est; // one column of

  vY_Tru = Eigen::VectorXd::Zero(dimO);
  vY_Est = Eigen::VectorXd::Zero(dimO);
  double time_obs = 0.0; // counter

  for (double t = 0; t < T_END + dt; t += dt) {
    nol.writefile(t, vSysTru, dcm_Tru, "R2True");
    nol.writefile(t, vSysEst, dcm_Est, "R2Est");
    // 真値計算
    nol.update_vW(); // make disturbance noise
    nol.propagation(dt, vSysTru, dcm_Tru, t, t + dt);
    // 推定値計算
    nol.del_noise();
    nol.propagation(dt, vSysEst, dcm_Est, t, t + dt); // estimate with no noise

    // Matrix更新
    UpdateMatrixA(nol);
    UpdateMatrixB(nol);
    UpdateMatrixPhi();
    UpdateMatrixGmm();
    mM = mPhi * mP * mPhi.transpose() + mGmm * mQ * mGmm.transpose();
    mP = mM;
    if (time_obs >= obs_interval) {
      std::cout << "observe time: " << t << std::endl;
      std::cout << "  vSysTru:" << vSysTru.transpose() << std::endl;
      std::cout << "  vSysEst:" << vSysEst.transpose() << std::endl;
      time_obs = 0.0;
      int seed = 0;

      int dcm_col = select_column(seed); // select column to observe
      std::cout << "  column observed:" << dcm_col << std::endl;
      UpdateMatrixH(dcm_col);
      nol.update_vV();                         // make observation noise
      vY_Tru = dcm_Tru.col(dcm_col) + nol.vV_; // real
      nol.del_noise();
      vY_Est = dcm_Est.col(dcm_col) + nol.vV_; // estimation
      vz = vY_Tru - vY_Est;
      std::cout << "  vz:" << vz.transpose() << std::endl;
      double err = CalcErrNorm(vSysTru, vSysEst);
      std::cout << "  Error Norm(before):" << err << std::endl;

      Eigen::MatrixXd mTemp; // temporary defined matrix used for calculation
      mTemp = Eigen::MatrixXd::Zero(3, 3);
      mTemp = mH * mM * mH.transpose() + mR;
      mP = mM - mM * mH.transpose() * mTemp.inverse() * mH * mM;
      mK = mP * mH.transpose() * mR.inverse();

      vSysEst = vSysEst + mK * vz; // update estimation
      err = CalcErrNorm(vSysTru, vSysEst);
      std::cout << "  Error Norm(after):" << err << std::endl;
      std::cout << "  Kalman Gain norm:" << mK.norm() << std::endl;
      for (int i = 0; i < dimS; i++) {
        std::cout << "  P" << i << ":" << sqrt(mP(i, i)) << " ";
        if (i == dimS - 1)
          std::cout << std::endl;
      }
      std::cout << std::endl;
      NormalizeQuaternion(vSysEst);
      write_file_k(t);
    }
    write_file_pe(t); // write mP,mK
    time_obs += dt;
  }
}

/* @Input
    vX: System states Vector
*/
double KalmanFilter::CalcQuaternionNorm(Eigen::VectorXd &vX) {
  double result = 0;
  for (int i = 0; i < 4; i++) {
    result += vX(i) * vX(i);
  }
  return sqrt(result);
}

/* @Input
    vX: System states Vector
*/
void KalmanFilter::NormalizeQuaternion(Eigen::VectorXd &vX) {
  double norm = CalcQuaternionNorm(vX);
  for (int i = 0; i < 4; i++) {
    vX(i) = vX(i) / norm;
  }
}

/* @Input
    vSysTru: System states Vector (True(propagated with errors))
    vSysEst: System states Vector (Estimated(propagated without errors))
*/
double KalmanFilter::CalcErrNorm(Eigen::VectorXd &vSysTru,
                                 Eigen::VectorXd &vSysEst) {
  Eigen::VectorXd vD = vSysTru - vSysEst;
  return vD.norm();
}

// select column od DCM (for observation(updating Matrix H))
int KalmanFilter::select_column(int seed) {
  int column = 0;
  if (seed == 0) {
    std::random_device seed_gen;
    std::mt19937 randMt(seed_gen());
    std::uniform_int_distribution<int> rand_dist(0, 2);
    column = rand_dist(randMt);
  } else {
    std::mt19937 randMt(seed);
    std::uniform_int_distribution<int> rand_dist(0, 2);
    column = rand_dist(randMt);
  }
  return column;
}

/*****************************************************************************
    Functions for writing result (Matrix P,K)
      write_file_pe()
      write_file_k()
*****************************************************************************/

void KalmanFilter::write_file_pe(double t) {
  // files to write result
  // Estimate error covariance
  std::fstream fout_p;
  std::string filename_p;
  filename_p = "../datafile/R2p.txt";
  // errors
  std::fstream fout_e;
  std::string filename_e;
  filename_e = "../datafile/R2e.txt";

  if (t < 2 * dt) {
    fout_p.open(filename_p, std::ios::out);
    fout_e.open(filename_e, std::ios::out);

  } else {
    fout_p.open(filename_p, std::ios::app);
    fout_e.open(filename_e, std::ios::app);
  }

  fout_p << t << " ";
  fout_e << t << " ";
  for (int i = 0; i < dimS; i++) {
    fout_p << sqrt(mP(i, i)) << " ";
    if (i == dimS - 1)
      fout_p << std::endl;
  }
  Eigen::VectorXd vD = vSysTru - vSysEst;
  for (int i = 0; i < dimS; i++) {
    fout_e << vD(i) << " ";
    if (i == dimS - 1)
      fout_e << std::endl;
  }
}

void KalmanFilter::write_file_k(double t) {
  // files to write result
  // kalman gain
  std::fstream fout_k;
  std::string filename_k;
  filename_k = "../datafile/R2k.txt";

  if (t < 2) {
    fout_k.open(filename_k, std::ios::out);

  } else {
    fout_k.open(filename_k, std::ios::app);
  }
  double norm = mK.norm();
  fout_k << t << " " << norm << std::endl;
}
