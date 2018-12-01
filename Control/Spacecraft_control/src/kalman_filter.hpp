/*********************************************************
 * Descreted Kalman filter Program
 *
 * < Linearlized Equation of Stateo>
 *   Δx' = mA * Δx  + mB * (u+vW)
 *   y = mH * Δx + vV
 *
 * <Descreted Linearlized Equation of State>
 *   Δx[k]' = mPhi * Δx[k-1] + mGmm * vW[k-1]
 *
 * <Matrices>
 *   mA - System dynamics matrix(State Matrix) (SxS)
 *   mB - input matrix (SxI)
 *   mH - Output matrix (OxS)
 *   mQ - Process noise covariance (IxI)
 *   mR - Measurement noise covariance (OxO)
 *   mP - Estimate error covariance (SxS)
 *   mK - Kalman Gain
 *
 * <Vectors>
 *   vSysTru - True State Vector
 *   vSysEst - State Vector for Estimated Models
 *   vW - disturbance noise(defined in NonLinear Class)
 *   vV - Observation noise(defined in NonLinear Class)
 *
 *  Non-Linear equation of state is written in class NonLinear.
 *  For propagation of state vectors in both True State and
 *  Estimated Model,propagation function in class NonLinear
 *  is used. In addition,system parameters in class NonLinear
 *  might be used to update A,B,H matrices.
 *
 **************************************************************/

#ifndef KFD_H
#define KFD_H
#include "Eigen/Core"
#include "Eigen/Geometry"
#include "non_linear.hpp"
#include <cstdio>
#include <iostream>

class KalmanFilter {

public:
  KalmanFilter(int dimI, int dimO, int dimS, const Eigen::VectorXd &vX0_Tru,
               const Eigen::VectorXd &vX0_EST, const Eigen::MatrixXd &Q,
               const Eigen::MatrixXd &R, const Eigen::MatrixXd &P0,
               const double dt)
      : dimI(dimI), dimO(dimO), dimS(dimS), vSysTru(vX0_Tru), vSysEst(vX0_EST),
        mQ(Q), mR(R), mP(P0), dt(dt) {
    // init Vectors and Matrix
    mA = Eigen::MatrixXd::Zero(dimS, dimS);
    mB = Eigen::MatrixXd::Zero(dimS, dimI);
    mH = Eigen::MatrixXd::Zero(dimO, dimS);
    vz = Eigen::VectorXd::Zero(dimO);

    // init
    mM = mP;
    mGmm = Eigen::MatrixXd::Zero(dimS, dimI);
    mK = Eigen::MatrixXd::Zero(dimS, dimO);
    // display Initial Matrix
    std::cout << "A" << std::endl << mA << std ::endl;
    std::cout << "B" << std::endl << mB << std ::endl;
    std::cout << "H" << std::endl << mH << std ::endl;
    std::cout << "Q" << std::endl << mQ << std ::endl;
    std::cout << "R" << std::endl << mR << std ::endl;
    std::cout << "P" << std::endl << mP << std ::endl;
    std::cout << "M" << std::endl << mM << std ::endl;
    std::cout << "K" << std::endl << mK << std ::endl;
    std::cout << "H" << std::endl << mH << std ::endl;
    std::cout << "Phi" << std::endl << mPhi << std ::endl;
    std::cout << "Gmm" << std::endl << mGmm << std ::endl;
  };

  // System Parameters (class SpinSat)
  // System Matrices
  int dimI, dimO, dimS; // Number of Input,Output,State
  Eigen::MatrixXd mA;   // State Matrix
  Eigen::MatrixXd mB;   // Input Matrix
  Eigen::MatrixXd mH;   // Output Matrix
  // Matrices for Descreted State Equation
  Eigen::MatrixXd mPhi, mGmm;
  // Vectors for states(vSys = x)
  Eigen::VectorXd vSysTru, vSysEst;
  Eigen::VectorXd vz; // observed Y - estimated Y
  // Matrices for errors
  Eigen::MatrixXd mQ, mR;
  Eigen::MatrixXd mP, mM, mK;
  // Discrete time step
  double dt;

  //@ method

  /* Update functions (need to rewrite when system is changed)*/
  void UpdateMatrixA(NonLinear &nol);
  void UpdateMatrixB(NonLinear &nol);
  void UpdateMatrixH(int col);
  // Calculate Phi,Gmm using Matrix A and Matrix B
  void UpdateMatrixPhi();
  void UpdateMatrixGmm();
  // Calculating State Vector Norms
  double CalcQuaternionNorm(Eigen::VectorXd &vX);
  void NormalizeQuaternion(Eigen::VectorXd &vX);
  double CalcErrNorm(Eigen::VectorXd &vSysTru, Eigen::VectorXd &vSysEst);

  // select column of DCM observed
  int select_column(int seed);

  /* Whole Update Sequence(need to rewrite when system is changed)*/
  void simulation(NonLinear &nol, double T_END, double obs_interval);

  // write files
  void write_file_pe(double t);
  void write_file_k(double t);
};

#endif
