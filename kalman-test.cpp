/**
 *  @file   ExampleAlgorithms/KalmanFilterAlgorithm.cc
 *
 *  @brief  Implementation of the list preparation algorithm class.
 *
 *  $Log: $
 */

#include "kalman-test.h"
#include <fstream>

//------------------------------------------------------------------------------------------------------------------------------------------


KalmanFilterAlgorithm::KalmanFilterAlgorithm(
    double kdt,
    const Eigen::MatrixXd& kA,
    const Eigen::MatrixXd& kC,
    const Eigen::MatrixXd& kQ,
    const Eigen::MatrixXd& kR,
    const Eigen::MatrixXd& kP)
  : kA(kA), kC(kC), kQ(kQ), kR(kR), kP0(kP),
    m(kC.rows()), n(kA.rows()), kdt(kdt), kinitialized(false),
    kI(n, n), kx_hat(n), kx_hat_new(n)
{
  kI.setIdentity();
}


void KalmanFilterAlgorithm::init(double t0, const Eigen::VectorXd& x0) {
  kx_hat = x0;
  kP = kP0;
  this->kt0 = t0;
  kt = t0;
  kinitialized = true;
}


void KalmanFilterAlgorithm::update(const Eigen::VectorXd& y) {

  if(!kinitialized)
    throw std::runtime_error("Filter is not initialized!");

  kx_hat_new = kA * kx_hat;
  kP = kA*kP*kA.transpose() + kQ;
  kK = kP*kC.transpose()*(kC*kP*kC.transpose() + kR).inverse();
  kx_hat_new += kK * (y - kC*kx_hat_new);
  kP = (kI - kK*kC)*kP;
  kx_hat = kx_hat_new;

  kt += kdt;
}

void KalmanFilterAlgorithm::update(const Eigen::VectorXd& y, double dt, const Eigen::MatrixXd A) {

  this->kA = A;
  this->kdt = dt;
  update(y);
}

