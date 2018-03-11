#include "kalman_filter.h"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

// Please note that the Eigen library does not initialize 
// VectorXd or MatrixXd objects with zeros upon creation.

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in) {
  x_ = x_in;
  P_ = P_in;
}

void KalmanFilter::Predict() {

  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  CalculateGainAndMakeEstimation(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // Normalize input data
  double in_rho = z[0];
  double in_theta = z[1];
  double in_drho = z[2];
  while (in_theta > M_PI) {
    in_theta -= 2 * M_PI;
  }
  while (in_theta < -M_PI) {
    in_theta += 2 * M_PI;
  }
  VectorXd z_in = VectorXd(3);
  z_in << in_rho, in_theta, in_drho;

  // Convert state vector to radar coordinate system
  double px = x_(0);
  double py = x_(1);
  double pvx = x_(2);
  double pvy = x_(3);

  double rho = sqrt(px * px + py * py);
  double theta = atan2(py, px);
  double drho = (px * pvx + py * pvy) / sqrt(px * px + py * py);
  VectorXd z_pred = VectorXd(3);
  z_pred << rho, theta, drho;

  VectorXd y = z_in - z_pred;

  CalculateGainAndMakeEstimation(y);
}

void KalmanFilter::CalculateGainAndMakeEstimation(const VectorXd &y) {

  // Calculate Kalman gain
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd PHt = P_ * Ht;
  MatrixXd K = PHt * Si;

  //new estimate
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}
