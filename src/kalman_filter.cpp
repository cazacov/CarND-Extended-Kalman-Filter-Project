#include "kalman_filter.h"

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
    auto Ft = F_.transpose();
    P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

    // calculate Kalman gain
    auto z_pred = H_ * x_;
    auto y = z - z_pred;
    auto Ht = H_.transpose();
    auto S = H_ * P_ * Ht + R_;
    auto Si = S.inverse();
    auto PHt = P_ * Ht;
    auto K = PHt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
}
