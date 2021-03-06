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

  // Apply state prediction function
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) {

  // Apply measurement function
  VectorXd z_pred = H_ * x_;
  VectorXd y = z - z_pred;

  CalculateGainAndMakeEstimation(y);
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {

  // Convert state vector to radar coordinate system
  double px = x_(0);
  double py = x_(1);
  double pvx = x_(2);
  double pvy = x_(3);

  // Apply measurement function
  double p_rho = sqrt(px * px + py * py);
  double p_theta = atan2(py, px);
  double p_rho_dot = (px * pvx + py * pvy) / sqrt(px * px + py * py);
  VectorXd z_pred = VectorXd(3);
  z_pred << p_rho, p_theta, p_rho_dot;


  double in_rho = z[0];
  double in_theta = z[1];
  double in_rho_dot = z[2];

  // Put input angle into range -PI..+PI
  while (in_theta > M_PI) {
    in_theta -= 2 * M_PI;
  }
  while (in_theta < -M_PI) {
    in_theta += 2 * M_PI;
  }

  /*****************************************************************************
   *  Make algorithm more stable if theta is near +PI or -PI
   ****************************************************************************/

  if (fabs(in_theta) > M_PI / 2 && fabs(p_theta) > M_PI / 2) {
    /* Both measured and predicted angles are in the rear semi-plane. That can be dangerous!
     *
     * Let's say predicted value is +178 degrees and the radar reports -170 degrees.
     * That makes only 12 degrees difference, but the algorithm will try to incrementally
     * * update the state vector moving it along the longer path +178 -> +90 -> 0 -> -90 ->170.
     */

    if (in_theta * p_theta < 0) {
      // Prediction and measurement have different signs

      // Add or subtract 2*PI to move the input theta closer to the predicted value
      if (p_theta > 0) {
        in_theta += 2 * M_PI;
      } else {
        in_theta -= 2 * M_PI;
      }
    }
  }

  VectorXd z_in = VectorXd(3);
  z_in << in_rho, in_theta, in_rho_dot;


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

  //new estimate becomes current state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);

  // update process covariance matrix
  P_ = (I - K * H_) * P_;
}
