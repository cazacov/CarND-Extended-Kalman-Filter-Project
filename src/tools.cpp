#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;

  // check the validity of the following inputs:
  //  * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if (estimations.size() != ground_truth.size()
      || estimations.size() == 0) {
    cout << "Invalid estimation or ground_truth data" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (unsigned int i = 0; i < estimations.size(); ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse / estimations.size();

  //calculate the squared root
  rmse = rmse.array().sqrt();

  //return the result
  return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd &x_state) {

  double x = x_state(0);
  double y = x_state(1);
  double vx = x_state(2);
  double vy = x_state(3);

  const double epsilon = 1e-4;

  double xy2 = x * x + y * y;
  double xy12 = sqrt(xy2);
  double xy32 = xy2 * xy12;

  if (fabs(xy2) < epsilon) {
    cout << "Cannot calculate Jacobian" << endl;
    return MatrixXd::Zero(4, 3);
  }

  MatrixXd jacobian = MatrixXd(3, 4);
  jacobian << x / xy12, y / xy12, 0, 0,
          -y / xy2, x / xy2, 0, 0,
          y * (vx * y - vy * x) / xy32, x * (vy * x - vx * y) / xy32, x / xy12, y / xy12;

  return jacobian;
}
