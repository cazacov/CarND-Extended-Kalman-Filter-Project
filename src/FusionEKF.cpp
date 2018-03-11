#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
    is_initialized_ = false;

    previous_timestamp_ = 0;
    iteration = 0;

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    float noise_ax = 9;
    float noise_ay = 9;

    A_covariance_  = MatrixXd(2,2);
    A_covariance_ << noise_ax, 0,
                    0, noise_ay;

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
            0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
            0, 0.0009, 0,
            0, 0, 0.09;

    H_laser_ << 1, 0, 0, 0,
                0, 1, 0, 0;

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {

        // first measurement

        double x;
        double y;
        double vx = 0;
        double vy = 0;

        previous_timestamp_ = measurement_pack.timestamp_;

        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
            Convert radar from polar to cartesian coordinates and initialize state.
            */
            double rho = measurement_pack.raw_measurements_[0];
            double theta = measurement_pack.raw_measurements_[1];

            x = rho * cos(theta);
            y = rho * sin(theta);
        } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
            Initialize state.
            */
            x = measurement_pack.raw_measurements_[0]; // X
            y = measurement_pack.raw_measurements_[1]; // Y
        }

        VectorXd state_X = VectorXd(4);
        state_X << x, y, vx, vy;

        MatrixXd state_covariance_P = MatrixXd(4,4);
        state_covariance_P <<   1, 0, 0, 0,
                                0, 1, 0, 0,
                                0, 0, 1000, 0,
                                0, 0, 0, 1000;

        ekf_.Init(state_X, state_covariance_P);

        // done initializing, no need to predict or update
        is_initialized_ = true;
        return;
    }

    iteration++;

    cout << iteration << " ----------------" << endl;
    cout << (measurement_pack.sensor_type_ == MeasurementPackage::RADAR ? "radar " : "laser ") << measurement_pack.raw_measurements_<< endl;


    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds

    MatrixXd state_transition_F = MatrixXd::Identity(4,4);
    state_transition_F(0, 2) = dt;
    state_transition_F(1, 3) = dt;
    ekf_.F_ = state_transition_F;

    // Estimate Q
    float dt_2 = dt * dt;
    MatrixXd g = MatrixXd(4,2);
    g << dt_2 / 2.0,  0,
           0,   dt_2 / 2.0,
           dt,  0,
           0,   dt;
    MatrixXd g_trans = g.transpose();

    MatrixXd noise_covariance_Q = g * A_covariance_ * g_trans;
    ekf_.Q_ = noise_covariance_Q;

    ekf_.Predict();
    float x = ekf_.x_(0);
    float y = ekf_.x_(1);
    cout << "After prediction:" << endl;
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;

    /*****************************************************************************
     *  Update
     ****************************************************************************/

    /**
     TODO:
       * Use the sensor type to perform the update step.
       * Update the state and covariance matrices.
     */

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
        // Radar
        ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
        ekf_.R_ = R_radar_;
        ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    } else {
        // Laser
        ekf_.H_ = H_laser_;
        ekf_.R_ = R_laser_;
        ekf_.Update(measurement_pack.raw_measurements_);
    }

    float xx = ekf_.x_(0);
    float yy = ekf_.x_(1);
    if (fabs(x - xx) > 3 || fabs(y - yy) > 3) {
        cout << "Big change" << endl;
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    // print the output
    cout << "After update:" << endl;
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;

}
