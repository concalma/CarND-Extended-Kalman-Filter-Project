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

    // initializing matrices
    R_laser_ = MatrixXd(2, 2);
    R_radar_ = MatrixXd(3, 3);
    H_laser_ = MatrixXd(2, 4);
    Hj_ = MatrixXd(3, 4);

    //measurement covariance matrix - laser
    R_laser_ << 0.0225, 0,
             0, 0.0225;

    //measurement covariance matrix - radar
    R_radar_ << 0.09, 0, 0,
             0, 0.0009, 0,
             0, 0, 0.09;

    /**
     * Finish initializing the FusionEKF.
     * Set the process and measurement noises
     */
    H_laser_ << 1,0,0,0,
             0,1,0,0;

    ekf_.Q_ = MatrixXd(4, 4);

    ekf_.F_ = MatrixXd(4,4);
    ekf_.F_ << 1,0,0,0,
        0,1,0,0,
        0,0,1,0,
        0,0,0,1;  // we'll update dts on 0,2  1,3 while we predict

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {

    // cout << measurement_pack.raw_measurements_ << endl;

    /*****************************************************************************
     *  Initialization
     ****************************************************************************/
    if (!is_initialized_) {
        /**
TODO:
         * Initialize the state ekf_.x_ with the first measurement.
         * Create the covariance matrix.
         * Remember: you'll need to convert radar from polar to cartesian coordinates.
         */


        // first measurement
        cout << "EKF: " << endl;
        ekf_.x_ = VectorXd(4);
        ekf_.x_ << 1, 1, 1, 1;

        // state covariance matrix P. We use placeholder values '1' for speed, '1000' for uncertainty in velocity
        ekf_.P_ = MatrixXd(4, 4);
        ekf_.P_ << 1, 0, 0, 0,
            0, 1, 0, 0,
            0, 0, 1000, 0,
            0, 0, 0, 1000;

        
        // Updating first state directly with first measurement
        if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
            /**
              Convert radar from polar to cartesian coordinates and initialize state.
              */
            float ro =    measurement_pack.raw_measurements_[0];
            float theta = measurement_pack.raw_measurements_[1];
            float rodot = measurement_pack.raw_measurements_[2];

            ekf_.x_ << ro*cos(theta), ro*sin(theta), 0, 0;  // we drop rodot info according to project tips

        }
        else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
            /**
              Initialize state.
              */
            //set the state with the initial location and zero velocity
            ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;

        }

        // done initializing, no need to predict or update
        is_initialized_ = true;
        previous_timestamp_ = measurement_pack.timestamp_;
        return;
    }

    /*****************************************************************************
     *  Prediction
     ****************************************************************************/

    /**
     * Update the state transition matrix F according to the new elapsed time.
     - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
     */

    // if no previous timestamp take measurement as is. Consider dt=0. This triggers F = Id and Q =0 first time around
    float dt = 0;
    if( previous_timestamp_ != 0 ) {
        //compute the time elapsed between the current and previous measurements
        dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
    }

    previous_timestamp_ = measurement_pack.timestamp_;

    float dt_2 = dt * dt;
    float dt_3 = dt_2 * dt;
    float dt_4 = dt_3 * dt;

    //Modify the F matrix so that the time is integrated
    // this assumes proper initialization to 1s in constructor
    ekf_.F_(0, 2) = dt;
    ekf_.F_(1, 3) = dt;

    //set the process covariance matrix Q
    float noise_ax =  9;
    float noise_ay =  9;

    ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
        0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
        dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
        0, dt_3/2*noise_ay, 0, dt_2*noise_ay;


    ekf_.Predict();

    /*****************************************************************************
     *  Update
     ****************************************************************************/
    // comment out desired mode
    //Mode mode = MODE_RADAR;
    //Mode mode = MODE_LASER;
    Mode mode = MODE_FUSION;

    /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
     */
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR && mode!=MODE_LASER ) {
        
            ekf_.R_ = R_radar_;
            // Jacobian calculation and h(x) is done inside KarmanFilter class
            ekf_.UpdateEKF( measurement_pack.raw_measurements_);

    } else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER && mode!=MODE_RADAR ) {
            ekf_.H_ = H_laser_;
            ekf_.R_ = R_laser_;
            ekf_.Update(measurement_pack.raw_measurements_);
    }


    // print the output
    cout << "x_ = " << ekf_.x_ << endl;
    cout << "P_ = " << ekf_.P_ << endl;
}
