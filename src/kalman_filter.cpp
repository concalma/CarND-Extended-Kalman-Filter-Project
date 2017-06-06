#include "kalman_filter.h"
#include "tools.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

#define EKF_PI   3.14159265358979323846264338327950288

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}


void KalmanFilter::Predict() {
  /**
   TODO
    * predict the state
  */
  // F_ must be updated before calling this
  x_ = F_ * x_;
  P_ = F_ * P_ * F_.transpose() + Q_;

}

void KalmanFilter::Update(const VectorXd &z) {
    /**
TODO:
     * update the state by using Kalman Filter equations
     */
    VectorXd z_pred = H_ * x_;
    VectorXd y = z - z_pred;
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

VectorXd KalmanFilter::EKF_h(VectorXd &x_state) {
    VectorXd ret = VectorXd(3);

    float px = x_state(0);
    float py = x_state(1);
    float vx = x_state(2);
    float vy = x_state(3);

    ret(0) = sqrt( px*px + py*py );
    ret(1) = atan2( py , px );
    ret(2) = (px*vx + py*vy) / ret(0); 


    return ret;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */

    MatrixXd Hj_ = tools.CalculateJacobian(x_);

    VectorXd z_pred = EKF_h(x_); // using actual measurement function h(X)
    VectorXd y = z - z_pred;
    // y = z - h(x)
    // h(x) = atan2(py/px) can only be -pi .. pi. z can only be -pi .. pi, therefore only need to check for over and under
    y(1) = y(1) > (EKF_PI)  ? y(1)-(2*EKF_PI) : y(1); 
    y(1) = y(1) < (-EKF_PI) ? y(1)+(2*EKF_PI) : y(1);

    MatrixXd Hjt = Hj_.transpose();
    MatrixXd S = Hj_ * P_ * Hjt + R_;
    MatrixXd Si = S.inverse();
    MatrixXd PHjt = P_ * Hjt;
    MatrixXd K = PHjt * Si;

    //new estimate
    x_ = x_ + (K * y);
    long x_size = x_.size();
    MatrixXd I = MatrixXd::Identity(x_size, x_size);
    P_ = (I - K * Hj_) * P_;
}
