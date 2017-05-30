#include "kalman_filter.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) 
{
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() 
{
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z) 
{
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
  VectorXd y = z - H_ * x_;
  
  MatrixXd Ht = H_.transpose();
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  x_ = x_ + (K * y);
  int x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) 
{
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float rho = sqrt(x_(0) * x_(0) + x_(1) * x_(1));
  float theta = atan2 (x_(1), x_(0));
  float rho_dot = (x_(0) * x_(2) + x_(1) * x_(3)) / rho;

  VectorXd h = VectorXd(3);
  h << rho, theta, rho_dot;
  
  MatrixXd Ht = H_.transpose();
  VectorXd y = z - h;
  MatrixXd S = H_ * P_ * Ht + R_;
  MatrixXd K = P_ * Ht * S.inverse();
  int size = x_.size();
  MatrixXd I = MatrixXd::Identity(size, size);
  
  // At this point, we normalise the angles such that the resulting
  // angles lie between -pi and pi since the Kalman Filter is expecting small
  // values between -pi and pi
  y[1] -= (2 * M_PI) * floor((y[1] + M_PI) / (2 * M_PI));

  x_ = x_ + K * y;
  P_ = (I - K * H_) * P_;
}