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
 * => Constructors are used to initialise private variables
 */
FusionEKF::FusionEKF() 
{
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  H_laser_ << 	1, 0, 0, 0,
  				0, 1, 0, 0;

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;

  /**
  TODO:
    * Finish initializing the FusionEKF.
    * Set the process and measurement noises
  */


}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) 
{
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

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
    {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
    	float rho = measurement_pack.raw_measurements_[0];
    	float phi = measurement_pack.raw_measurements_[1];
    	float rho_dot = measurement_pack.raw_measurements_[2];
    	/**
		Here we perform the calculations to convert polar coordinates to 
		cartesian coordinates using the Pythogoras theorem
    	*/
    	float x = rho * cos(phi);
    	float y = rho * sin(phi);
    	float vx = rho_dot * cos(phi);
    	float vy = rho_dot * sin(phi);

    	ekf_.x_ << x, y, vx, vy;
    }

    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) 
    {
      /**
      Initialize state.
      */ekf_.x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }

    if (fabs(ekf_.x_(0)) < 0.0001f and fabs(ekf_.x_(1)) < 0.0001f)
    {
    	ekf_.x_(0) = 0.0001f;
    	ekf_.x_(1) = 0.0001f;

    }

    // Initialise covariance matrix
    ekf_.P_ = MatrixXd(4, 4);
    ekf_.P_ << 	1, 0, 0, 0,
    			0, 1, 0, 0,
    			0, 0, 1000, 0,
    			0, 0, 0, 1000;

    previous_timestamp_ = measurement_pack.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
   TODO:
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

  // For better simplicity, the prediction step will be divided into two parts,
  // particularly; 
  // 1. Updating the state transition matrix
  // 2. Updating the noise vector

  /** 
	1. Update state transition matrix
  **/
  // First we calculate the change in time 
  // (current time, measurement_pack.timestamp_ - previous time, previous_timestamp_)

  // Change in time (delta time, dt) is calculated as follows:
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0f;
  previous_timestamp_ = measurement_pack.timestamp_;

  // Initialise the state transition matrix, F with a [ 4 X 4 ] matrix as specified in
  // the notes.
  ekf_.F_ = MatrixXd(4, 4);

  // Update the state transition matrix, F with the latest time difference as specified in
  // the notes
  ekf_.F_ << 	1, 0, dt, 0,
  				0, 1, 0, dt,
  				0, 0, 1, 0,
  				0, 0, 0, 1;

  /**
	2. Update the noise vector
  **/

  // Specify the noise values
  float noise_ax = 9.0f;
  float noise_ay = 9.0f;

  // Since the noise vector, Q = GQv(G)^T contains delta time, dt to the powers of 2 to 4,
  // we specify them here first
  float dt_power_of_2 = dt * dt;
  float dt_power_of_3 = dt_power_of_2 * dt;
  float dt_power_of_4 = dt_power_of_3 * dt;

  // The noise vector, Q also contains the division of the dt_power_of_n with numbers ranging
  // from 2 to 4, we specify them here first as well
  float dt_power_of_3_divided_by_2 = dt_power_of_3 / 2;
  float dt_power_of_4_divided_by_4 = dt_power_of_4 / 4;

  // Initialise the noise vector, Q with a [ 4 by 4 ] matrix as specified in
  // the notes
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << 	dt_power_of_4_divided_by_4 * noise_ax,		0,		dt_power_of_3_divided_by_2 * noise_ax,	0,
  				0,		dt_power_of_4_divided_by_4 * noise_ay,		0,	dt_power_of_3_divided_by_2 * noise_ay,
  				dt_power_of_3_divided_by_2 * noise_ax,		0,			dt_power_of_2 * noise_ax,			0,
  				0,		dt_power_of_3_divided_by_2 * noise_ay,		0,				dt_power_of_2 * noise_ay;						


  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) 
  {
    // Radar updates
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
    ekf_.R_ = R_radar_;

    ekf_.UpdateEKF(measurement_pack.raw_measurements_);
  } 

  else 
  {
    // Laser updates
    ekf_.H_ = H_laser_;
    ekf_.R_ = R_laser_;
    ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
