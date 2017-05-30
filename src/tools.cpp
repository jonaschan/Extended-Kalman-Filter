#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) 
{
  /**
  TODO:
    * Calculate the RMSE here.
  */

	VectorXd rmse(4);

	// Initialise rmse to zeroes matrix
	rmse << 0, 0, 0, 0;

	// Check the validity of the following inputs:
	// 1. Estimation vector size should not be zero
	// 2. Estimation vector size should be equal to the size of the ground truth vector

	if (estimations.size() == 0 || estimations.size() != ground_truth.size())
	{
		cout << "Input is either empty or not equal to the ground truth size" << endl;
		return rmse;
	}

	// Accumulate squared residuals
	for (unsigned int i = 0; i < estimations.size(); i++)
	{
		VectorXd residual = estimations[i] - ground_truth[i];

		// Coefficient wise multiplication
		residual = residual.array() * residual.array();
		rmse += residual;
	}

	//Calculate the mean
	rmse = rmse / estimations.size();

	//Calculate the square root of the mean
	rmse = rmse.array().sqrt();
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) 
{
  /**
  TODO:
    * Calculate a Jacobian here.
  */
	// Since the state consists of the combination of the position and
	// velocity states, we first break them down into individual variables
	// for easier calculation of the values individually

	float position_x = x_state(0);
	float position_y = x_state(1);
	float velocity_x = x_state(2);
	float velocity_y = x_state(3);

	// Create a Jacobian matrix variable, Hj of size [ 3 X 4 ]
	MatrixXd Hj(3, 4);

	if (fabs(position_x) < 0.0001f and fabs(position_y) < 0.0001f)
	{
		position_x = 0.0001f;
		position_y = 0.0001f;
	}

	// Since the Jacobian Matrix consists of certain constants such as 
	// position to the power of two and the square root of position, the values
	// are calculated here first
	float position_x_power_of_2 = position_x * position_x;
	float position_y_power_of_2 = position_y * position_y;

	float sum_of_position_power = position_x_power_of_2 + position_y_power_of_2;

	// Since we don't want the calculation to break due to division by 0,
	// we check the for this case and return the smallest value possible, e.g. 0.0000001
	if (fabs(sum_of_position_power) < 0.0000001f)
	{
		sum_of_position_power = 0.0000001f;
	}

	float square_root_of_position_powers = sqrt(sum_of_position_power);
	float sum_of_position_power_to_power_of_3_by_2 = square_root_of_position_powers * sum_of_position_power;

	// Here we calculate the Jacobian Matrix as specified in the notes
	Hj << 	position_x / square_root_of_position_powers,	position_y / square_root_of_position_powers,	0, 	0,
			-(position_y / sum_of_position_power), 			position_x / sum_of_position_power, 			0, 	0,
			position_y * (velocity_x * position_y - velocity_y * position_x) / sum_of_position_power_to_power_of_3_by_2,	
			position_x * (-velocity_x * position_y + velocity_y * position_x) / sum_of_position_power_to_power_of_3_by_2,
			position_x / square_root_of_position_powers, position_y / square_root_of_position_powers;

	return Hj;  
}
