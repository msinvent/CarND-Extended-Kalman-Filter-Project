/*
 * sensor.h
 * This class used to define a single sensor
 * Created on: Nov 4, 2018
 * Author: Manish Sharma
 */

#ifndef SENSOR_H_
#define SENSOR_H_

#include "Eigen/Dense"
#include <functional>

class sensor{
public :
	MatrixXd H_; // Measurement Matrix
	MatrixXd R_; // Measurement Covariance Matrix, we get it from the sensor datasheet
	bool linear_; // this flag is used to check the requirement to initialize the nonlinear h
	std::function<MatrixXd(MatrixXd)> h;//nonlinear map from state to measurement;
	std::function<MatrixXd(MatrixXd)> hj;//linear map from state to measurement;
};

#endif /* SENSOR_H_ */
