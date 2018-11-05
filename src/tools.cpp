#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
		const vector<VectorXd> &ground_truth) {
	/**
  TODO:
	 * Calculate the RMSE here.
	 */
	VectorXd rmse(4);
	rmse << 0,0,0,0;

	// TODO: YOUR CODE HERE

	// check the validity of the following inputs:
	//  * the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	// ... your code here

	//accumulate squared residuals
	for(int i=0; i < estimations.size(); ++i){
		VectorXd err = estimations[i]-ground_truth[i];
		rmse = rmse + (err.array()*err.array()).matrix();
		// std::cout<<err.array()*err.array()<<"\n\n"<<rmse<<"\n";
		// ... your code here

	}

	//calculate the mean
	rmse /= estimations.size();
	// ... your code here

	//calculate the squared root
	// ... your code here
	rmse = rmse.array().sqrt();
	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state) {
	/**
  TODO:
	 * Calculate a Jacobian here.
	 */
}
