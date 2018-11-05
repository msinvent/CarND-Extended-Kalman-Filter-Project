#ifndef FusionEKF_H_
#define FusionEKF_H_

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>
#include <unordered_map>
#include "kalman_filter.h"
#include "tools.h"
#include "sensor.h"

class FusionEKF {
public:
	///* state vector
	VectorXd x_;

	///* state covariance matrix
	MatrixXd P_;

	///* state transistion matrix
	MatrixXd F_;

	///* process covariance matrix
	MatrixXd Q_;

  /**
  * Constructor.
  */
  FusionEKF();

  /**
  * Destructor.
  */
  virtual ~FusionEKF();

  /**
  * Run the whole flow of the Kalman Filter from here.
  */
  void ProcessMeasurement(const MeasurementPackage &measurement_pack);


  void Predict();
  void Update(const sensor&,const VectorXd&);
  /**
  * Kalman Filter update and prediction math lives in here.
  */
  KalmanFilter ekf_;
  std::unordered_map<std::string,sensor> sensors;

private:
  // check whether the tracking toolbox was initialized or not (first measurement)
  bool is_initialized_;

  // previous timestamp
  long long previous_timestamp_;

	//acceleration noise components
	float noise_ax;
	float noise_ay;

  // tool object used to compute Jacobian and RMSE
  Tools tools;
  Eigen::MatrixXd R_laser_;
  Eigen::MatrixXd R_radar_;
  Eigen::MatrixXd H_laser_;
  Eigen::MatrixXd Hj_;
};

#endif /* FusionEKF_H_ */
