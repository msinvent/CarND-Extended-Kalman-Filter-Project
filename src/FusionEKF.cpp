#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
#include <cmath>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static const double PI = 3.141592653589;

// nonlinear map from state to measurement
VectorXd radar_h(const VectorXd& x_state){
	VectorXd hx(3);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = sqrt(px*px+py*py);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return VectorXd::Zero(3);
	}

	//check division by zero issue
	if(fabs(py) < 0.001){
		hx << c1,atan(py/px),PI/2,(px*vx+py*vy)/c1;
		return hx;
	}

	//compute the Jacobian matrix
	hx << c1,atan(py/px),(px*vx+py*vy)/c1;
	return hx;
}

// linear map from state to measurement
MatrixXd radar_hj(const VectorXd& x_state){
	MatrixXd Hj = MatrixXd::Zero(3,4);
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	//pre-compute a set of terms to avoid repeated calculation
	float c1 = px*px+py*py;
	float c2 = sqrt(c1);
	float c3 = (c1*c2);

	//check division by zero
	if(fabs(c1) < 0.0001){
		cout << "CalculateJacobian () - Error - Division by Zero" << endl;
		return Hj;
	}

	//compute the Jacobian matrix
	Hj << (px/c2), (py/c2), 0, 0,
			-(py/c1), (px/c1), 0, 0,
			py*(vx*py - vy*px)/c3, px*(px*vy - py*vx)/c3, px/c2, py/c2;

	return Hj;
}

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {

  is_initialized_ = false;
  previous_timestamp_ = 0;

	sensor laser;
	sensor radar;

	// defining laser sensor
	laser.linear_ = true;
	laser.R_ = MatrixXd::Zero(2,2);
	laser.R_ << 0.0225, 0,
	        0, 0.0225;
	laser.H_ = MatrixXd::Zero(2,4);
	laser.H_ << 1, 0, 0, 0,
			  0, 1, 0, 0;
	radar.h = nullptr;
	radar.hj = nullptr;
	sensors["laser"] = laser;

	// defining radar sensor
	radar.linear_ = false;
	laser.H_ = MatrixXd::Zero(3,4);
	radar.R_ = MatrixXd::Zero(3,3);
	radar.R_ << 0.09, 0, 0,
	        0, 0.0009, 0,
	        0, 0, 0.09;
	radar.h = radar_h;
	radar.hj = radar_hj;
	sensors["radar"] = radar;

	P_ = MatrixXd::Zero(4, 4);
	P_ << 1, 0, 0, 0,
				  0, 1, 0, 0,
				  0, 0, 1000, 0,
				  0, 0, 0, 1000;

  x_ = VectorXd::Zero(4);
  x_ << 1, 1, 1, 1;

	//the initial transition matrix F_
	F_ = MatrixXd::Zero(4, 4);
	F_ << 1, 0, 1, 0,
			  0, 1, 0, 1,
			  0, 0, 1, 0,
			  0, 0, 0, 1;

	//set the acceleration noise components
	noise_ax = 5;
	noise_ay = 5;
}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::Predict() {
	x_ = F_ * x_;
	P_ = F_ * P_ * F_.transpose() + Q_;
}

void FusionEKF::Update(const sensor& Sensor,const VectorXd& measurement) {
	// innovation calculation
	VectorXd y;
	MatrixXd H;

	if(Sensor.linear_){
		y = measurement-Sensor.H_*x_;
		H = Sensor.H_;
	}else{
		y = measurement-Sensor.h(x_);
		H = Sensor.hj(x_);
	}

	MatrixXd S = H * P_ * H.transpose() + Sensor.R_;
	MatrixXd K = P_ * H.transpose()* S.inverse();
	x_ = x_ + K*y;
	P_ = (MatrixXd::Identity(K.rows(),H.cols()) - K*H)*P_;
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {
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

    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
  		x_ << measurement_pack.raw_measurements_[0]*cos(measurement_pack.raw_measurements_[1]*(PI/180.0)),
  				measurement_pack.raw_measurements_[0]*sin(measurement_pack.raw_measurements_[1]*(PI/180.0)),
					0, 0;
    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
//    	std::cout<<measurement_pack.raw_measurements_<<"\n";
  		x_ << measurement_pack.raw_measurements_[0], measurement_pack.raw_measurements_[1], 0, 0;
    }
//
//    // done initializing, no need to predict or update
		previous_timestamp_ = measurement_pack.timestamp_;
    is_initialized_ = true;
    return;
  }

	//compute the time elapsed between the current and previous measurements
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  //Modify the F matrix so that the time is integrated
  F_(0, 2) = dt;
  F_(1, 3) = dt;

  //set the process covariance matrix Q
  Q_ = MatrixXd(4, 4);
  Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
  		0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
			dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
			0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  //predict
  Predict();
//  P_ = F_ * P_ * F_.transpose() + Q_;
//  std::cout<<"dt : \n"<<dt<<"\n\n";
//  std::cout<<"Q_ : \n"<<Q_<<"\n\n";
//  std::cout<<"x_ : \n"<<x_<<"\n\n";
//  std::cout<<"P_ : \n"<<P_<<"\n\n";
//  std::cout<<"F_ : \n"<<F_<<"\n\n";


//  std::cout<<"Manish Sharma";
  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
   TODO:
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
//  	std::cout<<"radar : \n"<<measurement_pack.raw_measurements_<<"\n\n";
  	Update(sensors["radar"],measurement_pack.raw_measurements_);
    // Radar updates
  } else if(measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
//  	std::cout<<"laser : \n"<<measurement_pack.raw_measurements_<<"\n\n";
  	Update(sensors["laser"],measurement_pack.raw_measurements_);
    // Laser updates
  }
//
//  // print the output
//  cout << "x_ = " << ekf_.x_ << endl;
//  cout << "P_ = " << ekf_.P_ << endl;
}
