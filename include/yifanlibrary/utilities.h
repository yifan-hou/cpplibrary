#pragma once
#ifndef _MATH_ULTILITIES_H_
#define _MATH_ULTILITIES_H_

#ifdef __cplusplus
    #include <cmath>
    #include <iostream>
#else
    #include <math.h>
    #include <stdio.h>
#endif

// Eigen
#include <Eigen/Geometry>
#include <Eigen/Dense>
// #include <Eigen/Core>
// #include <Eigen/SVD>
// #include <unsupported/Eigen/MatrixFunctions>

namespace UT
{
  /////////////////////////////////////////////////////////////////////////
  //                   types and static variables
  /////////////////////////////////////////////////////////////////////////

  typedef Eigen::Vector3d Vector3d;
  typedef Eigen::Matrix3d Matrix3d;
  typedef Eigen::Matrix4d Matrix4d;
  typedef Eigen::MatrixXd MatrixXd;
  typedef Eigen::VectorXd VectorXd;
  typedef Eigen::Matrix<double, 6, 1> Vector6d;
  typedef Eigen::Matrix<double, 6, 6> Matrix6d;
  typedef Eigen::Quaterniond Quaterniond;

  typedef Eigen::Vector3f Vector3f;
  typedef Eigen::Matrix3f Matrix3f;
  typedef Eigen::Matrix4f Matrix4f;
  typedef Eigen::MatrixXf MatrixXf;
  typedef Eigen::VectorXf VectorXf;
  typedef Eigen::Matrix<float, 6, 1> Vector6f;
  typedef Eigen::Matrix<float, 6, 6> Matrix6f;
  typedef Eigen::Quaternionf Quaternionf;

	/////////////////////////////////////////////////////////////////////////
	//                          iostream
	/////////////////////////////////////////////////////////////////////////
  void stream_array_in(std::ostream &st, double *array, int length);
  void stream_array_in(std::ostream &st, float *array, int length);
  void stream_array_in(std::ostream &st, int *array, int length);

	/////////////////////////////////////////////////////////////////////////
	//                          scalar
	/////////////////////////////////////////////////////////////////////////
  void truncate(double *ele, const double _min, const double _max);

	/////////////////////////////////////////////////////////////////////////
	//                          vector&array
	/////////////////////////////////////////////////////////////////////////

  void buf_insert(const double ele, const int size, double * buf);
  void copyArray(const float *src, float *dest, int dim);
  void copyArray(const double *src, double *dest, int dim);
  void setArray(float *array, float value, int dim);
  void truncate(float *array, float min, float max, int dim);
  double vec_max(const double * vec, const int size);
  double vec_min(const double * vec, const int size);
  double vec_mean(const double * vec, const int size);
  double vec_slope(const double * x, const double * y,const int size);
  // numerical differentiation with low pass filter
  // input x, calculate dx/dt
  // s/(as+1),
  double diff_LPF(const double xdold, const double xnew, const double xold, const double T,const double a);
  void truncate6f(Vector6f *v, float min, float max);
  void truncate6d(Vector6d *v, double min, double max);
  void truncate6d(Vector6d *v, const Vector6d &min, const Vector6d &max);
  void stream_array_in6f(std::ostream &st, const Vector6f &array);
  void stream_array_in6d(std::ostream &st, const Vector6d &array);

  /////////////////////////////////////////////////////////////////////////
  //                          Matrices
  /////////////////////////////////////////////////////////////////////////
  MatrixXd pseudoInverse(const MatrixXd &a,
    double epsilon = std::numeric_limits<double>::epsilon());

  /////////////////////////////////////////////////////////////////////////
  //                          Robotics
  /////////////////////////////////////////////////////////////////////////
  /*  Frames/spaces:
          W: world frame
          T: current tool frame
          So: set tool frame with offset
          Tf: transformed generalized space
      Quantities:
          SE3: 4x4 homogeneous coordinates
          se3: 6x1 twist coordinate of SE3
          spt: 6x1 special twist: 3x1 position, 3x1 exponential coordinate for rotation
          td: 6x1 time derivative of twist.
          v: 6x1 velocity, either spatial or body

          wrench: 6x1 wrench. Makes work with body velocity
  */

  Matrix3d wedge(const Vector3d &v);
  Matrix4d wedge6(const Vector6d &t);

  // Axis-angle to matrix
  // input:
  //   theta: scalar
  //   n: 3 x 1
  // output:
  // m: 3x3
  Eigen::Matrix3d aa2mat(const double theta, const Eigen::Vector3d n);
  Matrix3d quat2SO3(const Quaterniond &q);
  Matrix3d quat2SO3(double qw, double qx, double qy, double qz);
  Matrix3d so32SO3(const Vector3d &v);
  Vector3d SO32so3(const Matrix3d &R);
  void so32quat(const Vector3d &so3, double *q);
  void SO32quat(const Matrix3d &SO3, double *q);
  Matrix4d pose2SE3(const double *pose);
  Matrix4d posemm2SE3(const double *pose);
  Matrix4d se32SE3(const Vector6d &twist);
  Matrix4d spt2SE3(const Vector6d &spt);
  Matrix4d SE3Inv(const Matrix4d &SE3);
  Vector6d SE32se3(const Matrix4d &SE3);
  Vector6d SE32spt(const Matrix4d &SE3);
  Matrix6d SE32Adj(const Matrix4d &SE3);
  void SE32Pose(const Matrix4d &SE3, double *pose);
  void SE32Posemm(const Matrix4d &SE3, double *pose);
  Eigen::Quaternionf QuatMTimes(const Eigen::Quaternionf &q1,
    const Eigen::Quaternionf &q2);
  float angBTquat(Eigen::Quaternionf &q1, Eigen::Quaternionf &q2);
  Eigen::Matrix3f quat2m(const Eigen::Quaternionf &q);
  // Return the 6x6 jacobian matrix mapping from spt time derivative
  //  to body velocity.
  // Jac * spt time derivative = body velocity
  Matrix6d JacobianSpt2BodyV(const Matrix3d &R);
  void double2float(const double *array_in, float *array_out, int n);
  void float2double(const float *array_in, double *array_out, int n);

  /////////////////////////////////////////////////////////////////////////
  //                      Motion Planning
  /////////////////////////////////////////////////////////////////////////
  // remember to allocate 2rd arrays!
  // remember to delete pose_traj!
  void MotionPlanningLinear(const double *pose0, const double *pose_set,
      const int Nsteps, double **pose_traj);

  /**
   * 1D trapezodial interpolation from x0 = 0 to x_f, limited by
   * maximum acceleration @p a_max, maximum velocity @p v_max. The return
   * trajectory is sampled in @p Nsteps time steps. The total duration of the
   * trajectory is 2*t1 + t2.
   *
   * If @p Nsteps=0 (default), the function only computes time @p t1, @p t2,
   * do not generate the trajectory.
   *
   * @param[in]  x_f     The final position
   * @param[in]  a_max   Maximum acceleration
   * @param[in]  v_max   Maximum velocity
   * @param      t1      Time duration of the acceleration phase
   * @param      t2      Time duration of the constant speed phase
   * @param[in]  Nsteps  The number of sampling points
   * @param      x_traj  The interpolated trajectory, 1 x Nsteps array
   */
  void TrapezodialTrajectory(double x_f, double a_max, double v_max, double *t1,
    double *t2, int Nsteps = 0, double * x_traj = 0 );
  /**
   * Pose to pose trapezodial interpolation. Poses are represented as double
   * arrays of length 7 (x, y, z, qw, qx, qy, qz). Unit of translation is mm.
   * Number of frames is determined by the speed/acceleration limit and the
   * control rate @p rate.
   * @param[in]  pose0        The current pose
   * @param[in]  pose_set     The goal pose
   * @param[in]  a_max_trans  Maximum translational acceleration (mm/s^2)
   * @param[in]  v_max_trans  Maximum translational velocity (mm/s)
   * @param[in]  a_max_rot    Maximum rotational acceleration (rad/s^2)
   * @param[in]  v_max_rot    Maximum rotational velocity (rad/s)
   * @param[in]  rate         Control rate (Hz).
   * @param      pose_traj    Output trajectory, a 7 x Nsteps Eigen Matrix.
   */
  void MotionPlanningTrapezodial(const double *pose0, const double *pose_set,
    double a_max_trans, double v_max_trans, double a_max_rot, double v_max_rot,
    double rate, MatrixXd *pose_traj);

  /////////////////////////////////////////////////////////////////////////
  //                      Probability
  /////////////////////////////////////////////////////////////////////////
  double Gaussian(double x, double var);
}

#endif // _MATH_ULTILITIES_H_
