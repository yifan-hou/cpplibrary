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

#include <memory> // smart pointers
#include <vector>

// Eigen
#include <Eigen/Geometry>
#include <Eigen/Dense>

namespace RUT
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
  typedef Eigen::Quaterniond Quaterniond;

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
  /**
   * Find in the vector the element with maximum abs value.
   *
   * @param[in]  vec   The vector.
   * @param[in]  size  Size of the vector.
   * @param      id    Id of the maximum element
   *
   * @return     The maximum element
   */
  double vec_max_abs(const double * vec, const int size, int *id);
  // numerical differentiation with low pass filter
  // input x, calculate dx/dt
  // s/(as+1),
  double diff_LPF(const double xdold, const double xnew, const double xold, const double T,const double a);
  void truncate6f(Vector6f *v, float min, float max);
  void truncate6d(Vector6d *v, double min, double max);
  void truncate6d(Vector6d *v, const Vector6d &min, const Vector6d &max);
  void stream_array_in6f(std::ostream &st, const Vector6f &array);
  void stream_array_in6d(std::ostream &st, const Vector6d &array);
  void double2float(const double *array_in, float *array_out, int n);
  void float2double(const float *array_in, double *array_out, int n);

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

  //
  // Transformations
  //

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
  Eigen::Matrix3f quat2m(const Eigen::Quaternionf &q);
  Eigen::Matrix3d rotX(double angle_rad);
  Eigen::Matrix3d rotY(double angle_rad);
  Eigen::Matrix3d rotZ(double angle_rad);

  //
  // Quaternions
  //
  Eigen::Quaternionf QuatMTimes(const Eigen::Quaternionf &q1,
    const Eigen::Quaternionf &q2);
  Eigen::Quaterniond QuatMTimes(const Eigen::Quaterniond &q1,
      const Eigen::Quaterniond &q2);
  float angBTquat(const Eigen::Quaternionf &q1, const Eigen::Quaternionf &q2);
  double angBTquat(const Eigen::Quaterniond &q1, const Eigen::Quaterniond &q2);

  //
  // SE(3)
  //
  class CartesianPose
  {
  public:
    CartesianPose(){};
    /**
     * Construct the pose from a 1x7 vector.
     *
     * @param[in]  pose  Pose vector. [x y z qw qx qy qz]
     */
    CartesianPose(std::vector<double> pose);
    /**
     * Constructs the pose from a 4x4 homogeneous matrix.
     *
     * @param[in]  T     The 4x4 homogeneous matrix.
     */
    CartesianPose(Eigen::Matrix4d T);
    CartesianPose(const Eigen::Quaterniond &q, const Eigen::Vector3d &p);
    CartesianPose(const Eigen::Matrix3d &R, const Eigen::Vector3d &p);
    ~CartesianPose(){}

    static CartesianPose Identity();

    // types
    using Ptr = std::shared_ptr<CartesianPose>;
    using ConstPtr = std::shared_ptr<const CartesianPose>;

    // setters/getters
    void setQuaternion(const Eigen::Quaterniond &q);
    void setQuaternion(const std::vector<double> &q);
    void setXYZ(const Eigen::Vector3d &p);
    void setXYZ(const std::vector<double> &p);
    Eigen::Matrix3d getRotationMatrix() const;
    Eigen::Quaterniond getQuaternion() const;
    Eigen::Vector3d getXYZ() const;
    Eigen::Vector3d getXAxis() const;
    Eigen::Vector3d getYAxis() const;
    Eigen::Vector3d getZAxis() const;
    Eigen::Matrix4d getTransformMatrix() const;
    Eigen::Isometry3d getIsometry3d() const;

    // operators
    CartesianPose operator*(const CartesianPose &pose) const;
    CartesianPose inv() const;

    // Transformations
    Eigen::Vector3d transformVec(const Eigen::Vector3d &v) const;
    Eigen::Vector3d transformPoint(const Eigen::Vector3d &p) const;
    Eigen::Quaterniond transformQuat(const Eigen::Quaterniond &q) const;

    // metric

    /**
     * Distance between this pose and another pose
     *
     * @param[in]  pose   Another pose
     * @param[in]  ratio  scaling of rotation to length
     *
     * @return     angle*ratio + length
     */
    double distBTPose(const CartesianPose & pose, double ratio = 1.0) const;

    // MISC
    void print() const;
  private:
    Eigen::Quaterniond q_;
    Eigen::Vector3d p_;
    Eigen::Matrix3d R_;
  };

  //
  // Vector operations
  //
  /**
   * find the (signed) angle from vector x to vector b.
   *
   * @param[in]  x            initial vector.
   * @param[in]  b            final vector.
   * @param[in]  z            if specified, it is the rotation axis. It is used
   *                          to specify positive direction of rotation.
   * @param[in]  nonnegative  when z is present, nonnegative will decide the
   *                          range of return angle between [-pi, pi] or [0, 2pi]
   *
   * @return     The angle. If @p z is not given, range is [0, pi]. If @p z is
   *             given, range is [-pi, pi] (@p nonnegative = false) or
   *             [0, 2pi] (@p nonnegative = true).
   */
  double angBTVec(Eigen::Vector3d x, Eigen::Vector3d b,
      Eigen::Vector3d z = Eigen::Vector3d::Zero(), bool nonnegative = false);
  Eigen::MatrixXd transformByRAndP(const Eigen::MatrixXd &points_rowwise,
     const Eigen::Matrix3d &R, const Eigen::Vector3d &p);

  //
  // Others
  //

  // Return the 6x6 jacobian matrix mapping from spt time derivative
  //  to body velocity.
  // Jac * spt time derivative = body velocity
  Matrix6d JacobianSpt2BodyV(const Matrix3d &R);
  /////////////////////////////////////////////////////////////////////////
  //                      Motion Planning
  /////////////////////////////////////////////////////////////////////////

  /**
   * Linear interpolation between @p pose0 and @p pose_set.
   *
   * @param[in]  pose0      The initial pose. [x y z qw qx qy qz]
   * @param[in]  pose_set   The goal pose
   * @param[in]  Nsteps     The number of time steps
   * @param      pose_traj  The interpolated pose traj. 7 x Nsteps matrix
   */
  void MotionPlanningLinear(const double *pose0, const double *pose_set,
      const int Nsteps, MatrixXd *pose_traj);

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
