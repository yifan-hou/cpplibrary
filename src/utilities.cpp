#include "RobotUtilities/utilities.h"


// Eigen
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/SVD>
#include <unsupported/Eigen/MatrixFunctions>

namespace RUT
{
  /////////////////////////////////////////////////////////////////////////
  //                   types and variables
  /////////////////////////////////////////////////////////////////////////
  const static double kEpsilon = 1e-7;
  const static float PIf = 3.1416f;
  const static double PI = 3.1415926;

  /////////////////////////////////////////////////////////////////////////
  //                          iostream
  /////////////////////////////////////////////////////////////////////////
  void stream_array_in(std::ostream &st, double *array, int length)
  {
    for (int i = 0; i<length; i++)
    {
     st << array[i];
     st << "\t";
   }
 }

 void stream_array_in(std::ostream &st, float *array, int length)
 {
  for (int i = 0; i<length; i++)
  {
   st << array[i];
   st << "\t";
 }
}


void stream_array_in(std::ostream &st, int *array, int length)
{
  for (int i = 0; i<length; i++)
  {
   st << array[i];
   st << "\t";
 }
}

    /////////////////////////////////////////////////////////////////////////
    //                          scalar
    /////////////////////////////////////////////////////////////////////////

void truncate(double *ele, const double _min, const double _max)
{
  if ( (*ele) > _max)
    (*ele) = _max;
  else if ( (*ele) < _min)
    (*ele) = _min;
}

    /////////////////////////////////////////////////////////////////////////
    //                          vector&array
    /////////////////////////////////////////////////////////////////////////

void buf_insert(const double ele, const int size, double * buf)
{
  for (int i = 1; i < size; ++i)
  {
    buf[size - i] = buf[size - 1 - i];
  }
  buf[0] = ele;
}

void copyArray(const float *src, float *dest, int dim)
{
  for(int i = 0; i<dim; i++)
  {
    dest[i] = src[i];
  }
}

void copyArray(const double *src, double *dest, int dim)
{
  for(int i = 0; i<dim; i++)
  {
    dest[i] = src[i];
  }
}

void setArray(float *array, float value, int dim)
{
 for(int i=0; i<dim; i++)
 {
   array[i] = value;
 }
}

void truncate(float *array, float min, float max, int dim)
{
 for(int i=0; i<dim; i++)
 {
  array[i] = (array[i] > max)? max:array[i];
  array[i] = (array[i] < min)? min:array[i];
  }
}

double vec_max(const double * vec, const int size)
{
  double m = vec[0];
  double t1;
  for (int i = 0; i < size; ++i)
  {
   t1 = vec[i];
   if (t1 > m) m = t1;
 }
 return m;
}

double vec_min(const double * vec, const int size)
{
  double m = vec[0];
  double t1;
  for (int i = 0; i < size; ++i)
  {
   t1 = vec[i];
   if (t1 < m) m = t1;
 }
 return m;
}

double vec_max_abs(const double * vec, const int size, int *id)
  {
    double m = vec[0];
    double t1;
    for (int i = 0; i < size; ++i)
    {
      t1 = fabs(vec[i]);
      if (t1 > m)
      {
        *id = i;
        m = t1;
      }
    }
    return m;
  }

double vec_mean(const double * vec, const int size)
{
  double sum = 0;
  for (int i = 0; i < size; ++i)
  {
   sum += vec[i];
 }
 return sum/double(size);
}


double vec_slope(const double * x, const double * y,const int size)
{
 double avgX = vec_mean(x,size);
 double avgY = vec_mean(y,size);

 double numerator = 0.0;
 double denominator = 0.0;

 double xd = 0;
 for(int i=0; i<size; ++i)
 {
  xd = x[i] - avgX;
  numerator += (xd) * (y[i] - avgY);
  denominator += xd * xd;
}

return numerator / denominator;
}

    // numerical differentiation with low pass filter
    // input x, calculate dx/dt
    // s/(as+1),
double diff_LPF(const double xdold, const double xnew, const double xold, const double T,const double a)
{
  double As = exp(-T/a);
  return As*xdold + (1 - As)*((xnew - xold)/T);
}

void truncate6f(Vector6f *v, float min, float max)
{
  for(int i=0; i<6; i++)
  {
    (*v)[i] = ((*v)[i] > max)? max:(*v)[i];
    (*v)[i] = ((*v)[i] < min)? min:(*v)[i];
  }
}

void truncate6d(Vector6d *v, double min, double max)
{
  for(int i=0; i<6; i++)
  {
    (*v)[i] = ((*v)[i] > max)? max:(*v)[i];
    (*v)[i] = ((*v)[i] < min)? min:(*v)[i];
  }
}

void truncate6d(Vector6d *v, const Vector6d &min, const Vector6d &max)
{
  for(int i=0; i<6; i++)
  {
    (*v)[i] = ((*v)[i] > max[i])? max[i]:(*v)[i];
    (*v)[i] = ((*v)[i] < min[i])? min[i]:(*v)[i];
  }
}

void stream_array_in6f(std::ostream &st, const Vector6f &array)
{
  for (int i = 0; i<6; i++)
  {
    st << array(i);
    st << "\t";
  }
}
void stream_array_in6d(std::ostream &st, const Vector6d &array)
{
  for (int i = 0; i<6; i++)
  {
    st << array(i);
    st << "\t";
  }
}

/////////////////////////////////////////////////////////////////////////
//                          Matrices
/////////////////////////////////////////////////////////////////////////
MatrixXd pseudoInverse(const MatrixXd &a,
    double epsilon) {
  if (a.norm() < epsilon) {
    return  MatrixXd::Zero(a.cols(), a.rows());
  } else {
    Eigen::JacobiSVD< MatrixXd > svd(a ,Eigen::ComputeThinU | Eigen::ComputeThinV);
    double tolerance = epsilon * std::max(a.cols(), a.rows()) *svd.singularValues().array().abs()(0);
    return svd.matrixV() * (svd.singularValues().array().abs() > tolerance).select(svd.singularValues().array().inverse(), 0).matrix().asDiagonal() * svd.matrixU().adjoint();
  }
}

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

Matrix3d wedge(const Vector3d &v) {
  Matrix3d v_wedge;
  v_wedge << 0, -v(2), v(1),
  v(2), 0, -v(0),
  -v(1), v(0), 0;
  return v_wedge;
}

Matrix4d wedge6(const Vector6d &t) {
  Matrix4d t_wedge;
  t_wedge <<   0,   -t(5),   t(4),  t(0),
  t(5),     0,   -t(3),  t(1),
  -t(4),   t(3),     0,   t(2),
  0,     0,     0,     0;
  return t_wedge;
}

// Axis-angle to matrix
// input:
//   theta: scalar
//   n: 3 x 1
// output:
// m: 3x3
Eigen::Matrix3d aa2mat(const double theta, const Eigen::Vector3d n) {
  Eigen::Vector3d n_unit= n.normalized();
  Eigen::Matrix3d N = wedge(n_unit);
  // m = eye(3) + sin(theta)*N + (1-cos(theta))*N*N;
  Eigen::Matrix3d m = Eigen::Matrix3d::Identity() + std::sin(theta)*N +
      (1-std::cos(theta))*N*N;
  return m;
}

Matrix3d quat2SO3(const Quaterniond &q) {
  return q.normalized().toRotationMatrix();
}

Matrix3d quat2SO3(double qw, double qx, double qy, double qz) {
  Quaterniond q(qw, qx, qy, qz);
  return q.normalized().toRotationMatrix();
}

Matrix3d so32SO3(const Vector3d &v) {
  double theta = v.norm();
  if (theta > kEpsilon) {
    Vector3d vn = v/theta;
    Matrix3d v_wedge = wedge(v);
    Matrix3d SO3;
    SO3 = Matrix3d::Identity() + v_wedge*sin(theta) +
      v_wedge*v_wedge*(1.0 - cos(theta));
    return SO3;
  } else {
    return Matrix3d::Identity();
  }
}

Vector3d SO32so3(const Matrix3d &R) {
  Vector3d so3;
  double temp_arg_to_cos = (R.trace() - 1.0)/2.0;
  truncate(&temp_arg_to_cos, -1.0, 1.0);
  double theta = acos(temp_arg_to_cos);
  if(fabs(theta) < kEpsilon) {
    so3(0) = 1.0;
    so3(1) = 0.0;
    so3(2) = 0.0;
  } else {
    so3(0) = R(2,1)-R(1,2);
    so3(1) = R(0,2)-R(2,0);
    so3(2) = R(1,0)-R(0,1);
    so3 /= 2.0*sin(theta);
  }
  so3 *= theta;
  return so3;
}

void so32quat(const Vector3d &so3, double *q) {
  double theta = so3.norm();
  if (theta < kEpsilon) {
    q[0] = 1;
    q[1] = 0;
    q[2] = 0;
    q[3] = 0;
  } else {
            // q = [cos(theta/2); sin(theta/2)*so3/theta];
    double sin_theta = sin(theta/2.0)/theta;
    q[0] = cos(theta/2.0);
    q[1] = so3(0)*sin_theta;
    q[2] = so3(1)*sin_theta;
    q[3] = so3(2)*sin_theta;
  }
}
void SO32quat(const Matrix3d &SO3, double *q) {
  Quaterniond q_eigen(SO3);
  q_eigen.normalize();
  q[0] = q_eigen.w();
  q[1] = q_eigen.x();
  q[2] = q_eigen.y();
  q[3] = q_eigen.z();
}

Matrix4d pose2SE3(const double *pose) {
  Matrix4d SE3 = Matrix4d::Identity();
  SE3(0, 3) = pose[0];
  SE3(1, 3) = pose[1];
  SE3(2, 3) = pose[2];
  SE3.block<3,3>(0,0) = quat2SO3(pose[3], pose[4], pose[5], pose[6]);
  return SE3;
}

Matrix4d posemm2SE3(const double *pose) {
  Matrix4d SE3 = Matrix4d::Identity();
  SE3(0, 3) = pose[0]/1000.0;
  SE3(1, 3) = pose[1]/1000.0;
  SE3(2, 3) = pose[2]/1000.0;
  SE3.block<3,3>(0,0) = quat2SO3(pose[3], pose[4], pose[5], pose[6]);
  return SE3;
}

Matrix4d se32SE3(const Vector6d &twist) {
  Matrix4d SE3 = Matrix4d::Identity();
  double theta = twist.tail(3).norm();
  if ( theta < kEpsilon ) {
    // no rotation
    SE3(0, 3) = twist(0);
    SE3(1, 3) = twist(1);
    SE3(2, 3) = twist(2);
  } else {
    Vector3d v = twist.head(3);
    Vector3d w = twist.tail(3);
    Matrix3d R = so32SO3(w);
    v /= theta;
    w /= theta;
    SE3.block<3,3>(0, 0) = R;
    SE3.block<3,1>(0, 3) = (Matrix3d::Identity() - R)*(w.cross(v)) +
        w*w.transpose()*v*theta;
  }
  return SE3;
}

Matrix4d spt2SE3(const Vector6d &spt) {
  Matrix4d SE3 = Matrix4d::Identity();
  SE3.block<3, 3>(0, 0) = so32SO3(spt.tail(3));
  SE3.block<3, 1>(0, 3) = spt.head(3);
  return SE3;
}

Matrix4d SE3Inv(const Matrix4d &SE3) {
  Matrix4d SE3_inv = Matrix4d::Identity();
  SE3_inv.block<3,1>(0, 3) =
      -SE3.block<3,3>(0,0).transpose()*SE3.block<3,1>(0,3);
  SE3_inv.block<3,3>(0,0) = SE3.block<3,3>(0,0).transpose();
  return SE3_inv;
}

Vector6d SE32se3(const Matrix4d &SE3) {
  Vector3d p     = SE3.block<3,1>(0, 3);
  Vector3d omega = SO32so3(SE3.block<3,3>(0,0));
  double theta = omega.norm();
  if (theta < kEpsilon) {
    Vector6d se3;
    se3 << p(0), p(1), p(2), 0, 0, 0;
    return se3;
  } else {
    omega /= theta;
    Matrix3d M =
        (Matrix3d::Identity() - wedge(omega*theta).exp())*
        wedge(omega)+omega*omega.transpose()*theta;
    Vector6d se3;
    se3.head(3) = M.fullPivLu().solve(p);
    se3.tail(3) = omega;
    se3 *= theta;
    return se3;
  }
}

Vector6d SE32spt(const Matrix4d &SE3) {
  Vector6d spt;
  spt.head(3) = SE3.block<3, 1>(0, 3);
  spt.tail(3) = SO32so3(SE3.block<3, 3>(0, 0));
  return spt;
}

Matrix6d SE32Adj(const Matrix4d &SE3) {
  Matrix6d Adj = Matrix6d::Zero();
  Adj.topLeftCorner(3, 3)     = SE3.topLeftCorner(3, 3);
  Adj.bottomRightCorner(3, 3) = SE3.topLeftCorner(3, 3);
  Adj.topRightCorner(3, 3)    =
      wedge(SE3.block<3,1>(0, 3)) * SE3.topLeftCorner(3, 3);
  return Adj;
}

void SE32Pose(const Matrix4d &SE3, double *pose) {
  pose[0] = SE3(0, 3);
  pose[1] = SE3(1, 3);
  pose[2] = SE3(2, 3);
  SO32quat(SE3.block<3,3>(0,0), pose + 3);
}

void SE32Posemm(const Matrix4d &SE3, double *pose) {
  pose[0] = SE3(0, 3)*1000.0;
  pose[1] = SE3(1, 3)*1000.0;
  pose[2] = SE3(2, 3)*1000.0;
  SO32quat(SE3.block<3,3>(0,0), pose + 3);
}

Eigen::Matrix3f quat2m(const Eigen::Quaternionf &q) {
  float q11 = q.x()*q.x();
  float q22 = q.y()*q.y();
  float q33 = q.z()*q.z();
  float q01 = q.w()*q.x();
  float q02 = q.w()*q.y();
  float q03 = q.w()*q.z();
  float q12 = q.x()*q.y();
  float q13 = q.x()*q.z();
  float q23 = q.y()*q.z();

  Eigen::Matrix3f m;
  m << 1.0f - 2.0f*q22 - 2.0f*q33, 2.0f*(q12 - q03),      2.0f*(q13 + q02),
      2.0f*(q12 + q03),     1.0f - 2.0f*q11 - 2.0f*q33,  2.0f*(q23 - q01),
      2.0f*(q13 - q02),     2.0f*(q23 + q01),      1.0f - 2.0f*q11 - 2.0f*q22;
  return m;
}

Matrix6d JacobianSpt2BodyV(const Matrix3d &R) {
  Matrix6d Jac;
  Jac = Matrix6d::Identity();
  Jac(3, 3) = R(0,2)*R(0,2) + R(1,2)*R(1,2) + R(2,2)*R(2,2);
  Jac(3, 5) = -R(0,0)*R(0,2) - R(1,0)*R(1,2) - R(2,0)*R(2,2);
  Jac(4, 3) = -R(0,0)*R(0,1) - R(1,0)*R(1,1) - R(2,0)*R(2,1);
  Jac(4, 4) = R(0,0)*R(0,0) + R(1,0)*R(1,0) + R(2,0)*R(2,0);
  Jac(5, 4) = -R(0,1)*R(0,2) - R(1,1)*R(1,2) - R(2,1)*R(2,2);
  Jac(5, 5) = R(0,1)*R(0,1) + R(1,1)*R(1,1) + R(2,1)*R(2,1);

  return Jac;
}

Eigen::Matrix3d rotX(double angle_rad) {
  Eigen::Vector3d x;
  x << 1, 0, 0;
  return aa2mat(angle_rad, x);
}
Eigen::Matrix3d rotY(double angle_rad) {
  Eigen::Vector3d y;
  y << 0, 1, 0;
  return aa2mat(angle_rad, y);
}
Eigen::Matrix3d rotZ(double angle_rad) {
  Eigen::Vector3d z;
  z << 0, 0, 1;
  return aa2mat(angle_rad, z);
}

Eigen::Quaternionf QuatMTimes(const Eigen::Quaternionf &q1,
    const Eigen::Quaternionf &q2)  {
  float s1 = q1.w();
  Eigen::Vector3f v1(q1.x(), q1.y(), q1.z());

  float s2 = q2.w();
  Eigen::Vector3f v2(q2.x(), q2.y(), q2.z());

  float cr_v1 = v1(1)*v2(2) - v1(2)*v2(1);
  float cr_v2 = v1(2)*v2(0) - v1(0)*v2(2);
  float cr_v3 = v1(0)*v2(1) - v1(1)*v2(0);

  Eigen::Quaternionf qp;
  qp.w() = s1*s2 - v2.dot(v1);
  qp.x() = v2(0)*s1 + s2*v1(0) + cr_v1;
  qp.y() = v2(1)*s1 + s2*v1(1) + cr_v2;
  qp.z() = v2(2)*s1 + s2*v1(2) + cr_v3;

  return qp;
}

Eigen::Quaterniond QuatMTimes(const Eigen::Quaterniond &q1,
    const Eigen::Quaterniond &q2)  {
  double s1 = q1.w();
  Eigen::Vector3d v1(q1.x(), q1.y(), q1.z());

  double s2 = q2.w();
  Eigen::Vector3d v2(q2.x(), q2.y(), q2.z());

  double cr_v1 = v1(1)*v2(2) - v1(2)*v2(1);
  double cr_v2 = v1(2)*v2(0) - v1(0)*v2(2);
  double cr_v3 = v1(0)*v2(1) - v1(1)*v2(0);

  Eigen::Quaterniond qp;
  qp.w() = s1*s2 - v2.dot(v1);
  qp.x() = v2(0)*s1 + s2*v1(0) + cr_v1;
  qp.y() = v2(1)*s1 + s2*v1(1) + cr_v2;
  qp.z() = v2(2)*s1 + s2*v1(2) + cr_v3;

  return qp;
}

float angBTquat(Eigen::Quaternionf &q1, Eigen::Quaternionf &q2) {
  q1.normalize();
  q2.normalize();

  Eigen::Quaternionf q_ = QuatMTimes(q1.inverse(), q2);

  float ang = 2.0f*acos(q_.w()); // acos: [0, pi]

  if (ang > PIf){
    ang = 2.0f*PIf - ang;
  }
  return fabs(ang);
}

double angBTquat(Eigen::Quaterniond &q1, Eigen::Quaterniond &q2) {
  q1.normalize();
  q2.normalize();

  Eigen::Quaterniond q_ = QuatMTimes(q1.inverse(), q2);

  double ang = 2.0*acos(q_.w()); // acos: [0, pi]

  if (ang > PI){
    ang = 2.0*PI - ang;
  }
  return fabs(ang);
}

double angBTVec(Eigen::Vector3d x, Eigen::Vector3d b,
    Eigen::Vector3d z, bool nonnegative) {
  x.normalize();
  b.normalize();
  if (z.norm() < 1e-5) {
    return acos(x.dot(b));
  } else {
    z.normalize();
    double ang = atan2(x.cross(b).dot(z), x.dot(b));
    if (nonnegative) ang = (ang < 0)? 2*PI + ang : ang;
    return ang;
  }
}

// SE(3)
CartesianPose::CartesianPose(std::vector<double> pose) {
  p_[0] = pose[0];
  p_[1] = pose[1];
  p_[2] = pose[2];
  q_.w() = pose[3];
  q_.x() = pose[4];
  q_.y() = pose[5];
  q_.z() = pose[6];
  R_ = q_.toRotationMatrix();
}

CartesianPose::CartesianPose(Eigen::Matrix4d T) {
  p_ = T.block<3,1>(0,3);
  R_ = T.block<3,3>(0,0);
  q_ = Eigen::Quaterniond(R_);
}

CartesianPose::CartesianPose(const Eigen::Quaterniond &q, const Eigen::Vector3d &p) {
  p_ = p;
  q_ = q;
  R_ = q_.toRotationMatrix();
}

void CartesianPose::setQuaternion(const Eigen::Quaterniond &q) {
  q_ = q;
  R_ = q_.toRotationMatrix();
}

void CartesianPose::setQuaternion(const std::vector<double> &q) {
  q_.w() = q[0];
  q_.x() = q[1];
  q_.y() = q[2];
  q_.z() = q[3];
  R_ = q_.toRotationMatrix();
}

void CartesianPose::setXYZ(const Eigen::Vector3d &p) {
  p_ = p;
}

void CartesianPose::setXYZ(const std::vector<double> &p) {
  p_[0] = p[0];
  p_[1] = p[1];
  p_[2] = p[2];
}

Eigen::Matrix3d CartesianPose::getRotationMatrix() const {
  return R_;
}

Eigen::Quaterniond CartesianPose::getQuaternion() const {
  return q_;
}

Eigen::Vector3d CartesianPose::getXYZ() const {
  return p_;
}

Eigen::Vector3d CartesianPose::getXAxis() const {
  return R_.block<3,1>(0, 0);
}

Eigen::Vector3d CartesianPose::getYAxis() const {
  return R_.block<3,1>(0, 1);
}

Eigen::Vector3d CartesianPose::getZAxis() const {
  return R_.block<3,1>(0, 2);
}


Eigen::Matrix4d CartesianPose::getTransformMatrix() const{
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3, 3>(0, 0) = R_;
  T.block<3, 1>(0, 3) = p_;
  return T;
}

Eigen::Isometry3d CartesianPose::getIsometry3d() const {
  Eigen::Isometry3d transform = Eigen::Translation<double, 3>(
      p_[0], p_[1], p_[2]) * q_;
  return transform;
}

CartesianPose CartesianPose::operator*(const CartesianPose &pose) const {
  Eigen::Matrix4d T = getTransformMatrix()*pose.getTransformMatrix();
  return CartesianPose(T);
}

CartesianPose CartesianPose::inv() const {
  Eigen::Matrix4d T = Eigen::Matrix4d::Identity();
  T.block<3, 3>(0, 0) = R_.transpose();
  T.block<3, 1>(0, 3) = -R_.transpose() * p_;
  return CartesianPose(T);
}

Eigen::Vector3d CartesianPose::transformVec(const Eigen::Vector3d &v) const {
  return R_ * v;
}

Eigen::Vector3d CartesianPose::transformPoint(const Eigen::Vector3d &p) const {
  return R_*p + p_;
}

Eigen::Quaterniond CartesianPose::transformQuat(const Eigen::Quaterniond &q) const {
  Eigen::Matrix3d R = q.toRotationMatrix();
  return Eigen::Quaterniond(R_*R);
}


void CartesianPose::print() const{
  std::cout << "p:\n";
  std::cout << p_ << std::endl;
  std::cout << "R:\n";
  std::cout << R_ << std::endl;
  std::cout << "q:\n";
  std::cout << q_.w() << " " << q_.x() << " "
      << q_.y() << " " << q_.z() << std::endl;
}


Eigen::MatrixXd transformByRAndP(const Eigen::MatrixXd &points_rowwise, const Eigen::Matrix3d &R, const Eigen::Vector3d &p) {
  Eigen::MatrixXd points_transformed = R*points_rowwise.transpose() + p.replicate(1, points_rowwise.rows());
  return points_transformed.transpose();
}

void double2float(const double *array_in, float *array_out, int n) {
  for (int i = 0; i < n; ++i)
    array_out[i] = array_in[i];
}
void float2double(const float *array_in, double *array_out, int n) {
  for (int i = 0; i < n; ++i)
    array_out[i] = array_in[i];
}


/////////////////////////////////////////////////////////////////////////
//                      Motion Planning
/////////////////////////////////////////////////////////////////////////
void MotionPlanningLinear(const double *pose0, const double *pose_set, const int Nsteps,
    MatrixXd *pose_traj) {
  Eigen::Quaternionf q0(pose0[3], pose0[4], pose0[5],pose0[6]);
  Eigen::Quaternionf qset(pose_set[3], pose_set[4], pose_set[5],pose_set[6]);
  Eigen::Quaternionf q;

  pose_traj->resize(7, Nsteps);
  for (int i = 0; i < Nsteps; ++i) {
    q = q0.slerp(double(i+1)/double(Nsteps), qset);
    (*pose_traj)(0, i) = (pose0[0]*double(Nsteps-i-1) + pose_set[0]*double(i+1))/double(Nsteps);
    (*pose_traj)(1, i) = (pose0[1]*double(Nsteps-i-1) + pose_set[1]*double(i+1))/double(Nsteps);
    (*pose_traj)(2, i) = (pose0[2]*double(Nsteps-i-1) + pose_set[2]*double(i+1))/double(Nsteps);
    (*pose_traj)(3, i) = q.w();
    (*pose_traj)(4, i) = q.x();
    (*pose_traj)(5, i) = q.y();
    (*pose_traj)(6, i) = q.z();
  }
}

void TrapezodialTrajectory(double x_f, double a_max, double v_max, double *t1,
    double *t2, int Nsteps, double * x_traj) {
  assert(x_f > -1e-7);
  assert(a_max > 0);
  assert(v_max > 0);
  double delta_x1 = v_max*v_max/2.0/a_max;
  if (2.0*delta_x1 > x_f) {
    double v_max_actual = std::sqrt(x_f*a_max);
    *t1 = v_max_actual/a_max;
    *t2 = 0;
  } else {
    *t1 = v_max/a_max;
    double delta_x2 = x_f - 2*delta_x1;
    *t2 = delta_x2/v_max;
  }
  if (Nsteps == 0) return;

  double t_all = 2*(*t1) + (*t2);
  for (int i = 0; i < Nsteps; ++i) {
    double ratio = double(i)/double(Nsteps-1); // [0, 1]
    double t = ratio*t_all;
    if (t < *t1) {
      x_traj[i] = 0.5*a_max*t*t;
    } else if (t < *t1 + *t2) {
      x_traj[i] = 0.5*(*t1)*v_max + v_max*(t-(*t1));
    } else {
      x_traj[i] = x_f - 0.5*a_max*(t_all - t)*(t_all - t);
    }
  } // end for
}

void MotionPlanningTrapezodial(const double *pose0, const double *pose_set,
    double a_max_trans, double v_max_trans, double a_max_rot, double v_max_rot,
    double rate, MatrixXd *pose_traj) {
  Eigen::Vector3d p0, pf;
  p0 << pose0[0], pose0[1], pose0[2];
  pf << pose_set[0], pose_set[1], pose_set[2];

  Eigen::Quaterniond q0(pose0[3], pose0[4], pose0[5], pose0[6]);
  Eigen::Quaterniond qf(pose_set[3], pose_set[4], pose_set[5],pose_set[6]);
  double dist_trans = (p0 - pf).norm();
  double dist_rot = angBTquat(q0, qf);

  double t1_trans = 0;
  double t2_trans = 0;
  double t1_rot = 0;
  double t2_rot = 0;
  TrapezodialTrajectory(dist_trans, a_max_trans, v_max_trans,
      &t1_trans, &t2_trans);
  TrapezodialTrajectory(dist_rot, a_max_rot, v_max_rot,
      &t1_rot, &t2_rot);
  int Nsteps = (int)std::round(std::max(2.0*t1_trans + t2_trans,
      2.0*t1_rot + t2_rot)*rate);

  if (Nsteps > 2) {
    // need to do interpolation
    // get ratio
    double *r_trans = new double[Nsteps];
    double *r_rot = new double[Nsteps];
    TrapezodialTrajectory(dist_trans, a_max_trans, v_max_trans,
        &t1_trans, &t2_trans, Nsteps, r_trans);
    TrapezodialTrajectory(dist_rot, a_max_rot, v_max_rot,
        &t1_rot, &t2_rot, Nsteps, r_rot);

    // regularize r_trans and r_rot to [0, 1]
    // note here Nsteps > 2
    if (dist_trans > 1e-5)
      for (int i = 0; i < Nsteps; ++i)
        r_trans[i] /= dist_trans;
    else
      for (int i = 0; i < Nsteps; ++i)
        r_trans[i] = i/double(Nsteps-1);

    if (dist_rot > 1e-5)
      for (int i = 0; i < Nsteps; ++i)
        r_rot[i] /= dist_rot;
    else
      for (int i = 0; i < Nsteps; ++i)
        r_rot[i] = i/double(Nsteps-1);

    // compute the actual pose from r
    Eigen::Quaterniond q;
    pose_traj->resize(7, Nsteps);
    for (int i = 0; i < Nsteps; ++i) {
      q = q0.slerp(r_rot[i], qf);
      (*pose_traj)(0, i) = pose0[0]*(1.0 - r_trans[i]) + pose_set[0]*r_trans[i];
      (*pose_traj)(1, i) = pose0[1]*(1.0 - r_trans[i]) + pose_set[1]*r_trans[i];
      (*pose_traj)(2, i) = pose0[2]*(1.0 - r_trans[i]) + pose_set[2]*r_trans[i];
      (*pose_traj)(3, i) = q.w();
      (*pose_traj)(4, i) = q.x();
      (*pose_traj)(5, i) = q.y();
      (*pose_traj)(6, i) = q.z();
    }
    delete [] r_trans;
    delete [] r_rot;
  } else if (Nsteps == 2) {
    // no need for interpolation
    pose_traj->resize(7, 2);
    (*pose_traj)(0, 0) = pose0[0];
    (*pose_traj)(1, 0) = pose0[1];
    (*pose_traj)(2, 0) = pose0[2];
    (*pose_traj)(3, 0) = pose0[3];
    (*pose_traj)(4, 0) = pose0[4];
    (*pose_traj)(5, 0) = pose0[5];
    (*pose_traj)(6, 0) = pose0[6];
    (*pose_traj)(0, 1) = pose_set[0];
    (*pose_traj)(1, 1) = pose_set[1];
    (*pose_traj)(2, 1) = pose_set[2];
    (*pose_traj)(3, 1) = pose_set[3];
    (*pose_traj)(4, 1) = pose_set[4];
    (*pose_traj)(5, 1) = pose_set[5];
    (*pose_traj)(6, 1) = pose_set[6];
  } else {
    // no need for interpolation
    pose_traj->resize(7, 1);
    (*pose_traj)(0, 0) = pose_set[0];
    (*pose_traj)(1, 0) = pose_set[1];
    (*pose_traj)(2, 0) = pose_set[2];
    (*pose_traj)(3, 0) = pose_set[3];
    (*pose_traj)(4, 0) = pose_set[4];
    (*pose_traj)(5, 0) = pose_set[5];
    (*pose_traj)(6, 0) = pose_set[6];
  }
}

double Gaussian(double x, double var) {
  double k = 1.0/var/std::sqrt(2.0*PI);
  double exp = std::exp(-x*x/2.0/var/var);
  return k*exp;
}

}
