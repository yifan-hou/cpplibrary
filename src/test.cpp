#include "RobotUtilities/utilities.h"
#include <Eigen/Geometry>

using namespace RUT;

int main() {
  // Eigen::Vector3d p1, p2;
  // p1 << 0, 0, 0;
  // p2 << 1, 0, 0;
  // Eigen::Vector3d n1, n2;
  // n1 << 0, 0, 1;
  // n2 << 0, -1, 1;
  // RUT::Vector6d line1, line2;
  // line1 = RUT::getPluckerLine(p1, n1);
  // line2 = RUT::getPluckerLine(p2, n2);
  // double dist = RUT::distBTPluckerLines(line1, line2);
  // double ang = RUT::angleBTPluckerLines(line1, line2);
  // std::cout << "line1: " << line1.transpose() << std::endl;
  // std::cout << "line2: " << line2.transpose() << std::endl;
  // std::cout << "dist = " << dist << std::endl;
  // std::cout << "ang (deg) = " << ang*180.0/3.1415926 << std::endl;
  // Eigen::Vector3d p3;
  // p3 << 1, 1, 0;
  // double dist1 = distPoint2PluckerLine(p3, line1);
  // double dist2 = distPoint2PluckerLine(p3, line2);
  // std::cout << "point dist to line 1 = " << dist1 << std::endl;
  // std::cout << "point dist to line 2 = " << dist2 << std::endl;

  // std::vector<double> pose_vec = {1,2,3,0,1,0,0};
  // CartesianPose pose = CartesianPose(pose_vec);
  // Eigen::Isometry3d iso = pose.getIsometry3d();
  // CartesianPose pose2 = CartesianPose(iso);
  // std::cout << "Pose 1: " << pose.poseString() << std::endl;
  // std::cout << "Pose 2: " << pose2.poseString() << std::endl;
  // pose2.print();

  // Eigen::Quaterniond q = Eigen::Quaterniond(0.1, 0.9, 0.9, 0.1);
  // MatrixXd R = q.toRotationMatrix();
  // std::cout << "q: " << q.w() << ", " << q.x() << ", " << q.y() << ", " << q.z() << std::endl;
  // std::cout << "R:\n" << R << std::endl;
  // std::cout << "R'*R:\n" << R.transpose()*R << std::endl;
  // std::cout << "R norm: " << R.norm() << std::endl;

  std::vector<double> pose_vec1 = {0,0,0,1,0,0,0};
  std::vector<double> pose_vec2 = {0.405014, -0.348191, 0.330936, -0.017058, 0.998968, -0.019972, 0.037051};
  CartesianPose pose_WGs = CartesianPose(pose_vec1);
  CartesianPose pose_WT = CartesianPose(pose_vec2);
  CartesianPose pose_TGs = pose_WT.inv()*pose_WGs;
  std::cout << "pose_WGs: " << pose_WGs.poseString() << std::endl;
  std::cout << "pose_WT: " << pose_WT.poseString() << std::endl;
  std::cout << "pose_WT.inv(): " << pose_WT.inv().poseString() << std::endl;
  std::cout << "pose_TGs: " << pose_TGs.poseString() << std::endl;

  // Eigen::MatrixXd A(2,5);
  // A.setRandom();
  // // A << -4.3102,   14.3674,   -3.3332,
  // //            0,         0,         0,
  // //       4.3102,   14.3674,   -0.8333;
  // std::cout << "A: \n" << A << std::endl;
  // Eigen::MatrixXd Arowspace = A;
  // int r = rowSpace(&Arowspace, 1e-10);
  // std::cout << "rank: " << r << std::endl;
  // std::cout << "rowSpace(A): \n" << Arowspace << std::endl;
  // std::cout << "Arowspace*(Arowspace'): \n" << Arowspace*Arowspace.transpose() << std::endl;
  // Eigen::MatrixXd nullA;
  // Arowspace = A;
  // r = nullSpace(&Arowspace, &nullA);
  // std::cout << "rank: " << r << std::endl;
  // std::cout << "NullSpace(A): \n" << nullA << std::endl;
  // std::cout << "nullA*(nullA'): \n" << nullA*nullA.transpose() << std::endl;
  // std::cout << "A*NullA: \n" << A*nullA.transpose() << std::endl;

  /* test sampling */
  // int size = 2;
  // Eigen::MatrixXd covar(size,size);
  // covar << 1, .5,
  //         .5, 1;
  // normal_random_variable sample { covar };
  // std::cout << sample() << std::endl;
  // std::cout << sample() << std::endl;

  /* test wedgeright*/
  // Eigen::VectorXd v = Eigen::VectorXd::Random(6);
  // Eigen::VectorXd p = Eigen::VectorXd::Random(4);
  // std::cout << "wedge6(v)*p:\n" << wedge6(v)*p << std::endl;
  // std::cout << "wedgeRight6(p) * v:\n" << wedgeRight6(p) * v << std::endl;
  return 0;
}