#include <Eigen/Geometry>
#include <iostream>
#include <fstream>
#include <vector>
#include <libInterpolate/Interpolate.hpp>

// https://github.com/CD3/libInterpolate

#define KOMO_DT 0.01
#define UR_DT 0.002

int main() {
  Eigen::MatrixXd traj(100, 6);
  traj <<
       1.29937, -0.91983,  1.73309, 0.796842,  1.58584, 0.720023,
       1.29937, -0.91983,  1.73309, 0.796842,  1.58584, 0.720023,
       1.29937, -0.91983,  1.73309, 0.796842,  1.58584, 0.720023,
       1.29937, -0.91983,  1.73309, 0.796842,  1.58584, 0.720023,
       1.29937, -0.91983,  1.73309, 0.796842,  1.58584, 0.720023,
       1.29937, -0.91983,  1.73309, 0.796842,  1.58584, 0.720023,
       1.29935, -0.919863, 1.73308, 0.796849,  1.58585, 0.720023,
       1.29921, -0.920097, 1.733  , 0.796899,  1.58599, 0.720023,
       1.29883, -0.920733, 1.73279, 0.797033,  1.58636, 0.720023,
       1.2981 , -0.921971, 1.73238, 0.797294,  1.58707, 0.720023,
       1.2969 , -0.924011, 1.73171, 0.797724,  1.58825, 0.720023,
       1.29511, -0.927025, 1.73072, 0.79836 ,  1.58999, 0.720023,
       1.29272, -0.931062, 1.72939, 0.799217,  1.59234, 0.720023,
       1.28968, -0.936143, 1.72775, 0.800303,  1.5953 , 0.720023,
       1.28598, -0.942288, 1.72577, 0.801629,  1.5989 , 0.720023,
       1.2816 , -0.949517, 1.72348, 0.803202,  1.60316, 0.720023,
       1.27652, -0.957837, 1.72087, 0.805029,  1.60809, 0.720023,
       1.27075, -0.9672  , 1.71798, 0.807102,  1.61368, 0.720023,
       1.26431, -0.977548, 1.71481, 0.80941 ,  1.6199 , 0.720023,
       1.2572 , -0.988819, 1.71141, 0.811942,  1.62672, 0.720023,
       1.24946, -1.00095,  1.70779, 0.814686,  1.63414, 0.720023,
       1.2411 , -1.01389,  1.70399, 0.81763 ,  1.64211, 0.720023,
       1.23214, -1.02753,  1.70004, 0.820756,  1.65062, 0.720023,
       1.22261, -1.0418 ,  1.69597, 0.824046,  1.65961, 0.720023,
       1.21254, -1.0566 ,  1.69182, 0.827481,  1.66905, 0.720023,
       1.20195, -1.07184,  1.68763, 0.831043,  1.67891, 0.720023,
       1.19088, -1.08743,  1.68343, 0.834712,  1.68915, 0.720023,
       1.17934, -1.10328,  1.67926, 0.838471,  1.69973, 0.720023,
       1.16736, -1.11928,  1.67517, 0.842304,  1.71061, 0.720023,
       1.15497, -1.13534,  1.67119, 0.846193,  1.72177, 0.720023,
       1.14219, -1.15136,  1.66736, 0.850121,  1.73315, 0.720023,
       1.12905, -1.16725,  1.66374, 0.854071,  1.74473, 0.720023,
       1.11557, -1.18292,  1.66035, 0.858033,  1.75648, 0.720023,
       1.10176, -1.19828,  1.65726, 0.861995,  1.76836, 0.720023,
       1.08765, -1.21324,  1.65449, 0.865947,  1.78035, 0.720023,
       1.07325, -1.22771,  1.65209, 0.869877,  1.79241, 0.720023,
       1.05859, -1.24161,  1.65011, 0.873776,  1.80451, 0.720023,
       1.04367, -1.25487,  1.64859, 0.877638,  1.81664, 0.720023,
       1.02852, -1.26744,  1.64754, 0.881459,  1.82877, 0.720023,
       1.01314, -1.27924,  1.64702, 0.885235,  1.84089, 0.720023,
      0.997556, -1.29022,  1.64704, 0.888959,  1.85296, 0.720023,
      0.981775, -1.30033,  1.64765, 0.89263 ,  1.86498, 0.720023,
      0.965813, -1.30956,  1.64883, 0.896244,  1.87694, 0.720023,
      0.949681, -1.31792,  1.65057, 0.899798,  1.88882, 0.720023,
      0.933391, -1.32541,  1.65287, 0.903292,  1.90062, 0.720023,
      0.916956, -1.33204,  1.65571, 0.906722,  1.91233, 0.720023,
      0.900388, -1.33782,  1.65909, 0.910087,  1.92396, 0.720023,
      0.883698, -1.34278,  1.66298, 0.913386,  1.93549, 0.720023,
      0.866901, -1.34694,  1.66736, 0.916618,  1.94693, 0.720023,
      0.850007, -1.35033,  1.67219, 0.919782,  1.95827, 0.720023,
       0.83303,   -1.353,  1.67747, 0.922877,  1.96951, 0.720023,
      0.815982, -1.35496,  1.68315, 0.925902,  1.98065, 0.720023,
      0.798878, -1.35625,  1.68921, 0.928858,  1.99169, 0.720023,
      0.781731, -1.35691,  1.69563, 0.931743,  2.00262, 0.720023,
      0.764558, -1.35698,  1.70238, 0.934557,  2.01345, 0.720023,
      0.747371, -1.35649,  1.70942, 0.937299,  2.02416, 0.720023,
      0.730187, -1.35548,  1.71673,  0.93997,  2.03476, 0.720023,
      0.713024, -1.35398,  1.72427, 0.942568,  2.04525, 0.720023,
      0.695899, -1.35203,  1.73203, 0.945094,  2.05561, 0.720023,
      0.678833, -1.34967,  1.73997, 0.947548,  2.06585, 0.720023,
      0.661843, -1.34693,  1.74807, 0.949929,  2.07595, 0.720023,
      0.644949, -1.34384,  1.75629, 0.952238,  2.08592, 0.720023,
      0.628174, -1.34044,  1.76461, 0.954473,  2.09574, 0.720023,
      0.611543, -1.33677,    1.773, 0.956636,  2.10541, 0.720023,
       0.59508, -1.33286,  1.78144, 0.958727,  2.11492, 0.720023,
       0.57881, -1.32873,  1.78989, 0.960745,  2.12425, 0.720023,
      0.562758, -1.32443,  1.79834, 0.962691,  2.13341, 0.720023,
      0.546954, -1.31998,  1.80674, 0.964564,  2.14237, 0.720023,
      0.531429, -1.31541,  1.81509, 0.966366,  2.15113, 0.720023,
      0.516214, -1.31076,  1.82334, 0.968097,  2.15967, 0.720023,
      0.501339, -1.30606,  1.83147, 0.969756,  2.16798, 0.720023,
      0.486837, -1.30134,  1.83947, 0.971345,  2.17604, 0.720023,
      0.472742, -1.29662,  1.84729, 0.972862,  2.18385, 0.720023,
      0.459091, -1.29194,  1.85492, 0.974308,  2.19138, 0.720023,
      0.445919, -1.28732,  1.86232, 0.975684,  2.19862, 0.720023,
      0.433263, -1.28279,  1.86947, 0.976989,  2.20555, 0.720023,
      0.421157, -1.27839,  1.87635, 0.978223,  2.21216, 0.720023,
      0.409637, -1.27413,  1.88293, 0.979385,  2.21843, 0.720023,
      0.398737, -1.27004,  1.88919, 0.980475,  2.22434, 0.720023,
      0.388491, -1.26614,   1.8951, 0.981493,  2.22989, 0.720023,
      0.378935, -1.26247,  1.90064, 0.982437,  2.23505, 0.720023,
      0.370098, -1.25903,  1.90578, 0.983308,  2.23981, 0.720023,
      0.362004, -1.25586,  1.91051, 0.984104,  2.24416, 0.720023,
      0.354668, -1.25295,  1.91482, 0.984825,  2.24809, 0.720023,
      0.348109, -1.25033,  1.91868, 0.985469,   2.2516, 0.720023,
      0.342344, -1.24801,   1.9221, 0.986037,  2.25468, 0.720023,
      0.337385, -1.24601,  1.92505, 0.986528,  2.25732, 0.720023,
      0.333208, -1.24431,  1.92754, 0.986944,  2.25955, 0.720023,
      0.329785, -1.24291,  1.92959, 0.987286,  2.26137, 0.720023,
      0.327087,  -1.2418,  1.93122, 0.987557,   2.2628, 0.720023,
      0.325085, -1.24097,  1.93242, 0.987758,  2.26386, 0.720023,
       0.32373, -1.24041,  1.93324, 0.987895,  2.26458, 0.720023,
      0.322909, -1.24008,  1.93374, 0.987978,  2.26501, 0.720023,
      0.322487,  -1.2399,  1.93399,  0.98802,  2.26524, 0.720023,
      0.322332, -1.23984,  1.93409, 0.988036,  2.26532, 0.720023,
      0.322309, -1.23983,   1.9341, 0.988038,  2.26533, 0.720023,
      0.322309, -1.23983,   1.9341, 0.988038,  2.26533, 0.720023,
      0.322309, -1.23983,   1.9341, 0.988038,  2.26533, 0.720023,
      0.322309, -1.23982,   1.9341, 0.988038,  2.26533, 0.720023,
      0.322309, -1.23982,   1.9341, 0.988038,  2.26533, 0.720023;

      std::vector<double> J1, J2, J3, J4, J5, J6, time;
      int N = traj.rows();
      for (int i = 0; i < N; ++i) {
        time.push_back(i*0.01);
        J1.push_back(traj(i,0));
        J2.push_back(traj(i,1));
        J3.push_back(traj(i,2));
        J4.push_back(traj(i,3));
        J5.push_back(traj(i,4));
        J6.push_back(traj(i,5));
      }

      // Use a cubic spline interpolator to interpolate the data
      _1D::CubicSplineInterpolator<double> interp1;
      _1D::CubicSplineInterpolator<double> interp2;
      _1D::CubicSplineInterpolator<double> interp3;
      _1D::CubicSplineInterpolator<double> interp4;
      _1D::CubicSplineInterpolator<double> interp5;
      _1D::CubicSplineInterpolator<double> interp6;
      interp1.setData(time, J1);
      interp2.setData(time, J2);
      interp3.setData(time, J3);
      interp4.setData(time, J4);
      interp5.setData(time, J5);
      interp6.setData(time, J6);

      int N_UR = (N-1)*KOMO_DT/UR_DT;

      std::ofstream fp;
      fp.open("test.yaml");
      if (!fp.is_open()) {
        std::cout << "[Error] Cannot open file" << std::endl;
      }

      fp << "actions:\n";

      for (int i = 0; i < N_UR; ++i) {
        double t = UR_DT*i;
        std::cout << "t: " << t << ",\tJ1: " << interp1(t)
                                << ",\tJ2: " << interp2(t)
                                << ",\tJ3: " << interp3(t)
                                << ",\tJ4: " << interp4(t)
                                << ",\tJ5: " << interp5(t)
                                << ",\tJ6: " << interp6(t)
                                << std::endl;
        fp << "  - joint_position:\n      q:\n";
        fp << "        - " <<  interp1(t) << "\n";
        fp << "        - " <<  interp2(t) << "\n";
        fp << "        - " <<  interp3(t) << "\n";
        fp << "        - " <<  interp4(t) << "\n";
        fp << "        - " <<  interp5(t) << "\n";
        fp << "        - " <<  interp6(t) << "\n";
      }

      fp.close();
  return 0;
}