cmake_minimum_required(VERSION 2.8.12)
project(RobotUtilities)
set(CMAKE_BUILD_TYPE Release)

add_compile_options(-std=c++11)

include_directories(include/)
# include_directories("${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

find_package(libInterpolate REQUIRED)
add_executable( interpolation2 src/interpolation2.cpp )
target_link_libraries( interpolation2 libInterpolate::Interpolate )

###########
## Build ##
###########
add_library(TimerLinux SHARED
	src/TimerLinux.cpp
)

add_library(Utilities SHARED
	src/utilities.cpp
)

add_executable(test_utilities src/test.cpp)
target_link_libraries(test_utilities Utilities)

add_executable(interpolation src/interpolation.cpp)
#############
## Install ##
#############
# headers

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
    )

# binaries
install(TARGETS TimerLinux Utilities
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/
    )