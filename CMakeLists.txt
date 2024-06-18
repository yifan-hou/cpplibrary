cmake_minimum_required(VERSION 3.10)
project(RobotUtilities)
set(CMAKE_BUILD_TYPE Release)

add_compile_options(-std=c++11)

include_directories(include/)
# include_directories("${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}")

find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})

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
target_link_libraries(test_utilities TimerLinux Utilities)

#############
## Install ##
#############
# headers

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_PREFIX}/include/
    )

# binaries
install(TARGETS TimerLinux Utilities
    DESTINATION ${CMAKE_INSTALL_PREFIX}/lib/${PROJECT_NAME}/
    )