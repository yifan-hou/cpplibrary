cmake_minimum_required(VERSION 2.8.12)
project(RobotUtilities)
set(CMAKE_BUILD_TYPE Release)

add_compile_options(-std=c++11)

include_directories(include/)
# include_directories("${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}")

include_directories(/usr/include/eigen3)


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