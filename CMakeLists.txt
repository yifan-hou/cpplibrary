cmake_minimum_required(VERSION 2.8.12)
project(yifanlibrary)
set (CMAKE_CXX_FLAGS "--std=gnu++11 ${CMAKE_CXX_FLAGS}")

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}")



###########
## Build ##
###########
add_library(SINGLETON SHARED
	src/SINGLETON.cpp 
# 	"src/PracticalSocket.cpp" 
# 	"src/egm.pb.cc" 
# 	"src/TimerLinux.cpp"
)

add_library(TimerLinux SHARED
	src/TimerLinux.cpp
)

# target_link_libraries(SINGLETON	
#     ${Boost_FILESYSTEM_LIBRARY}
#     ${Boost_SYSTEM_LIBRARY}
#     ${PROTOBUF_LIBRARY}
# )


#############
## Install ##
#############
# headers
install(DIRECTORY include/
    DESTINATION /usr/local/include/
    )

# binaries
install(TARGETS SINGLETON TimerLinux
    DESTINATION /usr/local/lib/yifanlibrary
    )
