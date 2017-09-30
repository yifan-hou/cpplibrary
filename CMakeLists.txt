cmake_minimum_required(VERSION 2.8.12)
project(yifanlibrary)

include_directories("${PROJECT_SOURCE_DIR}/src")
include_directories("${PROJECT_SOURCE_DIR}/include")
include_directories("${PROJECT_SOURCE_DIR}/include/${PROJECT_NAME}")



###########
## Build ##
###########
# add_library(egm_comm 
# 	"src/EGMClass.cpp" 
# 	"src/PracticalSocket.cpp" 
# 	"src/egm.pb.cc" 
# 	"src/TimerLinux.cpp")

# target_link_libraries(egm_comm	
#     ${Boost_FILESYSTEM_LIBRARY}
#     ${Boost_SYSTEM_LIBRARY}
#     ${PROTOBUF_LIBRARY}
# )


#############
## Install ##
#############
install(DIRECTORY include/
    DESTINATION /usr/local/include/
    )

