cmake_minimum_required(VERSION 3.6 FATAL_ERROR)
project(irls_plane_fitting)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 11)
endif()

# Default to C++17
set(CMAKE_CXX_STANDARD 17)

# Compile options
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# include directories
include_directories(
    ${PROJECT_SOURCE_DIR}/include
    ${PROJECT_SOURCE_DIR}/src
)

file(GLOB SOURCES
    "${PROJECT_SOURCE_DIR}/include/*.h"
    "${PROJECT_SOURCE_DIR}/include/*.hpp"
    "${PROJECT_SOURCE_DIR}/src/*.cpp"
    "${PROJECT_SOURCE_DIR}/src/*.c"
)

# Eigen
find_package(Eigen3 3.4 REQUIRED NO_MODULE)
link_directories(${Eigen_INCLUDE_DIRS})

# # VTK
# find_package(VTK 9.2 REQUIRED)

# # PCL
# find_package(PCL 1.7.2 REQUIRED)
# include_directories(${PCL_INCLUDE_DIRS})
# link_directories(${PCL_LIBRARY_DIRS})
# add_definitions(${PCL_DEFINITIONS})

# executable and build
add_executable(${PROJECT_NAME} ${SOURCES})

# linking
target_link_libraries(${PROJECT_NAME} ${Eigen_LIBRARIES})
# target_link_libraries(${PROJECT_NAME} ${VTK_LIBRARIES})
# target_link_libraries(${PROJECT_NAME} ${PCL_LIBRARIES})

