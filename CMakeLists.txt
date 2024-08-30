# At LEAST 3.2 but newer is better
cmake_minimum_required(VERSION 3.9 FATAL_ERROR)

project(ParRidgeExtraction VERSION 0.1 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Export compilation database to use with auto-complete tools
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Set default build type to Release if not specified
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# Requirements: Eigen3 and OpenMP library
find_package(Eigen3 REQUIRED)
find_package(OpenMP REQUIRED)

# Create a static library for the common code
add_library(
  common_objects
  STATIC
  "${CMAKE_SOURCE_DIR}/src/experiments.cpp"
  "${CMAKE_SOURCE_DIR}/src/setup.cpp"
  "${CMAKE_SOURCE_DIR}/src/utils.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/marchingridges.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/particlesystem.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/datasets/abc.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/datasets/amiraflow.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/datasets/doublegyre.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/datasets/forcedduffing.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/datasets/movingcircleridges.cpp"
  "${CMAKE_SOURCE_DIR}/src/flowsim/datasets/tornado.cpp"
  "${CMAKE_SOURCE_DIR}/src/images/colormap.cpp"
  "${CMAKE_SOURCE_DIR}/src/images/texture.cpp"
  "${CMAKE_SOURCE_DIR}/src/images/visualization.cpp"
)

# Include directories and compile options for the common library
target_include_directories(
  common_objects
  PUBLIC "${CMAKE_SOURCE_DIR}/src"
  PUBLIC "${CMAKE_SOURCE_DIR}/external/inc"
)
target_link_libraries(
  common_objects
  PUBLIC stdc++fs
  PUBLIC Eigen3::Eigen
  PUBLIC OpenMP::OpenMP_CXX
  PUBLIC nlohmann_json::nlohmann_json
)

# Define executables and link them to the common static library
add_executable(execute "${CMAKE_SOURCE_DIR}/src/main.cpp")
add_executable(test-ridges "${CMAKE_SOURCE_DIR}/src/test-ridges.cpp")
target_link_libraries(execute PRIVATE common_objects)
target_link_libraries(test-ridges PRIVATE common_objects)

# Add subdirectory for the JSON library
add_subdirectory(external/json)

# Additional compiler flags
set(CMAKE_CXX_FLAGS "-Wall -Wextra -Wpedantic -march=native -Wno-c99-extensions")