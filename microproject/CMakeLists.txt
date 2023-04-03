cmake_minimum_required(VERSION 3.0)
project(microproject)

add_executable(test ./src/test.cpp)

# find_package(glm REQUIRED)
find_package(glfw3 REQUIRED)

# link_libraries(glm)
target_link_libraries(test glfw3)
