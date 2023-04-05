cmake_minimum_required(VERSION 3.20)
project(microproject)

set(SOURCES
        ./src/Swarm.cpp
)

set(HEADERS
        ./src/Swarm.h
)

add_executable(microproject ${SOURCES} ${HEADERS})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

target_link_libraries(microproject ${VTK_LIBRARIES})