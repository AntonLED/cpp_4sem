cmake_minimum_required(VERSION 3.20)
project(swarmit)

set(SOURCES
        ./src/Swarm.cpp
        ./src/main.cpp
)

set(HEADERS
        ./src/Swarm.h
)

add_executable(swarmit ${SOURCES} ${HEADERS})

find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

target_link_libraries(swarmit ${VTK_LIBRARIES})