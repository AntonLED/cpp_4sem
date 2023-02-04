#!/bin/bash
# bash script for more convenient building gmsh applications under UNIX OS
g++ -o executable $1 -lgmsh
./executable
rm ./executable