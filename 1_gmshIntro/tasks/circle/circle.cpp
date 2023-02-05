#include <gmsh.h>
#include <set>
#include <iostream>
#include "../../utils/utils.cpp"



int main(int argc, char **argv) {
    gmsh::initialize();
    gmsh::model::add("test1");
    double lc = 1e-1;

    // auto p1 = gmsh::model::geo::addPoint(
    //     0, 0, 0, lc, 1
    // );
    // auto p2 = gmsh::model::geo::addPoint(
    //     5, 0, 0, lc, 2
    // );
    // auto p3 = gmsh::model::geo::addPoint(
    //     10, 0, 0, lc, 3
    // );


    // // std::vector<int> circ = utils::_make_circle_surface(p1, p2, p3);


    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);

    gmsh::fltk::run();
    gmsh::finalize();

    return 0;
}