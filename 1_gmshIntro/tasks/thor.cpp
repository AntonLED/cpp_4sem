#include <gmsh.h>
#include <set>
#include <iostream>
#include "../utils/utils.cpp"



int main() {
    gmsh::initialize();
    gmsh::model::add("thor"); 
    double lc = 1; 
    
    auto p1 = gmsh::model::geo::addPoint(
        0 + 10, 0, 0, lc
    );
    auto p2 = gmsh::model::geo::addPoint(
        5 + 10, -5, 0, lc
    );
    auto p3 = gmsh::model::geo::addPoint(
        10 + 10, 0, 0, lc
    );
    auto p4 = gmsh::model::geo::addPoint(
        5 + 10, 5, 0, lc
    );
    auto p5 = gmsh::model::geo::addPoint(
        5 + 10, 0, 0, lc
    );

    auto p6 = gmsh::model::geo::addPoint(
        0 - 10, 0, 0, lc
    );
    auto p7 = gmsh::model::geo::addPoint(
        -5 - 10, -5, 0, lc
    );
    auto p8 = gmsh::model::geo::addPoint(
        -10 - 10, 0, 0, lc
    );
    auto p9 = gmsh::model::geo::addPoint(
        -5 - 10, 5, 0, lc
    );
    auto p10 = gmsh::model::geo::addPoint(
        -5 - 10, 0, 0, lc
    );
    
    int circ1 = utils::_make_circle_surface(p1, p2, p3, p4, p5);
    int circ2 = utils::_make_circle_surface(p6, p7, p8, p9, p10);

    std::vector<std::pair<int, int>> ov1; 
    gmsh::model::geo::revolve(
        {{2, circ1}},
        0, 0, 0,
        0, 1, 0, 
        M_PI / 2,
        ov1
    );
    std::vector<std::pair<int, int>> ov2; 
    gmsh::model::geo::revolve(
        {{2, circ1}},
        0, 0, 0,
        0, 1, 0, 
        -M_PI / 2,
        ov2
    );
    std::vector<std::pair<int, int>> ov3; 
    gmsh::model::geo::revolve(
        {{2, circ2}},
        0, 0, 0,
        0, 1, 0, 
        M_PI / 2,
        ov3
    );
    std::vector<std::pair<int, int>> ov4; 
    gmsh::model::geo::revolve(
        {{2, circ2}},
        0, 0, 0,
        0, 1, 0, 
        -M_PI / 2,
        ov4
    );

    gmsh::model::geo::addSurfaceLoop(
        {
            ov1[0].second,  ov2[0].second,  ov3[0].second,  ov4[0].second,
        }
    );

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(3);

    gmsh::fltk::run();
    gmsh::finalize();

    return 0; 
}