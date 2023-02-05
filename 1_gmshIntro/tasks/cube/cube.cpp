#include <set>
#include <cmath>
#include <gmsh.h> 
#include <vector>
#include "../../utils/utils.cpp"




int main() {
    gmsh::initialize();
    gmsh::model::add("cube");
    double lc = 1e-1; 

    std::vector<int> down_points = {
        gmsh::model::geo::addPoint(0, 0, 0, lc),
        gmsh::model::geo::addPoint(1, 0, 0, lc),
        gmsh::model::geo::addPoint(1, 1, 0, lc),
        gmsh::model::geo::addPoint(0, 1, 0, lc)
    };  

    std::vector<int> up_points = {
        gmsh::model::geo::addPoint(0, 0, 1, lc),
        gmsh::model::geo::addPoint(1, 0, 1, lc),
        gmsh::model::geo::addPoint(1, 1, 1, lc),
        gmsh::model::geo::addPoint(0, 1, 1, lc)
    };

    std::vector<int> down_lines; 
    for (int i = 0; i < 4; i++) {
        down_lines.push_back(
            gmsh::model::geo::addLine(down_points[i], down_points[(i+1)%4])
        );
    }

    std::vector<int> up_lines; 
    for (int i = 0; i < 4; i++) {
        up_lines.push_back(
            gmsh::model::geo::addLine(up_points[i], up_points[(i+1)%4])
        );
    }

    std::vector<int> side_lines; 
    for (int i = 0; i <  4; i++) {
        side_lines.push_back(
            gmsh::model::geo::addLine(down_points[i], up_points[i])
        );
    }

    std::vector<int> surfaces; 
    surfaces.push_back(
        utils::_make_rectangle_surface(down_lines[0], down_lines[1], down_lines[2], down_lines[3])
    );
    surfaces.push_back(
        utils::_make_rectangle_surface(up_lines[0], up_lines[1], up_lines[2], up_lines[3])
    );
    for (int i = 0; i < 4; i++) {
        surfaces.push_back(
            utils::_make_rectangle_surface(down_lines[i], side_lines[(i+1)%4], -up_lines[i], -side_lines[i])
        );
    }

    gmsh::model::geo::addVolume(
        {
            gmsh::model::geo::addSurfaceLoop(
                {surfaces[0], surfaces[1], surfaces[2], surfaces[3],  surfaces[4],  surfaces[5]}
            )
        }
    );

    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(3);

    gmsh::fltk::run();

    gmsh::finalize();

    return 0; 
}