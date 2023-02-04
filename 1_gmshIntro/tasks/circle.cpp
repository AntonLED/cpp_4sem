#include <gmsh.h>
#include <set>
#include <iostream>

int main(int argc, char **argv)
{
    gmsh::initialize();
    gmsh::model::add("test1");
    double lc = 1e-1;

    auto p1 = gmsh::model::geo::addPoint(
        0, 0, 0, lc, 1
    );
    auto p2 = gmsh::model::geo::addPoint(
        5, 0, 0, lc, 2
    );
    auto p3 = gmsh::model::geo::addPoint(
        10, 0, 0, lc, 3
    );

    auto c1 = gmsh::model::geo::addCircleArc(
        p1, p2, p3
    );
    auto c2 = gmsh::model::geo::addCircleArc(
        p3, p2, p1
    );

    auto cl1 = gmsh::model::geo::addCurveLoop(
        {c1, c2}, 1
    );
    auto serf1 = gmsh::model::geo::addPlaneSurface(
        {cl1}, 1
    );


    gmsh::model::geo::synchronize();
    gmsh::model::mesh::generate(2);

    gmsh::fltk::run();
    gmsh::finalize();

    return 0;
}