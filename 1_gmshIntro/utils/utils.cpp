#include <set>
#include <gmsh.h>

namespace utils {
    int _make_rectangle_surface(int line1, int line2, int line3, int line4) {
        return 
            gmsh::model::geo::addPlaneSurface(
                {
                    gmsh::model::geo::addCurveLoop(
                        {
                            line1, line2, line3, line4
                        }
                    )
                } 
            );
    }

    int _make_circle_surface(int start_point, int center_point, int end_point) { 
        return 
            gmsh::model::geo::addPlaneSurface(
                {
                    gmsh::model::geo::addCurveLoop(
                        {
                            gmsh::model::geo::addCircleArc(
                                start_point, center_point, end_point
                            ), 
                            gmsh::model::geo::addCircleArc(
                                end_point, center_point, start_point
                            )
                        }
                    )
                }
            );
    }
}
