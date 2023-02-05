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

    int _make_circle_curve(int left_point, int bottom_point, int right_point, int top_point, int center_point) {
        return
            gmsh::model::geo::addCurveLoop(
                {
                    gmsh::model::geo::addCircleArc(
                        left_point, center_point, bottom_point
                    ), 
                    gmsh::model::geo::addCircleArc(
                        bottom_point, center_point, right_point
                    ),
                    gmsh::model::geo::addCircleArc(
                        right_point, center_point, top_point
                    ),
                    gmsh::model::geo::addCircleArc(
                        top_point, center_point, left_point
                    )
                }
            );
    }

    int _make_circle_surface(int left_point, int bottom_point, int right_point, int top_point, int center_point) { 
        return 
            gmsh::model::geo::addPlaneSurface(
                {
                    gmsh::model::geo::addCurveLoop(
                        {
                            gmsh::model::geo::addCircleArc(
                                left_point, center_point, bottom_point
                            ), 
                            gmsh::model::geo::addCircleArc(
                                bottom_point, center_point, right_point
                            ),
                            gmsh::model::geo::addCircleArc(
                                right_point, center_point, top_point
                            ),
                            gmsh::model::geo::addCircleArc(
                                top_point, center_point, left_point
                            )
                        }
                    )
                }
            );
    }
}
