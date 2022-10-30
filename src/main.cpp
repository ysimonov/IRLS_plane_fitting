#include "irls_plane_fitting.hpp"

int main(int argc, char *argv[])
{

    Eigen::MatrixX3f points_xyz{{1, 3, 7}, {9, 5, 8}, {4, 10, 15}};
    PlaneFittingIRLS plane_fit(points_xyz);
    plane_fit.fitPlane();

    return 0;
}