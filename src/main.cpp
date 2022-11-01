#include "irls_plane_fitting.hpp"

#include <algorithm> // for std::copy
#include <fstream>
#include <iostream>
#include <iterator>
#include <vector>

struct Point
{
    float x;
    float y;
    float z;
};

int main()
{

    // read data
    std::vector<Point> points;
    std::ifstream fin("../data/dataset.txt");
    while (fin.good())
    {
        Point temp;
        fin >> temp.x >> temp.y >> temp.z;
        points.push_back(temp);
        // std::cout << temp.x << " " << temp.y << " " << temp.z << std::endl;
    }

    // convert to Eigen
    size_t number_of_lines = points.size();
    Eigen::MatrixXf points_xyz(number_of_lines, 3);
    for (size_t i = 0; i < number_of_lines; ++i)
    {
        points_xyz(i, 0) = points[i].x;
        points_xyz(i, 1) = points[i].y;
        points_xyz(i, 2) = points[i].z;
    }

    // Eigen::MatrixX3f points_xyz{{1, 3, 7}, {9, 5, 8}, {4, 10, 15}, {1.5, 3.3, 4.5}, {1.2, 10.5, 7.7}};
    PlaneFittingIRLS plane_fit(points_xyz);
    plane_fit.fitPlane();

    return 0;
}