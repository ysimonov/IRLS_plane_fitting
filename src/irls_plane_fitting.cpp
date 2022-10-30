#include "irls_plane_fitting.hpp"

PlaneFittingIRLS::PlaneFittingIRLS(PlaneFittingIRLS &&rhs) noexcept
{
    points_ = std::move(rhs.points_);
    number_of_points_ = std::move(rhs.number_of_points_);
}

PlaneFittingIRLS &PlaneFittingIRLS::operator=(const PlaneFittingIRLS &rhs)
{
    if (this != &rhs)
    {
        points_ = rhs.points_;
        number_of_points_ = rhs.number_of_points_;
    }
    return *this;
}

PlaneFittingIRLS &PlaneFittingIRLS::operator=(PlaneFittingIRLS &&rhs) noexcept
{
    if (this != &rhs)
    {
        points_ = std::move(rhs.points_);
        number_of_points_ = std::move(rhs.number_of_points_);
    }
    return *this;
}

void PlaneFittingIRLS::fitPlane()
{
    // get points
    Eigen::MatrixX3f points_xyz = this->getPoints();

    // calculate centroid
    Eigen::Vector3f centroid = points_xyz.colwise().mean();

    std::cout << "Centroid: " << centroid << std::endl;

    // calculate covariance matrix
    Eigen::MatrixXf centered = points_xyz.rowwise() - points_xyz.colwise().mean();
    Eigen::MatrixXf covariance_matrix = (centered.adjoint() * centered) / static_cast<double>(getNumberOfPoints() - 1);

    std::cout << "Covariance Matrix: " << covariance_matrix << std::endl;
}