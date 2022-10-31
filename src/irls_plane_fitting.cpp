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
    size_t number_of_points = this->getNumberOfPoints();

    // compute centroid vector
    Eigen::RowVector3f centroid = points_xyz.colwise().mean();

    std::cout << "Centroid: " << std::endl << centroid << std::endl;

    // compute covariance matrix
    Eigen::MatrixX3f centered = points_xyz.rowwise() - points_xyz.colwise().mean();
    Eigen::Matrix3f covariance_matrix = (centered.adjoint() * centered) / static_cast<double>(getNumberOfPoints() - 1);

    std::cout << "Covariance Matrix: " << std::endl << covariance_matrix << std::endl;

    // compute eigenvalues and eigenvectors of the covariance matrix
    // since covariance matrix is symmetric, it can only have real eigenvalues
    Eigen::EigenSolver<Eigen::Matrix3f> eigensolver;
    eigensolver.compute(covariance_matrix, /* computeEigenvectors */ true);
    Eigen::Vector3f eigenvalues = eigensolver.eigenvalues().real();
    Eigen::Matrix3f eigenvectors = eigensolver.eigenvectors().real();

    std::cout << "Eigenvalues: " << std::endl << eigenvalues << std::endl;
    std::cout << "Eigenvectors: " << std::endl << eigenvectors << std::endl;
}