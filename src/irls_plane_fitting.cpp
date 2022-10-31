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

void PlaneFittingIRLS::fitPlane(size_t max_iterations, float k_welsh, float conv_thresh)
{
    // get points
    Eigen::MatrixX3f points_xyz = this->getPoints();
    size_t number_of_points = this->getNumberOfPoints();

    std::cout << "Original Points: " << std::endl << points_xyz << std::endl;

    // compute centroid vector
    Eigen::RowVector3f centroid = points_xyz.colwise().mean();

    // copy values
    Eigen::RowVector3f centroid_prev;
    centroid_prev = centroid;

    Eigen::RowVector3f centroid_curr;
    centroid_curr = centroid;

    std::cout << "Centroid: " << std::endl << centroid << std::endl;

    // compute covariance matrix
    Eigen::MatrixX3f centered = points_xyz.rowwise() - points_xyz.colwise().mean();
    Eigen::Matrix3f covariance_matrix = (centered.adjoint() * centered) / static_cast<double>(getNumberOfPoints() - 1);

    std::cout << "Covariance Matrix: " << std::endl << covariance_matrix << std::endl;

    // compute eigenvalues and eigenvectors of the covariance matrix
    // since covariance matrix is symmetric, it can only have real eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver;
    eigensolver.compute(covariance_matrix);
    Eigen::RowVector3f eigenvalues = eigensolver.eigenvalues().transpose();
    Eigen::Matrix3f eigenvectors = eigensolver.eigenvectors();

    std::cout << "Eigenvalues: " << std::endl << eigenvalues << std::endl;
    std::cout << "Eigenvectors: " << std::endl << eigenvectors << std::endl;

    // repeat minimization procedure for max_iterations
    Eigen::Vector3f normal_vector = eigenvectors(Eigen::all, 0);
    Eigen::Vector3f old_normal_vector = Eigen::Vector3f::Zero();
    for (size_t iter = 0; iter < max_iterations; ++iter)
    {
        // update old normal vector
        old_normal_vector = normal_vector;

        // center points
        Eigen::MatrixX3f points_xyz_centered = points_xyz.rowwise() - centroid;

        std::cout << "Centered Points: " << std::endl << points_xyz_centered << std::endl;

        // calculate distances to the plane (dot product between points and normal vector)
        Eigen::VectorXf distances = points_xyz_centered * normal_vector;

        std::cout << "Distances: " << std::endl << distances << std::endl;

        // calculate Welsh function
        Eigen::RowVectorXf weights = Eigen::exp(-(distances.array() / k_welsh).pow(2.0));

        std::cout << "Weights: " << std::endl << weights << std::endl;

        // recalculate centroid
        centroid_curr = weights * (points_xyz_centered.rowwise() - centroid_prev) / weights.sum();

        Eigen::MatrixX3f shifted_xyz = points_xyz_centered.rowwise() - centroid_curr;

        std::cout << "shifted_xyz: " << std::endl << shifted_xyz << std::endl;

        // 1xN * Nx3
        std::cout << "Size of weights: " << weights.rows() << "x" << weights.cols() << std::endl;
        std::cout << "Size of shifted_xyz: " << shifted_xyz.rows() << "x" << shifted_xyz.cols() << std::endl;

        Eigen::Matrix3Xf temp = shifted_xyz.transpose().array().rowwise() * weights.array();
        std::cout << "Size of temp: " << temp.rows() << "x" << temp.cols() << std::endl;

        // covariance matrix
        covariance_matrix = temp * shifted_xyz;
        std::cout << "Covariance Matrix (Updated): " << std::endl << covariance_matrix << std::endl;

        // recalculate eigenvalues and eigenvectors
        eigensolver.compute(covariance_matrix);
        Eigen::RowVector3f eigenvalues = eigensolver.eigenvalues().transpose();
        Eigen::Matrix3f eigenvectors = eigensolver.eigenvectors();

        std::cout << "Eigenvalues: " << std::endl << eigenvalues << std::endl;
        std::cout << "Eigenvectors: " << std::endl << eigenvectors << std::endl;

        // check convergence criteria
        normal_vector = eigenvectors(Eigen::all, 0);
        float convg =
            ((old_normal_vector - normal_vector).array().abs() / (old_normal_vector.array().abs())).maxCoeff();

        std::cout << "convg: " << convg << std::endl;
        if (convg < conv_thresh)
        {
            std::cout << "Converged at " << iter << " iteration" << std::endl;
            break;
        }
    }
}