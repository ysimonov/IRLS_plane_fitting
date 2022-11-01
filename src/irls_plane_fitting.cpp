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

    auto start_time = std::chrono::high_resolution_clock::now();

    // get points
    Eigen::MatrixX3f points_xyz = this->getPoints(); // N x 3
    size_t number_of_points = points_xyz.rows();     // N

    // compute centroid vector
    Eigen::RowVector3f centroid = points_xyz.colwise().mean(); // 1 x 3

    // compute covariance matrix
    Eigen::MatrixX3f points_xyz_centered = points_xyz.rowwise() - centroid; // N x 3
    Eigen::Matrix3f covariance_matrix =
        (points_xyz_centered.transpose() * points_xyz_centered) / static_cast<double>(number_of_points - 1); // 3 x 3

    // eigendecomposition of the covariance matrix - returns ordered eigenvalues
    Eigen::SelfAdjointEigenSolver<Eigen::Matrix3f> eigensolver;
    eigensolver.compute(covariance_matrix);
    Eigen::Vector3f eigenvalues = eigensolver.eigenvalues();   // 3 x 1
    Eigen::Matrix3f eigenvectors = eigensolver.eigenvectors(); // 3 x 3

    // the principal component (first eigenvector) corresponds to a rotation matrix (direction of the normal vector)
    Eigen::Vector3f normal_curr = eigenvectors(Eigen::all, 0); // 3 x 1

    // normalize normal vector
    normal_curr.normalize();                   // 3 x 1
    Eigen::Vector3f normal_prev = normal_curr; // 3 x 1

    // get additional centered values
    Eigen::RowVector3f xk_curr, xk_prev;
    xk_prev = points_xyz_centered;
    xk_curr = xk_prev;

    // iterative reweighted least squares
    for (size_t iter = 0; iter < max_iterations; ++iter)
    {
        // update old normal vector
        normal_prev = normal_curr; // 3 x 1

        // calculate distances to the plane (dot product between points and normal vector)
        Eigen::VectorXf distances = points_xyz_centered * normal_curr; // N x 1

        // calculate Welsh function
        Eigen::VectorXf weights = Eigen::exp(-(distances.array() / k_welsh).pow(2.0f)); // N x 1

        // calculate xk_curr
        xk_curr = weights.transpose() * (points_xyz_centered.rowwise() - xk_prev) / weights.sum();
        xk_prev = xk_curr;

        // calculate covariance matrix
        Eigen::MatrixX3f temp = points_xyz_centered.rowwise() - xk_curr; // N x 3
        Eigen::Matrix3Xf temp_transposed = temp.transpose();             // 3 x N

        covariance_matrix = Eigen::Matrix3f::Zero(); // 3 x 3
        for (size_t i = 0; i < number_of_points; ++i)
        {
            // removed weight
            covariance_matrix += temp_transposed(Eigen::all, i) * temp(i, Eigen::all);
        }

        // recalculate eigenvalues and eigenvectors
        eigensolver.compute(covariance_matrix);
        eigenvalues = eigensolver.eigenvalues();   // 3 x 1
        eigenvectors = eigensolver.eigenvectors(); // 3 x 3

        // check convergence criteria
        normal_curr = eigenvectors(Eigen::all, 0); // 3 x 1

        // normalize normal vector
        normal_curr.normalize(); // 3 x 1

        float convg = ((normal_prev - normal_curr).array().abs() / (normal_prev.array().abs())).maxCoeff();

        std::cout << iter << ") Normal Vector: " << normal_curr.transpose() << " | Convg: " << convg << std::endl;

        if (convg < conv_thresh)
        {
            std::cout << "Converged at " << iter << " iteration" << std::endl;
            break;
        }
    }

    auto stop_time = std::chrono::high_resolution_clock::now();
    auto elapsed_time = std::chrono::duration_cast<std::chrono::milliseconds>(stop_time - start_time).count() / 1000.0;
    std::cout << "Elapsed time: " << elapsed_time << " seconds" << std::endl;
    std::cout << "Normal Vector: " << normal_curr.transpose() << std::endl;
}