#ifndef IRLS_PLANE_FITTING_HPP_
#define IRLS_PLANE_FITTING_HPP_

#include <chrono>
#include <eigen3/Eigen/Dense>
#include <iostream>
#include <memory>
#include <utility>
#include <vector>

class PlaneFittingIRLS
{
  private:
    Eigen::MatrixX3f points_;
    size_t number_of_points_;

  public:
    // Default Constructor
    PlaneFittingIRLS(const Eigen::MatrixX3f &points) : points_(points)
    {
        number_of_points_ = static_cast<size_t>(points.rows());
    };

    // Copy Constructor
    PlaneFittingIRLS(const PlaneFittingIRLS &rhs) : points_(rhs.points_), number_of_points_(rhs.number_of_points_){};

    // Move Constructor
    PlaneFittingIRLS(PlaneFittingIRLS &&) noexcept;

    // Destructor
    virtual ~PlaneFittingIRLS() = default;

    // Copy operator
    PlaneFittingIRLS &operator=(const PlaneFittingIRLS &);

    // Move operator
    PlaneFittingIRLS &operator=(PlaneFittingIRLS &&) noexcept;

    inline size_t getNumberOfPoints() const noexcept
    {
        return number_of_points_;
    }

    inline Eigen::MatrixX3f getPoints() const noexcept
    {
        return points_;
    }

    void fitPlane(size_t max_iterations = 1);
};

#endif /* IRLS_PLANE_FITTING_HPP_ */