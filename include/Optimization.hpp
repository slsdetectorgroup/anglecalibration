#include "aare/NDArray.hpp"
#include "logger.hpp"
#include <array>

namespace angcal {

/**
 * @brief calculates the gradient based on central finite differences
 * @param values_at_grid_points vector of function values at the grid points
 * values_at_grid_points[0,0] = f_x-delta_x,y-delta_y,
 *values_at_grid_points[0,1] = f_x,y-delta_y etc.
 * @param delta_x step size in x direction
 * @param delta_y step size in y direction
 **/
std::array<double, 2>
calculate_gradient(const aare::NDArray<double, 2> &values_at_grid_points,
                   const double delta_x, const double delta_y) {

    const double Dx_y =
        values_at_grid_points(1, 2) -
        values_at_grid_points(1, 0) / (2 * delta_x); // central difference at y0
    /*
    const double Dx_y_minus =
        values_at_grid_points(0, 2) -
        values_at_grid_points(0, 0) /
            (2 * delta_x); // central difference at y-delta_y
    const double Dx_y_plus =
        values_at_grid_points(2, 2) -
        values_at_grid_points(2, 0) /
            (2 * delta_x); // central difference at y+delta_y

    // TODO: does this make sense?
    const double Dx =
        (Dx_y + Dx_y_minus + Dx_y_plus) /
        3; // average of central differences at y0, y-delta_y and y+delta_y
    */

    const double Dy_x =
        values_at_grid_points(2, 1) -
        values_at_grid_points(0, 1) / (2 * delta_y); // central difference at x0
    /*
    const double Dy_x_minus =
        values_at_grid_points(2, 0) -
        values_at_grid_points(0, 0) /
            (2 * delta_y); // central difference at x-delta
    const double Dy_x_plus = values_at_grid_points(2, 2) -
                             values_at_grid_points(0, 2) /
                                 (2 * delta_y); // central difference at x+delta

    const double Dy =
        (Dy_x + Dy_x_minus + Dy_x_plus) /
        3; // average of central differences at x0, x-delta_x and x+delta_x
    */

    return {Dx_y, Dy_x}; // TODO: shouldnt I calculate it at a specific point?
}

/**
 * @brief calculates the Hessian matrix based on finite differences
 * @param values_at_grid_points vector of function values at the grid points
 * values_at_grid_points[0,0] = f_x-delta_x,y-delta_y,
 *values_at_grid_points[0,1] = f_x,y-delta_y etc.
 **/
aare::NDArray<double, 2>
calculate_Hessian(const aare::NDArray<double, 2> &values_at_grid_points,
                  const double delta_x, const double delta_y) {

    const double Dxx =
        values_at_grid_points(1, 2) - 2 * values_at_grid_points(1, 1) +
        values_at_grid_points(1, 0) /
            (delta_x * delta_x); // second order central difference at y0
    const double Dyy =
        values_at_grid_points(2, 1) - 2 * values_at_grid_points(1, 1) +
        values_at_grid_points(0, 1) /
            (delta_y * delta_y); // second order central difference at x0
    const double Dxy =
        (values_at_grid_points(2, 2) - values_at_grid_points(2, 0) -
         values_at_grid_points(0, 2) + values_at_grid_points(0, 0)) /
        (4 * delta_x *
         delta_y); // central difference for mixed partial derivative at x0,y0

    // TOD0 Antonio does some averaging !!

    aare::NDArray<double, 2> Hessian({2, 2});
    Hessian(0, 0) = Dxx;
    Hessian(1, 1) = Dyy;
    Hessian(0, 1) = Dxy;
    Hessian(1, 0) = Dxy;

    return Hessian;
}

/**
 * @brief calculates the regularized inverse of the Hessian matrix
 * @param Hessian Hessian matrix to be inverted
 * @return regularized inverse of the Hessian matrix
 */
aare::NDArray<double, 2> calculate_regularized_inverse_Hessian(
    const aare::NDView<double, 2> Hessian,
    const double regularization_parameter = 0.0) {

    // hessian + identity*regularization_parameter
    auto regularized_Hessian = Hessian;
    regularized_Hessian(0, 0) += regularization_parameter;
    regularized_Hessian(1, 1) += regularization_parameter;

    const double determinant =
        regularized_Hessian(0, 0) * regularized_Hessian(1, 1) -
        regularized_Hessian(0, 1) * regularized_Hessian(1, 0);
    aare::NDArray<double, 2> inverse_Hessian({2, 2});
    inverse_Hessian(0, 0) = regularized_Hessian(1, 1) / determinant;
    inverse_Hessian(1, 1) = regularized_Hessian(0, 0) / determinant;
    inverse_Hessian(0, 1) = -regularized_Hessian(0, 1) / determinant;
    inverse_Hessian(1, 0) = -regularized_Hessian(1, 0) / determinant;

    return inverse_Hessian;
}

std::array<double, 2>
calculate_eigenvalues(const aare::NDView<double, 2> matrix) {

    double determinant =
        matrix(0, 0) * matrix(1, 1) - matrix(0, 1) * matrix(1, 0);
    double trace = matrix(0, 0) + matrix(1, 1);
    const double eigenvalue1 =
        0.5 * trace - 0.5 * std::sqrt(trace * trace - 4 * determinant);
    const double eigenvalue2 =
        0.5 * trace + 0.5 * std::sqrt(trace * trace - 4 * determinant);

    return {eigenvalue1, eigenvalue2};
}

template <size_t Size, typename Func>
std::pair<std::array<double, Size>, double>
line_search(const std::array<double, Size> &initial_value,
            const double f_initial,
            const std::array<double, Size> &direction_of_steepest_descent,
            Func &&function) {

    double alpha = 1.0; // step size for line search

    std::array<double, Size> current_value = initial_value;

    double f_prev = f_initial;
    double f_next{};

    std::array<double, Size> next_value;
    std::transform(current_value.begin(), current_value.end(),
                   direction_of_steepest_descent.begin(), next_value.begin(),
                   [&alpha](double x, double d) { return x - alpha * d; });

    bool found_best_parameters = false;

    while (!found_best_parameters) {

        f_next = function(next_value);

        found_best_parameters = f_next < f_prev;

        if (f_next >= f_prev) {
            LOG(TLogLevel::logDEBUG1)
                << fmt::format("function value increased to {} ", f_next);

            // refine step size
            alpha *= 0.5;
            std::transform(
                current_value.begin(), current_value.end(),
                direction_of_steepest_descent.begin(), next_value.begin(),
                [&alpha](double x, double d) { return x - alpha * d; });
        } else {
            LOG(TLogLevel::logDEBUG1) << fmt::format(
                "found better parameters with function value: {}", f_next);
        }
    }

    return {next_value, f_next};
}

std::tuple<std::array<double, 2>, double, double> newton_method(
    const double delta_x, const double delta_y,
    const std::array<double, 2> &initial_parameters, const double f_initial,
    const std::function<double(const std::array<double, 2> &)> &function) {

    constexpr double tolerance1 =
        0.001; // tolerance for eigenvalues of Hessian, TODO: dont know if this
               // should be configurable

    NDArray<double, 2> function_at_grid_points({3, 3});
    // calculate function values at grid points
    function_at_grid_points(0, 0) = function(
        {initial_parameters[0] - delta_x, initial_parameters[1] + delta_y});
    function_at_grid_points(0, 1) =
        function({initial_parameters[0], initial_parameters[1] + delta_y});
    function_at_grid_points(0, 2) = function(
        {initial_parameters[0] + delta_x, initial_parameters[1] + delta_y});
    function_at_grid_points(1, 0) =
        function({initial_parameters[0] - delta_x, initial_parameters[1]});
    function_at_grid_points(1, 1) = f_initial;
    function_at_grid_points(1, 2) =
        function({initial_parameters[0] + delta_x, initial_parameters[1]});
    function_at_grid_points(2, 0) = function(
        {initial_parameters[0] - delta_x, initial_parameters[1] - delta_y});
    function_at_grid_points(2, 1) =
        function({initial_parameters[0], initial_parameters[1] - delta_y});
    function_at_grid_points(2, 2) = function(
        {initial_parameters[0] + delta_x, initial_parameters[1] - delta_y});

    // calculate gradient
    auto gradient =
        calculate_gradient(function_at_grid_points, delta_x, delta_y);

    // calculate Hessian
    auto Hessian = calculate_Hessian(function_at_grid_points, delta_x, delta_y);

    auto eigenvalues = calculate_eigenvalues(Hessian.view());

    const double regularization_term = std::max(
        0.0,
        tolerance1 -
            std::min(eigenvalues[0],
                     eigenvalues[1])); // eigenvalues need to be bigger than
                                       // tolerance ensures that eigenvalues
                                       // of regularized Hessian are
                                       // positive - thus ensuring that we
                                       // have a descent direction

    auto inverse_regularized_Hessian = calculate_regularized_inverse_Hessian(
        Hessian.view(), regularization_term);

    // Newton step x_k+1 = x_k - Hessian^(-1)*gradient
    // steepest descent in direction of Hessian^(-1)*gradient
    std::array<double, 2> steepest_descent{
        inverse_regularized_Hessian(0, 0) * gradient[0] +
            inverse_regularized_Hessian(0, 1) * gradient[1],
        inverse_regularized_Hessian(1, 0) * gradient[0] +
            inverse_regularized_Hessian(1, 1) * gradient[1]};

    LOG(TLogLevel::logDEBUG1)
        << fmt::format("direction of steepest descent: [{},{}]",
                       steepest_descent[0], steepest_descent[1]);

    auto [next_value, f_next] =
        line_search(initial_parameters, f_initial, steepest_descent, function);

    std::array<double, 2> change_in_parameter = {
        next_value[0] - initial_parameters[0],
        next_value[1] - initial_parameters[1]};

    double norm_change_in_parameter =
        std::sqrt(change_in_parameter[0] * change_in_parameter[0] +
                  change_in_parameter[1] * change_in_parameter[1]);

    double optimality_criterion =
        std::abs(change_in_parameter[0] * steepest_descent[0] +
                 change_in_parameter[1] * steepest_descent[1]) /
        std::max(
            std::numeric_limits<double>::epsilon(),
            norm_change_in_parameter); // |gradient dot change in parameter| /
                                       // norm of change in parameter

    return {next_value, f_next, optimality_criterion};
}

} // namespace angcal