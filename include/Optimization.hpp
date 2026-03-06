#include "aare/NDArray.hpp"
#include <array>

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
        values_at_grid_points(1, 0) / delta_x; // central difference at y0
    const double Dx_y_minus = values_at_grid_points(0, 2) -
                              values_at_grid_points(0, 0) /
                                  delta_x; // central difference at y-delta_y
    const double Dx_y_plus = values_at_grid_points(2, 2) -
                             values_at_grid_points(2, 0) /
                                 delta_x; // central difference at y+delta_y

    // TODO: does this make sense?
    const double Dx =
        (Dx_y + Dx_y_minus + Dx_y_plus) /
        3; // average of central differences at y0, y-delta_y and y+delta_y

    const double Dy_x =
        values_at_grid_points(2, 1) -
        values_at_grid_points(0, 1) / delta_y; // central difference at x0
    const double Dy_x_minus =
        values_at_grid_points(2, 0) -
        values_at_grid_points(0, 0) / delta_y; // central difference at x-delta
    const double Dy_x_plus =
        values_at_grid_points(2, 2) -
        values_at_grid_points(0, 2) / delta_y; // central difference at x+delta

    const double Dy =
        (Dy_x + Dy_x_minus + Dy_x_plus) /
        3; // average of central differences at x0, x-delta_x and x+delta_x

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
        (values_at_grid_points(2, 2) - values_at_grid_points(2, 1) -
         values_at_grid_points(1, 2) + values_at_grid_points(1, 1)) /
        (delta_x *
         delta_y); // central difference for mixed partial derivative at x0,y0

    // TOD0 Antonio does some averaging !!

    aare::NDArray<double, 2> Hessian({2, 2});
    Hessian(0, 0) = Dxx;
    Hessian(1, 1) = Dyy;
    Hessian(0, 1) = Dxy;
    Hessian(1, 0) = Dxy;

    return Hessian;

    /*
    double Dxx = ((sp[3] - 2 * sp[7] + sp[1]) + (sp[5] - 2 * sp[8] + sp[4]) +
                  (sp[2] - 2 * sp[6] + sp[0])) /
                 (shift_parameter1 * shift_parameter2 *
                  3); // (sp_x+1,y+1 - 2sp_x,y+1 + sp_x-1,y-1)/delta_x*delta_x
                      // etc.
    double Dyy = ((sp[3] - 2 * sp[5] + sp[2]) + (sp[7] - 2 * sp[8] + sp[6]) +
                  (sp[1] - 2 * sp[4] + sp[0])) /
                 (shift_parameter1 * shift_parameter2 * 3);

    double Dyx = ((sp[3] - sp[1]) - sp[2] + sp[0]) /
                 (4 * shift_parameter1 *
                  shift_parameter2); // ((sp_x+1,y+1 - sp_x-1,y+1)/(2*delta_x) -
                                     // (sp_x+1,y-1 -
                                     // sp_x-1,y-1)/(2*delta_x))/2delta_y
    */
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