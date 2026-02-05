#pragma once
#include <cmath>

/*
 * returns the propagated error e.g. combined variance of two independent
 * variables using standard error propagation formula z = f(x,y) give the
 * variance of z var(z) = (df/dx)^2 * var(x) + (df/dy)^2 * var(y)
 * @param derivative_1: derivative of z by first variable
 * @param derivative_2: derivative of z by second variable
 * @param variance_1: variance of first variable
 * @param variance_2: variance of second variable
 * @return propagated variance
 */
inline double error_propagation(const double derivative_1,
                                const double derivative_2,
                                const double variance_1,
                                const double variance_2) {
    return std::pow(derivative_1, 2) * variance_1 +
           std::pow(derivative_2, 2) * variance_2;
}
