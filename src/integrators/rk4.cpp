/**
 * @file rk4.cpp
 * @brief Implementation of RK4 integrator
 */

#include "physim/integrators/rk4.hpp"

namespace physim {
namespace integrators {

VecX RK4::step(double t, const VecX& y, double dt, const DerivativeFunction& f) {
    // Classic 4th order Runge-Kutta method
    //
    // Evaluates derivative at 4 points:
    //   1. Start of interval: k1 = f(t, y)
    //   2. Midpoint using k1: k2 = f(t + h/2, y + h*k1/2)
    //   3. Midpoint using k2: k3 = f(t + h/2, y + h*k2/2)
    //   4. End of interval: k4 = f(t + h, y + h*k3)
    //
    // Then combines with weights: y_new = y + h/6 * (k1 + 2*k2 + 2*k3 + k4)
    //
    // This provides 4th order accuracy: local error ~ O(h^5), global ~ O(h^4)

    const VecX k1 = f(t, y);
    const VecX k2 = f(t + 0.5 * dt, y + 0.5 * dt * k1);
    const VecX k3 = f(t + 0.5 * dt, y + 0.5 * dt * k2);
    const VecX k4 = f(t + dt, y + dt * k3);

    // Update function evaluation counter
    stats_.num_function_evals += 4;

    // Weighted combination (Simpson's rule)
    return y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4);
}

} // namespace integrators
} // namespace physim
