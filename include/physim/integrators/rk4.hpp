/**
 * @file rk4.hpp
 * @brief Classic 4th order Runge-Kutta integrator
 *
 * Implements the standard RK4 method, one of the most widely used
 * fixed-step ODE solvers. Provides good accuracy with moderate
 * computational cost (4 function evaluations per step).
 */

#pragma once

#include "physim/integrators/integrator.hpp"

namespace physim {
namespace integrators {

/**
 * @brief Classic 4th order Runge-Kutta integrator
 *
 * Implements the standard RK4 method:
 *   k1 = f(t, y)
 *   k2 = f(t + h/2, y + h*k1/2)
 *   k3 = f(t + h/2, y + h*k2/2)
 *   k4 = f(t + h, y + h*k3)
 *   y_{n+1} = y_n + h/6 * (k1 + 2*k2 + 2*k3 + k4)
 *
 * Properties:
 * - Order of accuracy: O(h^4) local, O(h^3) global
 * - Stability: Conditionally stable (CFL condition applies)
 * - Function evaluations: 4 per step
 * - Error: ~h^5 local truncation error
 *
 * Use cases:
 * - General purpose ODE integration with fixed step
 * - Smooth systems where adaptive stepping is not required
 * - Benchmarking and validation
 *
 * NOT recommended for:
 * - Stiff systems (use implicit methods)
 * - Systems requiring tight error control (use RK45/DOPRI)
 * - Long-term energy conservation (use symplectic methods)
 */
class RK4 : public Integrator {
public:
    /**
     * @brief Constructor
     */
    RK4() : Integrator("RK4") {}

    /**
     * @brief Take a single RK4 step
     *
     * @param t Current time [s]
     * @param y Current state vector
     * @param dt Time step [s]
     * @param f Derivative function dy/dt = f(t, y)
     * @return New state vector at time t + dt
     *
     * Performs one step of the classic 4th order Runge-Kutta method.
     * Requires 4 function evaluations.
     */
    VecX step(double t, const VecX& y, double dt,
              const DerivativeFunction& f) override;

    /**
     * @brief Get order of accuracy
     * @return 4 (4th order method)
     */
    int order() const override { return 4; }

    /**
     * @brief Check if adaptive stepping is supported
     * @return false (RK4 is fixed-step only)
     */
    bool supports_adaptive() const override { return false; }
};

} // namespace integrators
} // namespace physim
