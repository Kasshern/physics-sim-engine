/**
 * @file rk45.hpp
 * @brief Runge-Kutta-Fehlberg adaptive integrator
 *
 * Implements RK45 method with embedded error estimation for automatic
 * step size control. Provides 4th/5th order accuracy with adaptive
 * time stepping based on local error estimates.
 */

#pragma once

#include "physim/integrators/integrator.hpp"

namespace physim {
namespace integrators {

/**
 * @brief Runge-Kutta-Fehlberg adaptive integrator (RK45)
 *
 * Implements the Runge-Kutta-Fehlberg method with embedded 4th and 5th
 * order formulas for error estimation and automatic step size control.
 *
 * Algorithm:
 * - Uses 6 function evaluations per step
 * - Computes both 4th and 5th order estimates
 * - Error = |y5 - y4|
 * - Adjusts step size based on error vs tolerance
 *
 * Butcher tableau coefficients from Fehlberg (1969):
 *   0   |
 *   1/4 | 1/4
 *   3/8 | 3/32      9/32
 *   12/13 | 1932/2197  -7200/2197  7296/2197
 *   1   | 439/216    -8         3680/513    -845/4104
 *   1/2 | -8/27      2          -3544/2565  1859/4104  -11/40
 *
 * Properties:
 * - Order: 4th order with 5th order error estimate
 * - Function evaluations: 6 per accepted step
 * - Adaptive step size: Yes (automatic error control)
 * - FSAL: No (not First Same As Last, uses 6 evals)
 *
 * Step size control:
 * - Accepts step if error < tolerance
 * - Increases step if error << tolerance
 * - Decreases step if error > tolerance
 * - Safety factors prevent oscillation
 *
 * Use cases:
 * - General purpose with automatic error control
 * - Systems with varying dynamics (fast/slow regions)
 * - When error tolerance is more important than efficiency
 *
 * Advantages over fixed RK4:
 * - Automatic step size adaptation
 * - Guaranteed error tolerance (if conservative)
 * - Efficient in smooth regions (large steps)
 * - Robust in chaotic regions (small steps)
 */
class RK45 : public Integrator {
public:
    /**
     * @brief Constructor with default parameters
     */
    RK45();

    /**
     * @brief Take a single RK45 step (fixed step)
     *
     * @param t Current time [s]
     * @param y Current state vector
     * @param dt Time step [s]
     * @param f Derivative function
     * @return New state vector at t + dt
     *
     * Uses 5th order formula (ignores error estimate).
     * For adaptive stepping, use adaptive_step() instead.
     */
    VecX step(double t, const VecX& y, double dt,
              const DerivativeFunction& f) override;

    /**
     * @brief Take adaptive step with error control
     *
     * @param t Current time [s] (updated on success)
     * @param y Current state (updated on success)
     * @param dt Time step [s] (updated based on error)
     * @param f Derivative function
     * @param tol Error tolerance (relative)
     * @return true if step accepted, false if rejected
     *
     * Adjusts dt automatically based on local error estimate.
     * May reject step and suggest smaller dt.
     */
    bool adaptive_step(double& t, VecX& y, double& dt,
                      const DerivativeFunction& f,
                      double tol = 1e-12) override;

    /**
     * @brief Get order of accuracy
     * @return 5 (5th order method with 4th order error estimate)
     */
    int order() const override { return 5; }

    /**
     * @brief Check if adaptive stepping is supported
     * @return true (RK45 supports adaptive stepping)
     */
    bool supports_adaptive() const override { return true; }

    /**
     * @brief Set safety factor for step size adjustment
     * @param safety Safety factor (default 0.9, range [0.1, 1.0])
     *
     * Controls conservativeness of step size increase:
     * - Higher values: more aggressive (faster but less safe)
     * - Lower values: more conservative (slower but safer)
     */
    void set_safety_factor(double safety);

    /**
     * @brief Set minimum/maximum step size growth factors
     * @param min_factor Minimum shrink factor (default 0.2)
     * @param max_factor Maximum growth factor (default 5.0)
     */
    void set_step_bounds(double min_factor, double max_factor);

private:
    // Butcher tableau coefficients (Fehlberg 1969)
    static constexpr double A21 = 1.0 / 4.0;
    static constexpr double A31 = 3.0 / 32.0;
    static constexpr double A32 = 9.0 / 32.0;
    static constexpr double A41 = 1932.0 / 2197.0;
    static constexpr double A42 = -7200.0 / 2197.0;
    static constexpr double A43 = 7296.0 / 2197.0;
    static constexpr double A51 = 439.0 / 216.0;
    static constexpr double A52 = -8.0;
    static constexpr double A53 = 3680.0 / 513.0;
    static constexpr double A54 = -845.0 / 4104.0;
    static constexpr double A61 = -8.0 / 27.0;
    static constexpr double A62 = 2.0;
    static constexpr double A63 = -3544.0 / 2565.0;
    static constexpr double A64 = 1859.0 / 4104.0;
    static constexpr double A65 = -11.0 / 40.0;

    // Nodes
    static constexpr double C2 = 1.0 / 4.0;
    static constexpr double C3 = 3.0 / 8.0;
    static constexpr double C4 = 12.0 / 13.0;
    static constexpr double C5 = 1.0;
    static constexpr double C6 = 0.5;

    // Weights for 4th order solution
    static constexpr double B41 = 25.0 / 216.0;
    static constexpr double B43 = 1408.0 / 2565.0;
    static constexpr double B44 = 2197.0 / 4104.0;
    static constexpr double B45 = -1.0 / 5.0;

    // Weights for 5th order solution
    static constexpr double B51 = 16.0 / 135.0;
    static constexpr double B53 = 6656.0 / 12825.0;
    static constexpr double B54 = 28561.0 / 56430.0;
    static constexpr double B55 = -9.0 / 50.0;
    static constexpr double B56 = 2.0 / 55.0;

    double safety_factor_;     ///< Safety factor for step size adjustment (0.9)
    double min_step_factor_;   ///< Minimum step shrink factor (0.2)
    double max_step_factor_;   ///< Maximum step growth factor (5.0)
};

} // namespace integrators
} // namespace physim
