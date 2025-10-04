/**
 * @file rk45.cpp
 * @brief Implementation of RK45 adaptive integrator
 */

#include "physim/integrators/rk45.hpp"
#include "physim/core/logging.hpp"
#include <algorithm>
#include <cmath>

namespace physim {
namespace integrators {

RK45::RK45()
    : Integrator("RK45"),
      safety_factor_(0.9),
      min_step_factor_(0.2),
      max_step_factor_(5.0) {}

VecX RK45::step(double t, const VecX& y, double dt, const DerivativeFunction& f) {
    // Compute the 6 stages of RK45
    const VecX k1 = f(t, y);
    const VecX k2 = f(t + C2 * dt, y + dt * (A21 * k1));
    const VecX k3 = f(t + C3 * dt, y + dt * (A31 * k1 + A32 * k2));
    const VecX k4 = f(t + C4 * dt, y + dt * (A41 * k1 + A42 * k2 + A43 * k3));
    const VecX k5 = f(t + C5 * dt, y + dt * (A51 * k1 + A52 * k2 + A53 * k3 + A54 * k4));
    const VecX k6 = f(t + C6 * dt, y + dt * (A61 * k1 + A62 * k2 + A63 * k3 + A64 * k4 + A65 * k5));

    stats_.num_function_evals += 6;

    // Use 5th order formula for the step (more accurate)
    return y + dt * (B51 * k1 + B53 * k3 + B54 * k4 + B55 * k5 + B56 * k6);
}

bool RK45::adaptive_step(double& t, VecX& y, double& dt,
                        const DerivativeFunction& f,
                        double tol) {
    // Compute the 6 stages
    const VecX k1 = f(t, y);
    const VecX k2 = f(t + C2 * dt, y + dt * (A21 * k1));
    const VecX k3 = f(t + C3 * dt, y + dt * (A31 * k1 + A32 * k2));
    const VecX k4 = f(t + C4 * dt, y + dt * (A41 * k1 + A42 * k2 + A43 * k3));
    const VecX k5 = f(t + C5 * dt, y + dt * (A51 * k1 + A52 * k2 + A53 * k3 + A54 * k4));
    const VecX k6 = f(t + C6 * dt, y + dt * (A61 * k1 + A62 * k2 + A63 * k3 + A64 * k4 + A65 * k5));

    stats_.num_function_evals += 6;

    // Compute 4th and 5th order solutions
    const VecX y4 = y + dt * (B41 * k1 + B43 * k3 + B44 * k4 + B45 * k5);
    const VecX y5 = y + dt * (B51 * k1 + B53 * k3 + B54 * k4 + B55 * k5 + B56 * k6);

    // Error estimate: difference between 4th and 5th order
    const VecX error = y5 - y4;

    // Compute error norm (scaled by tolerance and state magnitude)
    // Use mixed absolute-relative error: err / (atol + rtol * |y|)
    const double atol = tol;  // Absolute tolerance
    const double rtol = tol;  // Relative tolerance

    double err_norm = 0.0;
    for (int i = 0; i < y.size(); ++i) {
        const double scale = atol + rtol * std::abs(y(i));
        const double scaled_err = std::abs(error(i)) / scale;
        err_norm = std::max(err_norm, scaled_err);
    }

    // Update max error stat
    stats_.max_error = std::max(stats_.max_error, err_norm);

    // Compute optimal step size using standard formula
    // h_new = h * safety * (tol / error)^(1/order)
    // For RK45, order = 5 for step size control
    const double exponent = 1.0 / 5.0;
    double step_factor = safety_factor_ * std::pow(1.0 / err_norm, exponent);

    // Clamp step factor to prevent wild oscillations
    step_factor = std::clamp(step_factor, min_step_factor_, max_step_factor_);

    // Decide whether to accept or reject step
    if (err_norm <= 1.0) {
        // Accept step
        y = y5;  // Use 5th order solution
        t += dt;
        stats_.num_steps++;

        // Update step size for next step (can grow)
        dt *= step_factor;

        PHYSIM_LOG_TRACE("Step accepted: t={:.6f}, error={:.3e}, new_dt={:.3e}",
                        t, err_norm, dt);

        return true;
    } else {
        // Reject step
        stats_.num_rejected_steps++;

        // Reduce step size (must shrink)
        dt *= step_factor;

        PHYSIM_LOG_TRACE("Step rejected: t={:.6f}, error={:.3e}, new_dt={:.3e}",
                        t, err_norm, dt);

        return false;
    }
}

void RK45::set_safety_factor(double safety) {
    if (safety <= 0.0 || safety > 1.0) {
        PHYSIM_LOG_WARN("Safety factor {} out of range [0, 1], clamping", safety);
        safety = std::clamp(safety, 0.1, 1.0);
    }
    safety_factor_ = safety;
    PHYSIM_LOG_DEBUG("RK45 safety factor set to {}", safety_factor_);
}

void RK45::set_step_bounds(double min_factor, double max_factor) {
    if (min_factor <= 0.0 || min_factor >= 1.0) {
        PHYSIM_LOG_WARN("Min step factor {} invalid, using 0.2", min_factor);
        min_factor = 0.2;
    }
    if (max_factor <= 1.0 || max_factor > 10.0) {
        PHYSIM_LOG_WARN("Max step factor {} invalid, using 5.0", max_factor);
        max_factor = 5.0;
    }
    min_step_factor_ = min_factor;
    max_step_factor_ = max_factor;
    PHYSIM_LOG_DEBUG("RK45 step bounds: [{}, {}]", min_step_factor_, max_step_factor_);
}

} // namespace integrators
} // namespace physim
