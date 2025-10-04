/**
 * @file integrator.cpp
 * @brief Implementation of base integrator functionality
 */

#include "physim/integrators/integrator.hpp"
#include "physim/core/logging.hpp"
#include <cmath>
#include <stdexcept>

namespace physim {
namespace integrators {

VecX Integrator::integrate(double t0, const VecX& y0, double tf, double dt,
                           const DerivativeFunction& f,
                           bool adaptive,
                           double tol) {
    if (tf <= t0) {
        throw std::invalid_argument("Final time must be greater than initial time");
    }

    if (dt <= 0.0) {
        throw std::invalid_argument("Time step must be positive");
    }

    if (adaptive && !supports_adaptive()) {
        PHYSIM_LOG_WARN("{} does not support adaptive stepping, using fixed step", name_);
        adaptive = false;
    }

    // Reset statistics
    stats_.reset();

    double t = t0;
    VecX y = y0;
    double h = dt;

    PHYSIM_LOG_DEBUG("Starting integration: t0={:.6f}, tf={:.6f}, dt={:.6e}, adaptive={}",
                    t0, tf, dt, adaptive);

    while (t < tf) {
        // Adjust step size to hit endpoint exactly
        if (t + h > tf) {
            h = tf - t;
        }

        if (adaptive) {
            // Adaptive step
            bool accepted = adaptive_step(t, y, h, f, tol);

            if (!accepted) {
                stats_.num_rejected_steps++;
                PHYSIM_LOG_TRACE("Step rejected at t={:.6f}, reducing dt to {:.6e}",
                                t, h);
            }

            // Prevent step size from becoming too small
            if (h < 1e-14) {
                PHYSIM_LOG_ERROR("Step size became too small (h={:.6e}) at t={:.6f}",
                               h, t);
                throw std::runtime_error("Integration failed: step size underflow");
            }
        } else {
            // Fixed step
            y = step(t, y, h, f);
            t += h;
            stats_.num_steps++;
        }
    }

    stats_.final_time = t;
    stats_.final_step_size = h;

    PHYSIM_LOG_DEBUG("Integration complete: {} steps, {} function evaluations, {} rejected",
                    stats_.num_steps, stats_.num_function_evals, stats_.num_rejected_steps);

    return y;
}

} // namespace integrators
} // namespace physim
