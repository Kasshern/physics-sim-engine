/**
 * @file integrator.hpp
 * @brief Base class and interfaces for ODE integrators
 *
 * Provides abstract interface for numerical integration of ordinary differential
 * equations (ODEs). Supports both fixed and adaptive time stepping with
 * configurable error tolerances.
 */

#pragma once

#include "physim/core/types.hpp"
#include <functional>
#include <memory>
#include <string>

namespace physim {
namespace integrators {

// ============================================================================
// Integration Statistics
// ============================================================================

/**
 * @brief Statistics collected during integration
 */
struct IntegrationStats {
    size_t num_steps;           ///< Number of integration steps taken
    size_t num_function_evals;  ///< Number of derivative evaluations
    size_t num_rejected_steps;  ///< Number of rejected steps (adaptive only)
    double final_time;          ///< Final simulation time reached [s]
    double final_step_size;     ///< Final step size used [s]
    double max_error;           ///< Maximum error estimate encountered

    IntegrationStats()
        : num_steps(0),
          num_function_evals(0),
          num_rejected_steps(0),
          final_time(0.0),
          final_step_size(0.0),
          max_error(0.0) {}

    /**
     * @brief Reset all statistics to zero
     */
    void reset() {
        num_steps = 0;
        num_function_evals = 0;
        num_rejected_steps = 0;
        final_time = 0.0;
        final_step_size = 0.0;
        max_error = 0.0;
    }
};

// ============================================================================
// Derivative Function Type
// ============================================================================

/**
 * @brief Function signature for computing state derivative
 *
 * @param t Current time [s]
 * @param y Current state vector
 * @return State derivative dy/dt
 *
 * This function computes dy/dt = f(t, y) for the ODE system.
 * For orbital mechanics, this includes velocity and acceleration.
 */
using DerivativeFunction = std::function<VecX(double t, const VecX& y)>;

// ============================================================================
// Integrator Base Class
// ============================================================================

/**
 * @brief Abstract base class for ODE integrators
 *
 * Provides interface for numerical integration of ODEs of the form:
 *   dy/dt = f(t, y)
 *
 * Supports both fixed and adaptive time stepping schemes.
 * Derived classes implement specific methods (RK4, RK45, DOPRI, etc.)
 */
class Integrator {
public:
    /**
     * @brief Constructor
     * @param name Human-readable name of the integration method
     */
    explicit Integrator(std::string name) : name_(std::move(name)) {}

    /**
     * @brief Virtual destructor
     */
    virtual ~Integrator() = default;

    /**
     * @brief Get integrator name
     */
    const std::string& name() const { return name_; }

    /**
     * @brief Get integration statistics
     */
    const IntegrationStats& stats() const { return stats_; }

    /**
     * @brief Reset statistics
     */
    void reset_stats() { stats_.reset(); }

    /**
     * @brief Take a single integration step (fixed step size)
     *
     * @param t Current time [s]
     * @param y Current state vector
     * @param dt Time step [s]
     * @param f Derivative function dy/dt = f(t, y)
     * @return New state vector at time t + dt
     *
     * Advances the state by one time step using fixed step size.
     * Updates statistics (num_steps, num_function_evals).
     */
    virtual VecX step(double t, const VecX& y, double dt,
                      const DerivativeFunction& f) = 0;

    /**
     * @brief Take a single integration step with adaptive step size
     *
     * @param t Current time [s] (will be updated)
     * @param y Current state vector (will be updated)
     * @param dt Suggested time step [s] (will be updated)
     * @param f Derivative function dy/dt = f(t, y)
     * @param tol Error tolerance (relative)
     * @return True if step was accepted, false if rejected
     *
     * Advances the state with automatic step size adjustment.
     * Updates t, y, and dt based on error estimate.
     * May reject step and retry with smaller dt.
     *
     * For fixed-step methods, this delegates to step() with dt unchanged.
     */
    virtual bool adaptive_step(double& t, VecX& y, double& dt,
                              const DerivativeFunction& f,
                              double tol = 1e-12) {
        // Default: fixed step (no adaptation)
        y = step(t, y, dt, f);
        t += dt;
        return true;
    }

    /**
     * @brief Integrate from t0 to tf
     *
     * @param t0 Initial time [s]
     * @param y0 Initial state vector
     * @param tf Final time [s]
     * @param dt Time step [s] (for fixed) or initial step (for adaptive)
     * @param f Derivative function
     * @param adaptive Use adaptive stepping if available
     * @param tol Error tolerance for adaptive methods
     * @return Final state vector at time tf
     *
     * Integrates the ODE from t0 to tf using repeated steps.
     * Automatically handles endpoint to exactly hit tf.
     */
    VecX integrate(double t0, const VecX& y0, double tf, double dt,
                   const DerivativeFunction& f,
                   bool adaptive = false,
                   double tol = 1e-12);

    /**
     * @brief Get order of accuracy
     * @return Order (e.g., 4 for RK4, 5 for RK45)
     */
    virtual int order() const = 0;

    /**
     * @brief Check if method supports adaptive stepping
     */
    virtual bool supports_adaptive() const { return false; }

protected:
    std::string name_;           ///< Integrator name
    IntegrationStats stats_;     ///< Integration statistics
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert State to vector for integration
 * @param state State structure
 * @return 6D vector [x, y, z, vx, vy, vz]
 */
inline VecX state_to_vector_6dof(const State& state) {
    VecX y(6);
    y.segment<3>(0) = state.position;
    y.segment<3>(3) = state.velocity;
    return y;
}

/**
 * @brief Convert vector back to State (position and velocity only)
 * @param y 6D vector
 * @param state State to update (only position and velocity modified)
 */
inline void vector_to_state_6dof(const VecX& y, State& state) {
    state.position = y.segment<3>(0);
    state.velocity = y.segment<3>(3);
}

/**
 * @brief Create derivative function for simple Newtonian dynamics
 * @param acceleration_func Function computing acceleration(t, pos, vel)
 * @return Derivative function dy/dt = [v, a(t, r, v)]
 *
 * For equations of motion: r'' = a(t, r, r')
 * Converts to first-order system: [r', v']^T = [v, a]^T
 */
inline DerivativeFunction newtonian_derivative(
    std::function<Vec3(double t, const Vec3& pos, const Vec3& vel)> acceleration_func) {
    return [acceleration_func](double t, const VecX& y) -> VecX {
        Vec3 pos = y.segment<3>(0);
        Vec3 vel = y.segment<3>(3);
        Vec3 acc = acceleration_func(t, pos, vel);

        VecX dydt(6);
        dydt.segment<3>(0) = vel;  // dr/dt = v
        dydt.segment<3>(3) = acc;  // dv/dt = a
        return dydt;
    };
}

} // namespace integrators
} // namespace physim
