/**
 * @file earth_moon_orbit.cpp
 * @brief Earth-Moon two-body problem example
 *
 * Demonstrates:
 * - RK4 and RK45 integrators
 * - Point mass gravity force
 * - Energy conservation tracking
 * - Orbital element calculation
 */

#include "physim/core/types.hpp"
#include "physim/core/constants.hpp"
#include "physim/core/logging.hpp"
#include "physim/integrators/integrator.hpp"
#include "physim/integrators/rk4.hpp"
#include "physim/integrators/rk45.hpp"
#include "physim/forces/force.hpp"
#include "physim/forces/gravity.hpp"

#include <iostream>
#include <iomanip>
#include <vector>
#include <cmath>

using namespace physim;

// Helper: Compute orbital elements from state
struct OrbitalElements {
    double a;       // Semi-major axis [m]
    double e;       // Eccentricity
    double i;       // Inclination [rad]
    double omega;   // Argument of periapsis [rad]
    double Omega;   // Right ascension of ascending node [rad]
    double nu;      // True anomaly [rad]
};

OrbitalElements compute_elements(const Vec3& r, const Vec3& v, double mu) {
    OrbitalElements elem;

    const double r_mag = r.norm();
    const double v_mag = v.norm();

    // Specific orbital energy
    const double epsilon = v_mag * v_mag / 2.0 - mu / r_mag;

    // Semi-major axis
    elem.a = -mu / (2.0 * epsilon);

    // Angular momentum vector
    const Vec3 h = r.cross(v);
    const double h_mag = h.norm();

    // Eccentricity vector
    const Vec3 e_vec = (v.cross(h) / mu) - r / r_mag;
    elem.e = e_vec.norm();

    // Inclination
    elem.i = std::acos(h(2) / h_mag);

    // Node vector
    const Vec3 n(-h(1), h(0), 0.0);
    const double n_mag = n.norm();

    // RAAN
    if (n_mag > 1e-10) {
        elem.Omega = std::acos(n(0) / n_mag);
        if (n(1) < 0.0) elem.Omega = constants::TWO_PI - elem.Omega;
    } else {
        elem.Omega = 0.0;
    }

    // Argument of periapsis
    if (n_mag > 1e-10 && elem.e > 1e-10) {
        elem.omega = std::acos(n.dot(e_vec) / (n_mag * elem.e));
        if (e_vec(2) < 0.0) elem.omega = constants::TWO_PI - elem.omega;
    } else {
        elem.omega = 0.0;
    }

    // True anomaly
    if (elem.e > 1e-10) {
        elem.nu = std::acos(e_vec.dot(r) / (elem.e * r_mag));
        if (r.dot(v) < 0.0) elem.nu = constants::TWO_PI - elem.nu;
    } else {
        elem.nu = 0.0;
    }

    return elem;
}

int main() {
    // Initialize logging
    logging::init(spdlog::level::info, false);
    PHYSIM_LOG_INFO("=== Earth-Moon Two-Body Problem ===");

    // Earth-Moon system parameters
    const double mu_earth = constants::gm::EARTH;
    const double mu_moon = constants::gm::MOON;

    // Initial conditions: Moon at perigee of its orbit
    // Semi-major axis: ~384,400 km
    // Eccentricity: ~0.0549
    const double a_moon = 384400e3;                    // m
    const double e_moon = 0.0549;
    const double r_perigee = a_moon * (1.0 - e_moon);  // Perigee distance

    Vec3 r0(r_perigee, 0, 0);                         // Position at perigee
    const double v_perigee = std::sqrt(mu_earth * (1.0 + e_moon) / r_perigee);
    Vec3 v0(0, v_perigee, 0);                         // Velocity at perigee

    PHYSIM_LOG_INFO("Initial conditions:");
    PHYSIM_LOG_INFO("  Position: [{:.3e}, {:.3e}, {:.3e}] m", r0(0), r0(1), r0(2));
    PHYSIM_LOG_INFO("  Velocity: [{:.3f}, {:.3f}, {:.3f}] m/s", v0(0), v0(1), v0(2));
    PHYSIM_LOG_INFO("  Distance: {:.3e} km", r0.norm() / 1000.0);
    PHYSIM_LOG_INFO("  Speed: {:.3f} m/s", v0.norm());

    // Compute initial orbital elements
    auto elem0 = compute_elements(r0, v0, mu_earth);
    PHYSIM_LOG_INFO("Initial orbital elements:");
    PHYSIM_LOG_INFO("  a = {:.3e} km", elem0.a / 1000.0);
    PHYSIM_LOG_INFO("  e = {:.6f}", elem0.e);
    PHYSIM_LOG_INFO("  i = {:.3f}°", elem0.i * constants::RAD_TO_DEG);

    // Create gravity force
    auto earth_gravity = forces::PointMassGravity::for_body("Earth");

    // Setup derivative function
    auto derivative = integrators::newtonian_derivative(
        [&earth_gravity](double t, const Vec3& pos, const Vec3& vel) {
            return earth_gravity->acceleration(t, pos, vel, 1.0);
        });

    // Initial state vector
    VecX y0 = integrators::state_to_vector_6dof(State(r0, v0, 1.0));

    // Simulation parameters
    const double orbital_period = 2.0 * constants::PI * std::sqrt(elem0.a * elem0.a * elem0.a / mu_earth);
    const double t0 = 0.0;
    const double tf = orbital_period;  // One complete orbit
    const double dt_rk4 = 60.0;        // 60 second steps for RK4
    const double dt_rk45 = 3600.0;     // 1 hour initial step for RK45

    PHYSIM_LOG_INFO("Orbital period: {:.2f} days", orbital_period / constants::SECONDS_PER_DAY);
    PHYSIM_LOG_INFO("Simulation duration: {:.2f} days", tf / constants::SECONDS_PER_DAY);

    // ========================================================================
    // Integration with RK4 (fixed step)
    // ========================================================================

    PHYSIM_LOG_INFO("\n--- RK4 Integration (fixed step: {}s) ---", dt_rk4);

    integrators::RK4 rk4;
    auto y_rk4 = rk4.integrate(t0, y0, tf, dt_rk4, derivative, false);

    const auto& stats_rk4 = rk4.stats();
    PHYSIM_LOG_INFO("RK4 Statistics:");
    PHYSIM_LOG_INFO("  Steps: {}", stats_rk4.num_steps);
    PHYSIM_LOG_INFO("  Function evaluations: {}", stats_rk4.num_function_evals);

    // Extract final state
    State final_state_rk4;
    integrators::vector_to_state_6dof(y_rk4, final_state_rk4);

    PHYSIM_LOG_INFO("RK4 Final state:");
    PHYSIM_LOG_INFO("  Position: [{:.3e}, {:.3e}, {:.3e}] m",
                   final_state_rk4.position(0), final_state_rk4.position(1), final_state_rk4.position(2));
    PHYSIM_LOG_INFO("  Velocity: [{:.3f}, {:.3f}, {:.3f}] m/s",
                   final_state_rk4.velocity(0), final_state_rk4.velocity(1), final_state_rk4.velocity(2));

    // Position error after one orbit
    const double pos_error_rk4 = (final_state_rk4.position - r0).norm();
    const double vel_error_rk4 = (final_state_rk4.velocity - v0).norm();

    PHYSIM_LOG_INFO("RK4 Error after one orbit:");
    PHYSIM_LOG_INFO("  Position error: {:.3f} km", pos_error_rk4 / 1000.0);
    PHYSIM_LOG_INFO("  Velocity error: {:.6f} m/s", vel_error_rk4);

    // Energy conservation check
    const double E0 = v0.squaredNorm() / 2.0 - mu_earth / r0.norm();
    const double E_final_rk4 = final_state_rk4.velocity.squaredNorm() / 2.0 -
                               mu_earth / final_state_rk4.position.norm();
    const double energy_error_rk4 = std::abs((E_final_rk4 - E0) / E0);

    PHYSIM_LOG_INFO("  Energy conservation (relative): {:.3e}", energy_error_rk4);

    // ========================================================================
    // Integration with RK45 (adaptive step)
    // ========================================================================

    PHYSIM_LOG_INFO("\n--- RK45 Integration (adaptive, initial dt: {}s, tol: 1e-12) ---",
                   dt_rk45);

    integrators::RK45 rk45;
    auto y_rk45 = rk45.integrate(t0, y0, tf, dt_rk45, derivative, true, 1e-12);

    const auto& stats_rk45 = rk45.stats();
    PHYSIM_LOG_INFO("RK45 Statistics:");
    PHYSIM_LOG_INFO("  Steps: {}", stats_rk45.num_steps);
    PHYSIM_LOG_INFO("  Function evaluations: {}", stats_rk45.num_function_evals);
    PHYSIM_LOG_INFO("  Rejected steps: {}", stats_rk45.num_rejected_steps);
    PHYSIM_LOG_INFO("  Final step size: {:.3f} s", stats_rk45.final_step_size);
    PHYSIM_LOG_INFO("  Max error: {:.3e}", stats_rk45.max_error);

    // Extract final state
    State final_state_rk45;
    integrators::vector_to_state_6dof(y_rk45, final_state_rk45);

    PHYSIM_LOG_INFO("RK45 Final state:");
    PHYSIM_LOG_INFO("  Position: [{:.3e}, {:.3e}, {:.3e}] m",
                   final_state_rk45.position(0), final_state_rk45.position(1), final_state_rk45.position(2));
    PHYSIM_LOG_INFO("  Velocity: [{:.3f}, {:.3f}, {:.3f}] m/s",
                   final_state_rk45.velocity(0), final_state_rk45.velocity(1), final_state_rk45.velocity(2));

    // Position error after one orbit
    const double pos_error_rk45 = (final_state_rk45.position - r0).norm();
    const double vel_error_rk45 = (final_state_rk45.velocity - v0).norm();

    PHYSIM_LOG_INFO("RK45 Error after one orbit:");
    PHYSIM_LOG_INFO("  Position error: {:.3f} km", pos_error_rk45 / 1000.0);
    PHYSIM_LOG_INFO("  Velocity error: {:.6f} m/s", vel_error_rk45);

    // Energy conservation check
    const double E_final_rk45 = final_state_rk45.velocity.squaredNorm() / 2.0 -
                               mu_earth / final_state_rk45.position.norm();
    const double energy_error_rk45 = std::abs((E_final_rk45 - E0) / E0);

    PHYSIM_LOG_INFO("  Energy conservation (relative): {:.3e}", energy_error_rk45);

    // ========================================================================
    // Comparison
    // ========================================================================

    PHYSIM_LOG_INFO("\n--- Comparison ---");
    PHYSIM_LOG_INFO("RK4:  {} steps, {} f-evals, {:.3f} km error, {:.3e} energy drift",
                   stats_rk4.num_steps, stats_rk4.num_function_evals,
                   pos_error_rk4 / 1000.0, energy_error_rk4);
    PHYSIM_LOG_INFO("RK45: {} steps, {} f-evals, {:.3f} km error, {:.3e} energy drift",
                   stats_rk45.num_steps, stats_rk45.num_function_evals,
                   pos_error_rk45 / 1000.0, energy_error_rk45);

    const double speedup = static_cast<double>(stats_rk4.num_function_evals) /
                          static_cast<double>(stats_rk45.num_function_evals);
    PHYSIM_LOG_INFO("RK45 efficiency: {:.2f}x fewer function evaluations", speedup);

    // Success criteria
    const bool success = (pos_error_rk45 < 1000.0) &&  // < 1 km error
                         (energy_error_rk45 < 1e-9);     // < 1e-9 relative

    if (success) {
        PHYSIM_LOG_INFO("\n✓ Example completed successfully!");
        PHYSIM_LOG_INFO("  Position error: {:.3f} km (< 1 km target)", pos_error_rk45 / 1000.0);
        PHYSIM_LOG_INFO("  Energy drift: {:.3e} (< 1e-9 target)", energy_error_rk45);
    } else {
        PHYSIM_LOG_ERROR("\n✗ Example failed accuracy criteria");
    }

    logging::shutdown();

    return success ? 0 : 1;
}
