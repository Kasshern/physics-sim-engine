/**
 * @file two_body_nbody.cpp
 * @brief Two-body problem validation using N-body system
 *
 * Validates N-body implementation against analytical two-body solution.
 * Tests energy and momentum conservation for Earth-Moon system.
 */

#include "physim/nbody/nbody_system.hpp"
#include "physim/nbody/particle.hpp"
#include "physim/integrators/rk45.hpp"
#include "physim/core/constants.hpp"
#include "physim/core/logging.hpp"

#include <iostream>
#include <iomanip>

using namespace physim;
using namespace physim::nbody;

int main() {
    // Initialize logging
    logging::init(spdlog::level::info, false);
    PHYSIM_LOG_INFO("=== Two-Body N-Body Validation ===");

    // Create N-body system
    NBodySystem system;

    // Earth-Moon system parameters
    const double m_earth = constants::gm::EARTH / constants::GRAVITATIONAL_CONSTANT;
    const double m_moon = constants::gm::MOON / constants::GRAVITATIONAL_CONSTANT;

    PHYSIM_LOG_INFO("System masses:");
    PHYSIM_LOG_INFO("  Earth: {:.6e} kg", m_earth);
    PHYSIM_LOG_INFO("  Moon:  {:.6e} kg", m_moon);

    // Orbital parameters (Moon at perigee)
    const double a_moon = 384400e3;  // Semi-major axis [m]
    const double e_moon = 0.0549;     // Eccentricity
    const double r_perigee = a_moon * (1.0 - e_moon);

    // Calculate barycentric positions
    const double total_mass = m_earth + m_moon;
    const double r_earth_bary = r_perigee * m_moon / total_mass;
    const double r_moon_bary = r_perigee * m_earth / total_mass;

    PHYSIM_LOG_INFO("Barycentric separation: {:.3e} km", r_perigee / 1000.0);
    PHYSIM_LOG_INFO("  Earth from barycenter: {:.3e} km", r_earth_bary / 1000.0);
    PHYSIM_LOG_INFO("  Moon from barycenter:  {:.3e} km", r_moon_bary / 1000.0);

    // Circular velocity at perigee
    const double mu_total = constants::gm::EARTH + constants::gm::MOON;
    const double v_perigee = std::sqrt(mu_total * (1.0 + e_moon) / r_perigee);

    // Barycentric velocities
    const double v_earth_bary = v_perigee * m_moon / total_mass;
    const double v_moon_bary = v_perigee * m_earth / total_mass;

    PHYSIM_LOG_INFO("Orbital velocity at perigee: {:.3f} m/s", v_perigee);
    PHYSIM_LOG_INFO("  Earth: {:.3f} m/s", v_earth_bary);
    PHYSIM_LOG_INFO("  Moon:  {:.3f} m/s", v_moon_bary);

    // Create particles at barycentric frame
    nbody::Particle earth(Vec3(-r_earth_bary, 0, 0), Vec3(0, -v_earth_bary, 0), m_earth, "Earth");
    nbody::Particle moon(Vec3(r_moon_bary, 0, 0), Vec3(0, v_moon_bary, 0), m_moon, "Moon");

    system.add_particle(earth);
    system.add_particle(moon);

    // Initial statistics
    auto stats0 = system.compute_stats();
    system.reset_energy_tracking();

    PHYSIM_LOG_INFO("\nInitial system state:");
    PHYSIM_LOG_INFO("  Total mass: {:.6e} kg", stats0.total_mass);
    PHYSIM_LOG_INFO("  Total energy: {:.6e} J", stats0.total_energy);
    PHYSIM_LOG_INFO("  Kinetic energy: {:.6e} J", stats0.kinetic_energy);
    PHYSIM_LOG_INFO("  Potential energy: {:.6e} J", stats0.potential_energy);
    PHYSIM_LOG_INFO("  Center of mass: [{:.3e}, {:.3e}, {:.3e}] m",
                   stats0.center_of_mass(0), stats0.center_of_mass(1), stats0.center_of_mass(2));
    PHYSIM_LOG_INFO("  Total momentum: [{:.3e}, {:.3e}, {:.3e}] kg·m/s",
                   stats0.total_momentum(0), stats0.total_momentum(1), stats0.total_momentum(2));
    PHYSIM_LOG_INFO("  Total angular momentum: [{:.3e}, {:.3e}, {:.3e}] kg·m²/s",
                   stats0.total_angular_momentum(0), stats0.total_angular_momentum(1), stats0.total_angular_momentum(2));

    // Compute orbital period
    const double period = 2.0 * constants::PI * std::sqrt(a_moon * a_moon * a_moon / mu_total);
    PHYSIM_LOG_INFO("\nOrbital period: {:.3f} days", period / constants::SECONDS_PER_DAY);

    // Propagate for one full orbit
    const double t_final = period;
    const double dt = 3600.0;  // 1 hour time step
    const double tolerance = 1e-12;

    PHYSIM_LOG_INFO("\nPropagating for one orbital period...");
    PHYSIM_LOG_INFO("  Duration: {:.3f} days", t_final / constants::SECONDS_PER_DAY);
    PHYSIM_LOG_INFO("  Initial dt: {:.1f} s", dt);
    PHYSIM_LOG_INFO("  Tolerance: {:.0e}", tolerance);

    integrators::RK45 integrator;
    system.propagate(t_final, dt, integrator, true, tolerance);

    // Final statistics
    auto stats_f = system.compute_stats();

    PHYSIM_LOG_INFO("\nIntegration statistics:");
    PHYSIM_LOG_INFO("  Steps: {}", integrator.stats().num_steps);
    PHYSIM_LOG_INFO("  Function evaluations: {}", integrator.stats().num_function_evals);
    PHYSIM_LOG_INFO("  Rejected steps: {}", integrator.stats().num_rejected_steps);

    PHYSIM_LOG_INFO("\nFinal system state:");
    PHYSIM_LOG_INFO("  Total energy: {:.6e} J", stats_f.total_energy);
    PHYSIM_LOG_INFO("  Energy error: {:.3e} (relative)", stats_f.energy_error);
    PHYSIM_LOG_INFO("  Center of mass: [{:.3e}, {:.3e}, {:.3e}] m",
                   stats_f.center_of_mass(0), stats_f.center_of_mass(1), stats_f.center_of_mass(2));
    PHYSIM_LOG_INFO("  Total momentum: [{:.3e}, {:.3e}, {:.3e}] kg·m/s",
                   stats_f.total_momentum(0), stats_f.total_momentum(1), stats_f.total_momentum(2));

    // Verify conservation laws
    const double com_drift = stats_f.center_of_mass.norm();
    const double momentum_error = stats_f.total_momentum.norm();

    PHYSIM_LOG_INFO("\nConservation law verification:");
    PHYSIM_LOG_INFO("  Energy conservation: {:.3e} (relative error)", stats_f.energy_error);
    PHYSIM_LOG_INFO("  Center of mass drift: {:.3e} m", com_drift);
    PHYSIM_LOG_INFO("  Linear momentum error: {:.3e} kg·m/s", momentum_error);

    // Verify particles returned to near-initial positions
    const nbody::Particle& earth_f = system.particle(0);
    const nbody::Particle& moon_f = system.particle(1);

    const double pos_error_earth = (earth_f.position() - earth.position()).norm();
    const double pos_error_moon = (moon_f.position() - moon.position()).norm();

    PHYSIM_LOG_INFO("\nPosition errors after one orbit:");
    PHYSIM_LOG_INFO("  Earth: {:.3f} km", pos_error_earth / 1000.0);
    PHYSIM_LOG_INFO("  Moon:  {:.3f} km", pos_error_moon / 1000.0);

    // Success criteria
    const bool success = (stats_f.energy_error < 1e-8) &&
                         (com_drift < 1000.0) &&
                         (momentum_error < 1e-3) &&
                         (pos_error_moon < 10000.0);  // < 10 km

    if (success) {
        PHYSIM_LOG_INFO("\n✓ Two-body validation PASSED");
        PHYSIM_LOG_INFO("  Energy conservation: excellent ({:.3e})", stats_f.energy_error);
        PHYSIM_LOG_INFO("  Momentum conservation: excellent ({:.3e})", momentum_error);
        PHYSIM_LOG_INFO("  Position accuracy: {:.3f} km", pos_error_moon / 1000.0);
    } else {
        PHYSIM_LOG_ERROR("\n✗ Two-body validation FAILED");
        if (stats_f.energy_error >= 1e-8) {
            PHYSIM_LOG_ERROR("  Energy error too large: {:.3e}", stats_f.energy_error);
        }
        if (pos_error_moon >= 10000.0) {
            PHYSIM_LOG_ERROR("  Position error too large: {:.3f} km", pos_error_moon / 1000.0);
        }
    }

    logging::shutdown();

    return success ? 0 : 1;
}
