/**
 * @file inner_solar_system.cpp
 * @brief Inner solar system N-body simulation
 *
 * Simulates Sun, Mercury, Venus, Earth, and Mars using direct N-body integration.
 * Demonstrates multi-body dynamics with realistic orbital parameters.
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
    PHYSIM_LOG_INFO("=== Inner Solar System N-Body Simulation ===");

    // Create N-body system
    NBodySystem system;

    // Add Sun at origin with zero velocity
    auto sun = create_solar_system_particle("Sun", Vec3::Zero(), Vec3::Zero());
    system.add_particle(sun);

    // Simplified orbital parameters (circular orbits for demonstration)
    // In reality, orbits are elliptical - this is just for testing

    // Mercury: 57.9 million km, orbital period ~88 days
    const double r_mercury = 57.9e9;  // m
    const double v_mercury = std::sqrt(constants::gm::SUN / r_mercury);
    auto mercury = create_solar_system_particle("Mercury",
                                                Vec3(r_mercury, 0, 0),
                                                Vec3(0, v_mercury, 0));
    system.add_particle(mercury);

    // Venus: 108.2 million km, orbital period ~225 days
    const double r_venus = 108.2e9;
    const double v_venus = std::sqrt(constants::gm::SUN / r_venus);
    auto venus = create_solar_system_particle("Venus",
                                             Vec3(r_venus, 0, 0),
                                             Vec3(0, v_venus, 0));
    system.add_particle(venus);

    // Earth: 149.6 million km (1 AU), orbital period 365.25 days
    const double r_earth = constants::ASTRONOMICAL_UNIT;
    const double v_earth = std::sqrt(constants::gm::SUN / r_earth);
    auto earth = create_solar_system_particle("Earth",
                                              Vec3(r_earth, 0, 0),
                                              Vec3(0, v_earth, 0));
    system.add_particle(earth);

    // Mars: 227.9 million km, orbital period ~687 days
    const double r_mars = 227.9e9;
    const double v_mars = std::sqrt(constants::gm::SUN / r_mars);
    auto mars = create_solar_system_particle("Mars",
                                             Vec3(r_mars, 0, 0),
                                             Vec3(0, v_mars, 0));
    system.add_particle(mars);

    // Initial statistics
    auto stats0 = system.compute_stats();
    system.reset_energy_tracking();

    PHYSIM_LOG_INFO("\nSystem initialized with {} bodies", system.size());
    PHYSIM_LOG_INFO("  Total mass: {:.6e} kg ({:.3f}% is Sun)",
                   stats0.total_mass,
                   100.0 * sun->mass() / stats0.total_mass);
    PHYSIM_LOG_INFO("  Total energy: {:.6e} J", stats0.total_energy);
    PHYSIM_LOG_INFO("  Kinetic energy: {:.6e} J", stats0.kinetic_energy);
    PHYSIM_LOG_INFO("  Potential energy: {:.6e} J", stats0.potential_energy);

    // Print orbital velocities
    PHYSIM_LOG_INFO("\nOrbital velocities:");
    PHYSIM_LOG_INFO("  Mercury: {:.3f} km/s", v_mercury / 1000.0);
    PHYSIM_LOG_INFO("  Venus:   {:.3f} km/s", v_venus / 1000.0);
    PHYSIM_LOG_INFO("  Earth:   {:.3f} km/s", v_earth / 1000.0);
    PHYSIM_LOG_INFO("  Mars:    {:.3f} km/s", v_mars / 1000.0);

    // Simulate for 1 Earth year
    const double t_final = 365.25 * constants::SECONDS_PER_DAY;
    const double dt = 86400.0;  // 1 day time step
    const double tolerance = 1e-10;

    PHYSIM_LOG_INFO("\nSimulation parameters:");
    PHYSIM_LOG_INFO("  Duration: 1 Earth year (365.25 days)");
    PHYSIM_LOG_INFO("  Initial dt: {:.1f} days", dt / constants::SECONDS_PER_DAY);
    PHYSIM_LOG_INFO("  Tolerance: {:.0e}", tolerance);

    PHYSIM_LOG_INFO("\nStarting propagation...");

    integrators::RK45 integrator;
    system.propagate(t_final, dt, integrator, true, tolerance);

    PHYSIM_LOG_INFO("Propagation complete!");

    // Final statistics
    auto stats_f = system.compute_stats();

    PHYSIM_LOG_INFO("\nIntegration statistics:");
    PHYSIM_LOG_INFO("  Steps: {}", integrator.stats().num_steps);
    PHYSIM_LOG_INFO("  Function evaluations: {}", integrator.stats().num_function_evals);
    PHYSIM_LOG_INFO("  Rejected steps: {}", integrator.stats().num_rejected_steps);
    PHYSIM_LOG_INFO("  Final step size: {:.3f} days",
                   integrator.stats().final_step_size / constants::SECONDS_PER_DAY);

    PHYSIM_LOG_INFO("\nFinal system state:");
    PHYSIM_LOG_INFO("  Total energy: {:.6e} J", stats_f.total_energy);
    PHYSIM_LOG_INFO("  Energy error: {:.3e} (relative)", stats_f.energy_error);

    // Print planet positions after 1 year
    PHYSIM_LOG_INFO("\nPlanet positions after 1 Earth year:");
    for (size_t i = 0; i < system.size(); ++i) {
        const auto& p = system.particle(i);
        const double r = p.position().norm();
        PHYSIM_LOG_INFO("  {}: {:.3e} km (initial: {:.3e} km)",
                       p.name(),
                       r / 1000.0,
                       (i == 0 ? 0.0 : (i == 1 ? r_mercury : i == 2 ? r_venus : i == 3 ? r_earth : r_mars)) / 1000.0);
    }

    // Earth should be back near its starting position after 1 year
    const nbody::Particle& earth_final = system.particle(3);
    const double earth_pos_error = (earth_final.position() - Vec3(r_earth, 0, 0)).norm();

    PHYSIM_LOG_INFO("\nEarth position error after 1 year: {:.3e} km", earth_pos_error / 1000.0);

    // Conservation checks
    PHYSIM_LOG_INFO("\nConservation verification:");
    PHYSIM_LOG_INFO("  Energy conservation: {:.3e} (relative error)", stats_f.energy_error);
    PHYSIM_LOG_INFO("  Momentum: [{:.3e}, {:.3e}, {:.3e}] kg·m/s",
                   stats_f.total_momentum(0), stats_f.total_momentum(1), stats_f.total_momentum(2));

    // Success criteria (looser for 5-body problem)
    const bool success = (stats_f.energy_error < 1e-6) &&
                         (earth_pos_error < 1e9);  // < 1 million km

    if (success) {
        PHYSIM_LOG_INFO("\n✓ Solar system simulation completed successfully");
        PHYSIM_LOG_INFO("  Energy conservation: {:.3e}", stats_f.energy_error);
        PHYSIM_LOG_INFO("  Earth orbital accuracy: {:.3e} km", earth_pos_error / 1000.0);
    } else {
        PHYSIM_LOG_WARN("\n⚠ Solar system simulation completed with warnings");
        if (stats_f.energy_error >= 1e-6) {
            PHYSIM_LOG_WARN("  Energy error: {:.3e} (expected < 1e-6)", stats_f.energy_error);
        }
    }

    logging::shutdown();

    return success ? 0 : 1;
}
