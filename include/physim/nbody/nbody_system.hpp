/**
 * @file nbody_system.hpp
 * @brief N-body gravitational system with direct summation
 *
 * Manages collections of particles and propagates them forward in time
 * using direct O(N²) force summation with configurable integrators.
 */

#pragma once

#include "physim/nbody/particle.hpp"
#include "physim/integrators/integrator.hpp"
#include "physim/core/types.hpp"
#include <vector>
#include <memory>
#include <functional>

namespace physim {
namespace nbody {

/**
 * @brief Statistics for N-body system
 */
struct NBodyStats {
    double total_mass;           ///< Total system mass [kg]
    double total_energy;         ///< Total energy (kinetic + potential) [J]
    double kinetic_energy;       ///< Total kinetic energy [J]
    double potential_energy;     ///< Total gravitational potential energy [J]
    Vec3 center_of_mass;         ///< Center of mass position [m]
    Vec3 total_momentum;         ///< Total linear momentum [kg⋅m/s]
    Vec3 total_angular_momentum; ///< Total angular momentum [kg⋅m²/s]
    double energy_error;         ///< Relative energy error since t=0
    size_t num_particles;        ///< Number of particles in system
    double current_time;         ///< Current simulation time [s]
};

/**
 * @brief N-body gravitational system using direct summation
 *
 * Manages a collection of gravitationally interacting particles and
 * propagates them forward in time using direct O(N²) force calculation.
 * Integrates with RK4, RK45, or other ODE integrators.
 */
class NBodySystem {
public:
    /**
     * @brief Construct empty N-body system
     * @param use_softening Enable gravitational softening (default: false)
     * @param softening_length Softening length ε [m] (default: 0)
     */
    explicit NBodySystem(bool use_softening = false, double softening_length = 0.0);

    /**
     * @brief Add particle to system
     * @param particle Particle to add
     */
    void add_particle(const Particle& particle);

    /**
     * @brief Add particle to system (shared_ptr version)
     * @param particle Particle to add
     */
    void add_particle(std::shared_ptr<Particle> particle);

    /**
     * @brief Remove particle at index
     * @param index Particle index
     */
    void remove_particle(size_t index);

    /**
     * @brief Clear all particles
     */
    void clear();

    /**
     * @brief Get number of particles
     * @return Number of particles in system
     */
    size_t size() const { return particles_.size(); }

    /**
     * @brief Check if system is empty
     * @return true if no particles
     */
    bool empty() const { return particles_.empty(); }

    /**
     * @brief Get particle at index (const)
     * @param index Particle index
     * @return Reference to particle
     */
    const Particle& particle(size_t index) const;

    /**
     * @brief Get particle at index (mutable)
     * @param index Particle index
     * @return Reference to particle
     */
    Particle& particle(size_t index);

    /**
     * @brief Get all particles (const)
     * @return Vector of particles
     */
    const std::vector<Particle>& particles() const { return particles_; }

    /**
     * @brief Propagate system forward in time
     * @param dt Time step [s]
     * @param integrator Integrator to use (default: RK4)
     * @param adaptive Use adaptive stepping if integrator supports it
     * @param tolerance Error tolerance for adaptive stepping
     */
    void step(double dt,
              integrators::Integrator& integrator,
              bool adaptive = false,
              double tolerance = 1e-12);

    /**
     * @brief Propagate system over time interval
     * @param t_final Final time [s]
     * @param dt Time step [s]
     * @param integrator Integrator to use
     * @param adaptive Use adaptive stepping
     * @param tolerance Error tolerance for adaptive stepping
     */
    void propagate(double t_final,
                   double dt,
                   integrators::Integrator& integrator,
                   bool adaptive = false,
                   double tolerance = 1e-12);

    /**
     * @brief Compute accelerations for all particles (direct summation)
     *
     * Updates acceleration for each particle using Newton's law of gravitation:
     * a_i = Σ(j≠i) G * m_j * (r_j - r_i) / |r_j - r_i|³
     *
     * Complexity: O(N²)
     */
    void compute_accelerations();

    /**
     * @brief Compute total kinetic energy
     * @return Sum of (1/2)m_i v_i² [J]
     */
    double kinetic_energy() const;

    /**
     * @brief Compute total gravitational potential energy
     * @return Sum of -G m_i m_j / r_ij [J]
     */
    double potential_energy() const;

    /**
     * @brief Compute total energy (kinetic + potential)
     * @return Total energy [J]
     */
    double total_energy() const;

    /**
     * @brief Compute center of mass position
     * @return r_cm = Σ(m_i r_i) / Σ(m_i) [m]
     */
    Vec3 center_of_mass() const;

    /**
     * @brief Compute center of mass velocity
     * @return v_cm = Σ(m_i v_i) / Σ(m_i) [m/s]
     */
    Vec3 center_of_mass_velocity() const;

    /**
     * @brief Compute total linear momentum
     * @return p = Σ(m_i v_i) [kg⋅m/s]
     */
    Vec3 total_momentum() const;

    /**
     * @brief Compute total angular momentum about origin
     * @return L = Σ(r_i × p_i) [kg⋅m²/s]
     */
    Vec3 total_angular_momentum() const;

    /**
     * @brief Compute total angular momentum about arbitrary point
     * @param origin Reference point
     * @return L = Σ((r_i - r_0) × p_i) [kg⋅m²/s]
     */
    Vec3 total_angular_momentum(const Vec3& origin) const;

    /**
     * @brief Compute total mass
     * @return M = Σ(m_i) [kg]
     */
    double total_mass() const;

    /**
     * @brief Get comprehensive system statistics
     * @return NBodyStats structure
     */
    NBodyStats compute_stats() const;

    /**
     * @brief Get current simulation time
     * @return Current time [s]
     */
    double time() const { return current_time_; }

    /**
     * @brief Set current simulation time
     * @param t Time [s]
     */
    void set_time(double t) { current_time_ = t; }

    /**
     * @brief Reset energy tracking (sets initial energy for error calculation)
     */
    void reset_energy_tracking();

    /**
     * @brief Enable/disable gravitational softening
     * @param enable Enable softening
     * @param epsilon Softening length [m]
     */
    void set_softening(bool enable, double epsilon = 0.0);

private:
    /**
     * @brief Convert system state to ODE state vector
     *
     * State vector layout: [x1, y1, z1, vx1, vy1, vz1, x2, y2, z2, ...]
     * Size: 6N elements
     */
    VecX to_state_vector() const;

    /**
     * @brief Update system from ODE state vector
     * @param state_vec State vector of size 6N
     */
    void from_state_vector(const VecX& state_vec);

    /**
     * @brief ODE derivative function for integration
     *
     * Computes dy/dt where y = [positions, velocities]
     * Returns dydt = [velocities, accelerations]
     */
    integrators::DerivativeFunction derivative_function();

    std::vector<Particle> particles_;  ///< Particle storage
    double current_time_;              ///< Current simulation time [s]
    double initial_energy_;            ///< Initial energy for error tracking
    bool energy_tracking_reset_;       ///< Whether energy tracking initialized
    bool use_softening_;               ///< Enable softening
    double softening_length_;          ///< Softening parameter ε [m]
};

} // namespace nbody
} // namespace physim
