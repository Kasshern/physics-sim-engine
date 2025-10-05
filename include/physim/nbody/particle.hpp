/**
 * @file particle.hpp
 * @brief Particle class for N-body simulations
 *
 * Simplified point mass representation for gravitational N-body problems.
 */

#pragma once

#include "physim/core/types.hpp"
#include <string>
#include <memory>

namespace physim {
namespace nbody {

/**
 * @brief Point mass particle for N-body simulations
 *
 * Represents a gravitational point mass with position, velocity, and mass.
 * Designed for efficient N-body force calculations.
 */
class Particle {
public:
    /**
     * @brief Default constructor - creates particle at origin with zero mass
     */
    Particle();

    /**
     * @brief Construct particle with initial state
     * @param position Initial position [m]
     * @param velocity Initial velocity [m/s]
     * @param mass Particle mass [kg]
     * @param name Optional particle name
     */
    Particle(const Vec3& position, const Vec3& velocity, double mass,
             const std::string& name = "");

    /**
     * @brief Construct particle from State object
     * @param state State containing position, velocity, mass
     * @param name Optional particle name
     */
    explicit Particle(const State& state, const std::string& name = "");

    // Accessors
    const Vec3& position() const { return position_; }
    const Vec3& velocity() const { return velocity_; }
    const Vec3& acceleration() const { return acceleration_; }
    double mass() const { return mass_; }
    const std::string& name() const { return name_; }

    // Mutators
    void set_position(const Vec3& pos) { position_ = pos; }
    void set_velocity(const Vec3& vel) { velocity_ = vel; }
    void set_acceleration(const Vec3& acc) { acceleration_ = acc; }
    void set_mass(double m);
    void set_name(const std::string& n) { name_ = n; }

    /**
     * @brief Compute kinetic energy
     * @return Kinetic energy T = (1/2)mv² [J]
     */
    double kinetic_energy() const;

    /**
     * @brief Compute linear momentum
     * @return Momentum p = mv [kg⋅m/s]
     */
    Vec3 momentum() const;

    /**
     * @brief Compute angular momentum about origin
     * @return Angular momentum L = r × p [kg⋅m²/s]
     */
    Vec3 angular_momentum() const;

    /**
     * @brief Compute angular momentum about arbitrary point
     * @param origin Reference point
     * @return Angular momentum L = (r - r₀) × p [kg⋅m²/s]
     */
    Vec3 angular_momentum(const Vec3& origin) const;

    /**
     * @brief Get distance to another particle
     * @param other Other particle
     * @return Distance between particles [m]
     */
    double distance_to(const Particle& other) const;

    /**
     * @brief Get squared distance to another particle (faster)
     * @param other Other particle
     * @return Squared distance [m²]
     */
    double distance_squared_to(const Particle& other) const;

    /**
     * @brief Convert particle state to State object
     * @return State object with position, velocity, mass
     */
    State to_state() const;

    /**
     * @brief Update particle state from State object
     * @param state New state
     */
    void from_state(const State& state);

    /**
     * @brief Convert to 6-DOF state vector [x, y, z, vx, vy, vz]
     * @return 6-element state vector
     */
    VecX to_vector() const;

    /**
     * @brief Update from 6-DOF state vector
     * @param state_vec 6-element vector [x, y, z, vx, vy, vz]
     */
    void from_vector(const VecX& state_vec);

    /**
     * @brief Apply velocity update (for integration)
     * @param dv Change in velocity [m/s]
     */
    void apply_velocity_update(const Vec3& dv) { velocity_ += dv; }

    /**
     * @brief Apply position update (for integration)
     * @param dr Change in position [m]
     */
    void apply_position_update(const Vec3& dr) { position_ += dr; }

private:
    Vec3 position_;       ///< Position [m]
    Vec3 velocity_;       ///< Velocity [m/s]
    Vec3 acceleration_;   ///< Acceleration [m/s²] (computed by NBodySystem)
    double mass_;         ///< Mass [kg]
    std::string name_;    ///< Optional particle name
};

/**
 * @brief Named particle factory for solar system bodies
 * @param body_name Name of celestial body (e.g., "Sun", "Earth")
 * @param position Initial position [m]
 * @param velocity Initial velocity [m/s]
 * @return Particle with correct mass
 */
std::shared_ptr<Particle> create_solar_system_particle(
    const std::string& body_name,
    const Vec3& position = Vec3::Zero(),
    const Vec3& velocity = Vec3::Zero());

} // namespace nbody
} // namespace physim
