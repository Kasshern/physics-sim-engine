/**
 * @file gravity.hpp
 * @brief Point mass gravitational force
 *
 * Implements Newton's law of universal gravitation for point masses.
 * Used for N-body simulations and basic orbital mechanics.
 */

#pragma once

#include "physim/forces/force.hpp"
#include "physim/core/constants.hpp"

namespace physim {
namespace forces {

/**
 * @brief Point mass gravitational force
 *
 * Computes gravitational acceleration using Newton's law:
 *   a = -GM/r² * r̂
 *
 * where:
 * - G is the gravitational constant
 * - M is the mass of the attracting body
 * - r is the distance between bodies
 * - r̂ is the unit vector from attracting body to test particle
 *
 * This is the fundamental force in orbital mechanics and N-body simulations.
 *
 * Properties:
 * - Conservative (potential energy: U = -GMm/r)
 * - Velocity-independent
 * - Time-independent (for fixed central body)
 * - Singular at r = 0 (requires regularization for close encounters)
 *
 * Limitations:
 * - Assumes point masses (no extended body effects)
 * - Ignores relativistic corrections (use GR for compact objects)
 * - Ignores aspherical gravity harmonics (use J2Force for Earth)
 *
 * Use cases:
 * - Two-body problem (planet orbiting star)
 * - N-body gravitational dynamics
 * - Solar system simulation
 * - Interplanetary trajectories
 */
class PointMassGravity : public Force {
public:
    /**
     * @brief Constructor
     * @param GM Gravitational parameter μ = GM [m³/s²]
     * @param center_pos Position of attracting body [m] (default: origin)
     * @param center_vel Velocity of attracting body [m/s] (default: stationary)
     *
     * For most applications (planet orbiting star), center is at origin.
     * For N-body simulations, center can move.
     */
    explicit PointMassGravity(double GM,
                             const Vec3& center_pos = Vec3::Zero(),
                             const Vec3& center_vel = Vec3::Zero());

    /**
     * @brief Convenience constructor for solar system bodies
     * @param body_name Name of body ("sun", "earth", "moon", etc.)
     * @param center_pos Position of body [m]
     * @param center_vel Velocity of body [m/s]
     * @return Shared pointer to configured gravity force
     */
    static std::shared_ptr<PointMassGravity> for_body(
        const std::string& body_name,
        const Vec3& center_pos = Vec3::Zero(),
        const Vec3& center_vel = Vec3::Zero());

    /**
     * @brief Compute gravitational acceleration
     *
     * @param t Current time [s] (unused, gravity is time-independent)
     * @param pos Position of test particle [m]
     * @param vel Velocity of test particle [m/s] (unused)
     * @param mass Mass of test particle [kg] (unused, cancels in a = F/m)
     * @return Acceleration [m/s²]
     *
     * Computes: a = -μ/r² * r̂ where r̂ = (pos - center_pos)/r
     */
    Vec3 acceleration(double t, const Vec3& pos, const Vec3& vel,
                     double mass) const override;

    /**
     * @brief Compute gravitational potential energy
     *
     * @param t Current time [s] (unused)
     * @param pos Position of test particle [m]
     * @param mass Mass of test particle [kg]
     * @return Potential energy U = -GMm/r [J]
     */
    double potential_energy(double t, const Vec3& pos, double mass) const override;

    /**
     * @brief Check if velocity-dependent (always false)
     */
    bool is_velocity_dependent() const override { return false; }

    /**
     * @brief Check if time-dependent (false for stationary center)
     */
    bool is_time_dependent() const override { return center_is_moving_; }

    /**
     * @brief Get gravitational parameter μ = GM
     */
    double mu() const { return mu_; }

    /**
     * @brief Get center position
     */
    const Vec3& center_position() const { return center_pos_; }

    /**
     * @brief Get center velocity
     */
    const Vec3& center_velocity() const { return center_vel_; }

    /**
     * @brief Set center position (for moving bodies in N-body)
     * @param pos New center position [m]
     */
    void set_center_position(const Vec3& pos) {
        center_pos_ = pos;
    }

    /**
     * @brief Set center velocity (for moving bodies in N-body)
     * @param vel New center velocity [m/s]
     */
    void set_center_velocity(const Vec3& vel) {
        center_vel_ = vel;
        center_is_moving_ = (vel.norm() > 1e-10);
    }

    /**
     * @brief Compute orbital velocity for circular orbit at radius r
     * @param r Orbital radius [m]
     * @return Circular orbit velocity [m/s]: v = sqrt(μ/r)
     */
    double circular_velocity(double r) const;

    /**
     * @brief Compute escape velocity at radius r
     * @param r Radius [m]
     * @return Escape velocity [m/s]: v = sqrt(2μ/r)
     */
    double escape_velocity(double r) const;

    /**
     * @brief Compute orbital period for circular orbit at radius r
     * @param r Orbital radius [m]
     * @return Orbital period [s]: T = 2π sqrt(r³/μ)
     */
    double orbital_period(double r) const;

private:
    double mu_;             ///< Gravitational parameter GM [m³/s²]
    Vec3 center_pos_;       ///< Position of attracting body [m]
    Vec3 center_vel_;       ///< Velocity of attracting body [m/s]
    bool center_is_moving_; ///< Whether center velocity is non-zero
};

} // namespace forces
} // namespace physim
