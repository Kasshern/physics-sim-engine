/**
 * @file j2_gravity.hpp
 * @brief J2 gravitational perturbation (Earth oblateness)
 *
 * Implements the dominant aspherical gravity effect due to Earth's
 * equatorial bulge. Essential for accurate LEO orbit propagation.
 */

#pragma once

#include "physim/forces/force.hpp"
#include "physim/core/constants.hpp"

namespace physim {
namespace forces {

/**
 * @brief J2 gravitational perturbation force
 *
 * Models the largest aspherical component of Earth's gravity field
 * due to its oblate shape (equatorial bulge). J2 is the second-degree
 * zonal harmonic coefficient.
 *
 * Effects:
 * - Nodal regression: Ω̇ (orbit plane precesses)
 * - Apsidal precession: ω̇ (argument of periapsis rotates)
 * - Changes in inclination: small but measurable
 *
 * The J2 perturbation acceleration is:
 *   a_J2 = (3/2) * J2 * μ * R²/r⁵ * [
 *     (5z²/r² - 1) * x̂ +
 *     (5z²/r² - 1) * ŷ +
 *     (5z²/r² - 3) * ẑ
 *   ]
 *
 * where:
 * - J2 = 1.08263×10⁻³ (Earth's J2 coefficient)
 * - μ = GM (gravitational parameter)
 * - R = equatorial radius
 * - r = distance from center
 * - z = component along rotation axis
 *
 * Magnitude:
 * - ~1000× smaller than point mass gravity
 * - Dominant perturbation for LEO satellites
 * - Causes ~5°/day nodal regression for ISS
 *
 * Properties:
 * - Conservative (has potential energy)
 * - Velocity-independent
 * - Time-independent (Earth rotation axis fixed in inertial frame)
 * - Axially symmetric (no longitude dependence)
 *
 * Use cases:
 * - LEO satellite orbit propagation
 * - Sun-synchronous orbit design
 * - Nodal regression calculations
 * - GPS orbit maintenance
 *
 * Limitations:
 * - Only includes J2 (ignores J3, J4, ..., tesserals)
 * - Assumes fixed rotation axis (no precession/nutation)
 * - For high-precision, use full EGM2008 gravity model
 */
class J2Gravity : public Force {
public:
    /**
     * @brief Constructor with default Earth parameters
     * @param mu Gravitational parameter GM [m³/s²] (default: Earth)
     * @param R Reference radius [m] (default: Earth equatorial radius)
     * @param J2 J2 coefficient (dimensionless, default: Earth J2)
     * @param rotation_axis Unit vector along rotation axis (default: z-axis)
     */
    explicit J2Gravity(double mu = constants::gm::EARTH,
                      double R = constants::earth::WGS84_A,
                      double J2 = constants::earth::J2,
                      const Vec3& rotation_axis = Vec3::UnitZ());

    /**
     * @brief Compute J2 perturbation acceleration
     *
     * @param t Current time [s] (unused)
     * @param pos Position vector [m] in inertial frame
     * @param vel Velocity vector [m/s] (unused)
     * @param mass Mass [kg] (unused)
     * @return Acceleration [m/s²]
     *
     * Computes J2 perturbation using spherical harmonics formula.
     */
    Vec3 acceleration(double t, const Vec3& pos, const Vec3& vel,
                     double mass) const override;

    /**
     * @brief Compute J2 potential energy
     *
     * @param t Current time [s] (unused)
     * @param pos Position vector [m]
     * @param mass Mass [kg]
     * @return Potential energy [J]
     *
     * U_J2 = -(μm/r) * (R/r)² * (J2/2) * (3sin²φ - 1)
     * where φ is latitude
     */
    double potential_energy(double t, const Vec3& pos, double mass) const override;

    /**
     * @brief Check if velocity-dependent (always false)
     */
    bool is_velocity_dependent() const override { return false; }

    /**
     * @brief Check if time-dependent (always false)
     */
    bool is_time_dependent() const override { return false; }

    /**
     * @brief Get J2 coefficient
     */
    double j2() const { return J2_; }

    /**
     * @brief Get reference radius
     */
    double reference_radius() const { return R_; }

    /**
     * @brief Get gravitational parameter
     */
    double mu() const { return mu_; }

    /**
     * @brief Compute nodal regression rate for circular orbit
     *
     * @param a Semi-major axis [m]
     * @param i Inclination [radians]
     * @return Nodal regression rate Ω̇ [rad/s]
     *
     * Formula: Ω̇ = -(3/2) * n * J2 * (R/a)² * cos(i)
     * where n = sqrt(μ/a³) is mean motion
     */
    double nodal_regression_rate(double a, double i) const;

    /**
     * @brief Compute apsidal precession rate for circular orbit
     *
     * @param a Semi-major axis [m]
     * @param e Eccentricity
     * @param i Inclination [radians]
     * @return Apsidal precession rate ω̇ [rad/s]
     *
     * Formula: ω̇ = -(3/4) * n * J2 * (R/a)² * (5cos²i - 1)
     */
    double apsidal_precession_rate(double a, double e, double i) const;

    /**
     * @brief Compute inclination for sun-synchronous orbit
     *
     * @param a Semi-major axis [m]
     * @return Inclination [radians] for sun-synchronous precession
     *
     * Sun-synchronous orbits have Ω̇ = 360°/year to maintain constant
     * solar angle. Requires retrograde orbit (i > 90°).
     */
    double sun_synchronous_inclination(double a) const;

private:
    double mu_;             ///< Gravitational parameter GM [m³/s²]
    double R_;              ///< Reference radius [m]
    double J2_;             ///< J2 coefficient (dimensionless)
    Vec3 rotation_axis_;    ///< Rotation axis (unit vector)
};

} // namespace forces
} // namespace physim
