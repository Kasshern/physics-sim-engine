/**
 * @file j2_gravity.cpp
 * @brief Implementation of J2 gravitational perturbation
 */

#include "physim/forces/j2_gravity.hpp"
#include "physim/core/logging.hpp"
#include <cmath>
#include <stdexcept>

namespace physim {
namespace forces {

J2Gravity::J2Gravity(double mu, double R, double J2, const Vec3& rotation_axis)
    : Force("J2Gravity"),
      mu_(mu),
      R_(R),
      J2_(J2),
      rotation_axis_(rotation_axis.normalized()) {
    if (mu_ <= 0.0) {
        throw std::invalid_argument("Gravitational parameter must be positive");
    }
    if (R_ <= 0.0) {
        throw std::invalid_argument("Reference radius must be positive");
    }

    PHYSIM_LOG_DEBUG("Created J2 gravity: μ={:.6e} m³/s², R={:.3e} m, J2={:.6e}",
                    mu_, R_, J2_);
}

Vec3 J2Gravity::acceleration(double t, const Vec3& pos, const Vec3& vel,
                             double mass) const {
    (void)t;    // Time-independent
    (void)vel;  // Velocity-independent
    (void)mass; // Mass cancels

    const double r = pos.norm();

    if (r < 1e-10) {
        PHYSIM_LOG_ERROR("J2 singularity: r = {:.3e} m", r);
        throw std::runtime_error("J2 gravity singular at origin");
    }

    // Component along rotation axis
    const double z = pos.dot(rotation_axis_);

    // Precompute powers
    const double r2 = r * r;
    const double r5 = r2 * r2 * r;
    const double z2 = z * z;
    const double R2 = R_ * R_;

    // Common factor: (3/2) * J2 * μ * R²/r⁵
    const double factor = 1.5 * J2_ * mu_ * R2 / r5;

    // J2 acceleration components:
    //   a_⊥ = factor * (5z²/r² - 1) * r_⊥
    //   a_z = factor * (5z²/r² - 3) * z
    //
    // Split position into perpendicular and parallel components:
    //   r = r_⊥ + z * ẑ
    //   r_⊥ = r - z * ẑ

    const double coeff_perp = 5.0 * z2 / r2 - 1.0;
    const double coeff_z = 5.0 * z2 / r2 - 3.0;

    // Perpendicular component
    const Vec3 r_perp = pos - z * rotation_axis_;

    // Total acceleration
    return factor * (coeff_perp * r_perp + coeff_z * z * rotation_axis_);
}

double J2Gravity::potential_energy(double t, const Vec3& pos, double mass) const {
    (void)t;  // Time-independent

    const double r = pos.norm();

    if (r < 1e-10) {
        PHYSIM_LOG_WARN("J2 potential singular at r={:.3e} m", r);
        return -std::numeric_limits<double>::infinity();
    }

    // Component along rotation axis
    const double z = pos.dot(rotation_axis_);

    // Latitude: sin(φ) = z/r
    const double sin_phi = z / r;
    const double sin2_phi = sin_phi * sin_phi;

    // J2 potential: U = -(μm/r) * (R/r)² * (J2/2) * (3sin²φ - 1)
    const double R_over_r = R_ / r;
    const double R_over_r2 = R_over_r * R_over_r;

    return -(mu_ * mass / r) * R_over_r2 * (J2_ / 2.0) * (3.0 * sin2_phi - 1.0);
}

double J2Gravity::nodal_regression_rate(double a, double i) const {
    if (a <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }

    // Mean motion: n = sqrt(μ/a³)
    const double n = std::sqrt(mu_ / (a * a * a));

    // Nodal regression: Ω̇ = -(3/2) * n * J2 * (R/a)² * cos(i)
    const double R_over_a = R_ / a;
    const double R_over_a2 = R_over_a * R_over_a;

    return -1.5 * n * J2_ * R_over_a2 * std::cos(i);
}

double J2Gravity::apsidal_precession_rate(double a, double e, double i) const {
    if (a <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }
    if (e < 0.0 || e >= 1.0) {
        throw std::invalid_argument("Eccentricity must be in [0, 1)");
    }

    // Mean motion: n = sqrt(μ/a³)
    const double n = std::sqrt(mu_ / (a * a * a));

    // Apsidal precession: ω̇ = -(3/4) * n * J2 * (R/a)² * (5cos²i - 1)
    const double R_over_a = R_ / a;
    const double R_over_a2 = R_over_a * R_over_a;
    const double cos_i = std::cos(i);
    const double cos2_i = cos_i * cos_i;

    return -0.75 * n * J2_ * R_over_a2 * (5.0 * cos2_i - 1.0);
}

double J2Gravity::sun_synchronous_inclination(double a) const {
    if (a <= 0.0) {
        throw std::invalid_argument("Semi-major axis must be positive");
    }

    // Sun-synchronous requirement: Ω̇ = 360°/year = 2π/(365.25 days)
    const double omega_dot_ss = constants::TWO_PI / constants::TROPICAL_YEAR;

    // From nodal regression formula:
    //   Ω̇ = -(3/2) * n * J2 * (R/a)² * cos(i)
    //
    // Solve for cos(i):
    //   cos(i) = Ω̇ / [-(3/2) * n * J2 * (R/a)²]

    const double n = std::sqrt(mu_ / (a * a * a));
    const double R_over_a = R_ / a;
    const double R_over_a2 = R_over_a * R_over_a;

    const double cos_i = omega_dot_ss / (-1.5 * n * J2_ * R_over_a2);

    // Check if physically realizable
    if (std::abs(cos_i) > 1.0) {
        PHYSIM_LOG_ERROR("Sun-synchronous orbit not possible at a={:.3e} m: "
                        "requires |cos(i)|={:.3f} > 1",
                        a, std::abs(cos_i));
        throw std::runtime_error("Sun-synchronous orbit not realizable at this altitude");
    }

    // Return inclination (retrograde orbit: i > 90°)
    return std::acos(cos_i);
}

} // namespace forces
} // namespace physim
