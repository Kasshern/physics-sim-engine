/**
 * @file gravity.cpp
 * @brief Implementation of point mass gravity
 */

#include "physim/forces/gravity.hpp"
#include "physim/core/logging.hpp"
#include <cmath>
#include <stdexcept>
#include <algorithm>

namespace physim {
namespace forces {

PointMassGravity::PointMassGravity(double GM,
                                   const Vec3& center_pos,
                                   const Vec3& center_vel)
    : Force("PointMassGravity"),
      mu_(GM),
      center_pos_(center_pos),
      center_vel_(center_vel),
      center_is_moving_(center_vel.norm() > 1e-10) {
    if (mu_ <= 0.0) {
        throw std::invalid_argument("Gravitational parameter GM must be positive");
    }

    PHYSIM_LOG_DEBUG("Created point mass gravity: μ={:.6e} m³/s²", mu_);
}

std::shared_ptr<PointMassGravity> PointMassGravity::for_body(
    const std::string& body_name,
    const Vec3& center_pos,
    const Vec3& center_vel) {

    // Convert body name to lowercase for case-insensitive matching
    std::string name_lower = body_name;
    std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(),
                  [](unsigned char c) { return std::tolower(c); });

    double mu = 0.0;

    // Match against known solar system bodies
    if (name_lower == "sun") {
        mu = constants::gm::SUN;
    } else if (name_lower == "mercury") {
        mu = constants::gm::MERCURY;
    } else if (name_lower == "venus") {
        mu = constants::gm::VENUS;
    } else if (name_lower == "earth") {
        mu = constants::gm::EARTH;
    } else if (name_lower == "moon") {
        mu = constants::gm::MOON;
    } else if (name_lower == "mars") {
        mu = constants::gm::MARS;
    } else if (name_lower == "jupiter") {
        mu = constants::gm::JUPITER;
    } else if (name_lower == "saturn") {
        mu = constants::gm::SATURN;
    } else if (name_lower == "uranus") {
        mu = constants::gm::URANUS;
    } else if (name_lower == "neptune") {
        mu = constants::gm::NEPTUNE;
    } else if (name_lower == "pluto") {
        mu = constants::gm::PLUTO;
    } else {
        throw std::invalid_argument("Unknown body name: " + body_name);
    }

    auto gravity = std::make_shared<PointMassGravity>(mu, center_pos, center_vel);
    gravity->name_ = "Gravity_" + body_name;

    PHYSIM_LOG_INFO("Created gravity force for {}: μ={:.6e} m³/s²",
                   body_name, mu);

    return gravity;
}

Vec3 PointMassGravity::acceleration(double t, const Vec3& pos, const Vec3& vel,
                                    double mass) const {
    (void)t;    // Time-independent for stationary center
    (void)vel;  // Velocity-independent
    (void)mass; // Mass cancels in a = F/m

    // Relative position vector from center to particle
    const Vec3 r_vec = pos - center_pos_;
    const double r = r_vec.norm();

    // Handle singularity at r = 0
    if (r < 1e-10) {
        PHYSIM_LOG_ERROR("Gravitational singularity: r = {:.3e} m (too close to center)", r);
        throw std::runtime_error("Gravitational singularity: particle too close to center");
    }

    // Newton's law: a = -μ/r² * r̂
    // r̂ = r_vec/r, so: a = -μ/r³ * r_vec
    const double r3 = r * r * r;
    return -mu_ / r3 * r_vec;
}

double PointMassGravity::potential_energy(double t, const Vec3& pos, double mass) const {
    (void)t;  // Time-independent

    const Vec3 r_vec = pos - center_pos_;
    const double r = r_vec.norm();

    if (r < 1e-10) {
        PHYSIM_LOG_WARN("Potential energy singular at r={:.3e} m", r);
        return -std::numeric_limits<double>::infinity();
    }

    // Gravitational potential: U = -GMm/r
    return -mu_ * mass / r;
}

double PointMassGravity::circular_velocity(double r) const {
    if (r <= 0.0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }

    // Circular orbit velocity: v = sqrt(μ/r)
    return std::sqrt(mu_ / r);
}

double PointMassGravity::escape_velocity(double r) const {
    if (r <= 0.0) {
        throw std::invalid_argument("Radius must be positive");
    }

    // Escape velocity: v = sqrt(2μ/r)
    return std::sqrt(2.0 * mu_ / r);
}

double PointMassGravity::orbital_period(double r) const {
    if (r <= 0.0) {
        throw std::invalid_argument("Orbital radius must be positive");
    }

    // Kepler's third law: T = 2π sqrt(r³/μ)
    return constants::TWO_PI * std::sqrt(r * r * r / mu_);
}

} // namespace forces
} // namespace physim
