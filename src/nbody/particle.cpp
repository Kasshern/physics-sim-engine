/**
 * @file particle.cpp
 * @brief Implementation of Particle class for N-body simulations
 */

#include "physim/nbody/particle.hpp"
#include "physim/core/constants.hpp"
#include "physim/core/logging.hpp"
#include <stdexcept>
#include <cmath>

namespace physim {
namespace nbody {

Particle::Particle()
    : position_(Vec3::Zero())
    , velocity_(Vec3::Zero())
    , acceleration_(Vec3::Zero())
    , mass_(0.0)
    , name_("") {
}

Particle::Particle(const Vec3& position, const Vec3& velocity, double mass,
                   const std::string& name)
    : position_(position)
    , velocity_(velocity)
    , acceleration_(Vec3::Zero())
    , mass_(mass)
    , name_(name) {

    if (mass <= 0.0) {
        throw std::invalid_argument("Particle mass must be positive");
    }
}

Particle::Particle(const State& state, const std::string& name)
    : position_(state.position)
    , velocity_(state.velocity)
    , acceleration_(Vec3::Zero())
    , mass_(state.mass)
    , name_(name) {

    if (mass_ <= 0.0) {
        throw std::invalid_argument("Particle mass must be positive");
    }
}

void Particle::set_mass(double m) {
    if (m <= 0.0) {
        throw std::invalid_argument("Particle mass must be positive");
    }
    mass_ = m;
}

double Particle::kinetic_energy() const {
    return 0.5 * mass_ * velocity_.squaredNorm();
}

Vec3 Particle::momentum() const {
    return mass_ * velocity_;
}

Vec3 Particle::angular_momentum() const {
    return position_.cross(mass_ * velocity_);
}

Vec3 Particle::angular_momentum(const Vec3& origin) const {
    const Vec3 r_rel = position_ - origin;
    return r_rel.cross(mass_ * velocity_);
}

double Particle::distance_to(const Particle& other) const {
    return (position_ - other.position_).norm();
}

double Particle::distance_squared_to(const Particle& other) const {
    return (position_ - other.position_).squaredNorm();
}

State Particle::to_state() const {
    State state;
    state.position = position_;
    state.velocity = velocity_;
    state.mass = mass_;
    state.time = 0.0;  // Time is managed by NBodySystem
    return state;
}

void Particle::from_state(const State& state) {
    position_ = state.position;
    velocity_ = state.velocity;

    if (state.mass > 0.0) {
        mass_ = state.mass;
    }
}

VecX Particle::to_vector() const {
    VecX vec(6);
    vec << position_, velocity_;
    return vec;
}

void Particle::from_vector(const VecX& state_vec) {
    if (state_vec.size() != 6) {
        throw std::invalid_argument("State vector must have 6 elements");
    }

    position_ = state_vec.segment<3>(0);
    velocity_ = state_vec.segment<3>(3);
}

// ============================================================================
// Solar System Particle Factory
// ============================================================================

std::shared_ptr<Particle> create_solar_system_particle(
    const std::string& body_name,
    const Vec3& position,
    const Vec3& velocity) {

    double mass = 0.0;

    // Get mass from GM values
    if (body_name == "Sun") {
        mass = constants::gm::SUN / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Mercury") {
        mass = constants::gm::MERCURY / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Venus") {
        mass = constants::gm::VENUS / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Earth") {
        mass = constants::gm::EARTH / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Moon") {
        mass = constants::gm::MOON / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Mars") {
        mass = constants::gm::MARS / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Jupiter") {
        mass = constants::gm::JUPITER / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Saturn") {
        mass = constants::gm::SATURN / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Uranus") {
        mass = constants::gm::URANUS / constants::GRAVITATIONAL_CONSTANT;
    } else if (body_name == "Neptune") {
        mass = constants::gm::NEPTUNE / constants::GRAVITATIONAL_CONSTANT;
    } else {
        PHYSIM_LOG_ERROR("Unknown solar system body: {}", body_name);
        throw std::invalid_argument("Unknown solar system body: " + body_name);
    }

    auto particle = std::make_shared<Particle>(position, velocity, mass, body_name);

    PHYSIM_LOG_INFO("Created particle '{}': mass={:.3e} kg", body_name, mass);

    return particle;
}

} // namespace nbody
} // namespace physim
