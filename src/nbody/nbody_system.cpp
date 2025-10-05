/**
 * @file nbody_system.cpp
 * @brief Implementation of N-body gravitational system
 */

#include "physim/nbody/nbody_system.hpp"
#include "physim/core/constants.hpp"
#include "physim/core/logging.hpp"
#include <stdexcept>
#include <cmath>

namespace physim {
namespace nbody {

NBodySystem::NBodySystem(bool use_softening, double softening_length)
    : current_time_(0.0)
    , initial_energy_(0.0)
    , energy_tracking_reset_(false)
    , use_softening_(use_softening)
    , softening_length_(softening_length) {

    if (use_softening && softening_length <= 0.0) {
        throw std::invalid_argument("Softening length must be positive");
    }
}

void NBodySystem::add_particle(const Particle& particle) {
    particles_.push_back(particle);
    energy_tracking_reset_ = false;  // Need to reset energy tracking
}

void NBodySystem::add_particle(std::shared_ptr<Particle> particle) {
    if (!particle) {
        throw std::invalid_argument("Cannot add null particle");
    }
    particles_.push_back(*particle);
    energy_tracking_reset_ = false;
}

void NBodySystem::remove_particle(size_t index) {
    if (index >= particles_.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    particles_.erase(particles_.begin() + index);
    energy_tracking_reset_ = false;
}

void NBodySystem::clear() {
    particles_.clear();
    current_time_ = 0.0;
    energy_tracking_reset_ = false;
}

const Particle& NBodySystem::particle(size_t index) const {
    if (index >= particles_.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    return particles_[index];
}

Particle& NBodySystem::particle(size_t index) {
    if (index >= particles_.size()) {
        throw std::out_of_range("Particle index out of range");
    }
    return particles_[index];
}

// ============================================================================
// Direct Summation Force Calculation
// ============================================================================

void NBodySystem::compute_accelerations() {
    const size_t n = particles_.size();
    const double G = constants::GRAVITATIONAL_CONSTANT;

    // Reset all accelerations to zero
    for (size_t i = 0; i < n; ++i) {
        particles_[i].set_acceleration(Vec3::Zero());
    }

    // O(N²) direct summation
    for (size_t i = 0; i < n; ++i) {
        Vec3 acc_i = Vec3::Zero();

        for (size_t j = 0; j < n; ++j) {
            if (i == j) continue;

            const Vec3 r_ij = particles_[j].position() - particles_[i].position();
            const double r2 = r_ij.squaredNorm();

            // Apply softening if enabled: F ∝ 1/(r² + ε²)^(3/2)
            const double r_soft2 = use_softening_ ?
                (r2 + softening_length_ * softening_length_) : r2;

            if (r_soft2 < 1e-20) {
                PHYSIM_LOG_WARN("Particle collision detected: r = {:.3e} m", std::sqrt(r2));
                continue;
            }

            const double r_soft3 = r_soft2 * std::sqrt(r_soft2);
            const double m_j = particles_[j].mass();

            // a_i += G * m_j * r_ij / r³
            acc_i += (G * m_j / r_soft3) * r_ij;
        }

        particles_[i].set_acceleration(acc_i);
    }
}

// ============================================================================
// Energy and Momentum Calculations
// ============================================================================

double NBodySystem::kinetic_energy() const {
    double T = 0.0;
    for (const auto& p : particles_) {
        T += p.kinetic_energy();
    }
    return T;
}

double NBodySystem::potential_energy() const {
    const size_t n = particles_.size();
    const double G = constants::GRAVITATIONAL_CONSTANT;
    double U = 0.0;

    // U = -Σ(i<j) G m_i m_j / r_ij
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            const double r_ij = particles_[i].distance_to(particles_[j]);

            if (r_ij < 1e-10) {
                PHYSIM_LOG_WARN("Near-collision in potential energy: r = {:.3e} m", r_ij);
                continue;
            }

            const double r_soft = use_softening_ ?
                std::sqrt(r_ij * r_ij + softening_length_ * softening_length_) : r_ij;

            U -= G * particles_[i].mass() * particles_[j].mass() / r_soft;
        }
    }

    return U;
}

double NBodySystem::total_energy() const {
    return kinetic_energy() + potential_energy();
}

Vec3 NBodySystem::center_of_mass() const {
    Vec3 r_cm = Vec3::Zero();
    double M = 0.0;

    for (const auto& p : particles_) {
        r_cm += p.mass() * p.position();
        M += p.mass();
    }

    return (M > 0.0) ? (r_cm / M) : Vec3::Zero();
}

Vec3 NBodySystem::center_of_mass_velocity() const {
    Vec3 v_cm = Vec3::Zero();
    double M = 0.0;

    for (const auto& p : particles_) {
        v_cm += p.mass() * p.velocity();
        M += p.mass();
    }

    return (M > 0.0) ? (v_cm / M) : Vec3::Zero();
}

Vec3 NBodySystem::total_momentum() const {
    Vec3 p_total = Vec3::Zero();
    for (const auto& p : particles_) {
        p_total += p.momentum();
    }
    return p_total;
}

Vec3 NBodySystem::total_angular_momentum() const {
    Vec3 L = Vec3::Zero();
    for (const auto& p : particles_) {
        L += p.angular_momentum();
    }
    return L;
}

Vec3 NBodySystem::total_angular_momentum(const Vec3& origin) const {
    Vec3 L = Vec3::Zero();
    for (const auto& p : particles_) {
        L += p.angular_momentum(origin);
    }
    return L;
}

double NBodySystem::total_mass() const {
    double M = 0.0;
    for (const auto& p : particles_) {
        M += p.mass();
    }
    return M;
}

NBodyStats NBodySystem::compute_stats() const {
    NBodyStats stats;

    stats.num_particles = particles_.size();
    stats.current_time = current_time_;
    stats.total_mass = total_mass();
    stats.kinetic_energy = kinetic_energy();
    stats.potential_energy = potential_energy();
    stats.total_energy = stats.kinetic_energy + stats.potential_energy;
    stats.center_of_mass = center_of_mass();
    stats.total_momentum = total_momentum();
    stats.total_angular_momentum = total_angular_momentum();

    // Compute energy error if tracking enabled
    if (energy_tracking_reset_) {
        const double energy_diff = std::abs(stats.total_energy - initial_energy_);
        stats.energy_error = (std::abs(initial_energy_) > 1e-10) ?
            (energy_diff / std::abs(initial_energy_)) : energy_diff;
    } else {
        stats.energy_error = 0.0;
    }

    return stats;
}

void NBodySystem::reset_energy_tracking() {
    initial_energy_ = total_energy();
    energy_tracking_reset_ = true;
    PHYSIM_LOG_DEBUG("Energy tracking reset: E0 = {:.6e} J", initial_energy_);
}

void NBodySystem::set_softening(bool enable, double epsilon) {
    if (enable && epsilon <= 0.0) {
        throw std::invalid_argument("Softening length must be positive");
    }
    use_softening_ = enable;
    softening_length_ = epsilon;
}

// ============================================================================
// State Vector Conversion
// ============================================================================

VecX NBodySystem::to_state_vector() const {
    const size_t n = particles_.size();
    VecX state(6 * n);

    for (size_t i = 0; i < n; ++i) {
        state.segment<3>(6 * i) = particles_[i].position();
        state.segment<3>(6 * i + 3) = particles_[i].velocity();
    }

    return state;
}

void NBodySystem::from_state_vector(const VecX& state_vec) {
    const size_t n = particles_.size();

    if (state_vec.size() != 6 * n) {
        throw std::invalid_argument("State vector size mismatch");
    }

    for (size_t i = 0; i < n; ++i) {
        particles_[i].set_position(state_vec.segment<3>(6 * i));
        particles_[i].set_velocity(state_vec.segment<3>(6 * i + 3));
    }
}

// ============================================================================
// ODE Integration
// ============================================================================

integrators::DerivativeFunction NBodySystem::derivative_function() {
    // Capture 'this' to access compute_accelerations()
    return [this](double t, const VecX& y) -> VecX {
        (void)t;  // Time not explicitly used

        // Update particle states from y
        this->from_state_vector(y);

        // Compute accelerations via direct summation
        this->compute_accelerations();

        const size_t n = this->particles_.size();
        VecX dydt(6 * n);

        // dy/dt = [velocities, accelerations]
        for (size_t i = 0; i < n; ++i) {
            dydt.segment<3>(6 * i) = this->particles_[i].velocity();
            dydt.segment<3>(6 * i + 3) = this->particles_[i].acceleration();
        }

        return dydt;
    };
}

void NBodySystem::step(double dt,
                       integrators::Integrator& integrator,
                       bool adaptive,
                       double tolerance) {
    if (particles_.empty()) {
        PHYSIM_LOG_WARN("Cannot step empty N-body system");
        return;
    }

    // Get current state
    VecX y0 = to_state_vector();

    // Create derivative function
    auto f = derivative_function();

    // Integrate one step
    VecX yf;
    if (adaptive) {
        yf = integrator.integrate(current_time_, y0, current_time_ + dt,
                                  dt, f, true, tolerance);
    } else {
        yf = integrator.step(current_time_, y0, dt, f);
    }

    // Update system state
    from_state_vector(yf);
    current_time_ += dt;

    // Recompute accelerations for consistency
    compute_accelerations();
}

void NBodySystem::propagate(double t_final,
                            double dt,
                            integrators::Integrator& integrator,
                            bool adaptive,
                            double tolerance) {
    if (particles_.empty()) {
        PHYSIM_LOG_WARN("Cannot propagate empty N-body system");
        return;
    }

    if (t_final <= current_time_) {
        throw std::invalid_argument("Final time must be greater than current time");
    }

    // Reset energy tracking at start of propagation
    if (!energy_tracking_reset_) {
        reset_energy_tracking();
    }

    PHYSIM_LOG_INFO("Propagating N-body system: {} particles, t={:.3f}s -> {:.3f}s",
                    particles_.size(), current_time_, t_final);

    // Get initial state
    VecX y0 = to_state_vector();

    // Create derivative function
    auto f = derivative_function();

    // Integrate over full time interval
    VecX yf = integrator.integrate(current_time_, y0, t_final, dt, f,
                                    adaptive, tolerance);

    // Update system state
    from_state_vector(yf);
    current_time_ = t_final;

    // Recompute accelerations
    compute_accelerations();

    PHYSIM_LOG_INFO("Propagation complete: t={:.3f}s", current_time_);
}

} // namespace nbody
} // namespace physim
