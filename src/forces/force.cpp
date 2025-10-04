/**
 * @file force.cpp
 * @brief Implementation of force base class and container
 */

#include "physim/forces/force.hpp"
#include "physim/core/logging.hpp"
#include <stdexcept>

namespace physim {
namespace forces {

// ============================================================================
// ForceModel Implementation
// ============================================================================

void ForceModel::add_force(std::shared_ptr<Force> force) {
    if (!force) {
        throw std::invalid_argument("Cannot add null force to model");
    }
    forces_.push_back(force);
    PHYSIM_LOG_DEBUG("Added force '{}' to model ({} total forces)",
                    force->name(), forces_.size());
}

std::shared_ptr<Force> ForceModel::get_force(size_t index) const {
    if (index >= forces_.size()) {
        throw std::out_of_range("Force index out of range");
    }
    return forces_[index];
}

Vec3 ForceModel::total_acceleration(double t, const Vec3& pos, const Vec3& vel,
                                   double mass) const {
    Vec3 total_accel = Vec3::Zero();

    for (const auto& force : forces_) {
        if (force->is_enabled()) {
            total_accel += force->acceleration(t, pos, vel, mass);
        }
    }

    return total_accel;
}

double ForceModel::total_potential_energy(double t, const Vec3& pos, double mass) const {
    double total_energy = 0.0;

    for (const auto& force : forces_) {
        if (force->is_enabled()) {
            total_energy += force->potential_energy(t, pos, mass);
        }
    }

    return total_energy;
}

bool ForceModel::set_force_enabled(const std::string& name, bool enabled) {
    for (auto& force : forces_) {
        if (force->name() == name) {
            force->set_enabled(enabled);
            PHYSIM_LOG_INFO("Force '{}' {}", name, enabled ? "enabled" : "disabled");
            return true;
        }
    }

    PHYSIM_LOG_WARN("Force '{}' not found in model", name);
    return false;
}

std::vector<std::string> ForceModel::get_force_names() const {
    std::vector<std::string> names;
    names.reserve(forces_.size());

    for (const auto& force : forces_) {
        names.push_back(force->name());
    }

    return names;
}

} // namespace forces
} // namespace physim
