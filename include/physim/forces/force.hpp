/**
 * @file force.hpp
 * @brief Base class for force models in orbital mechanics
 *
 * Provides abstract interface for computing forces and accelerations
 * in N-body gravitational systems and spacecraft dynamics.
 */

#pragma once

#include "physim/core/types.hpp"
#include "physim/core/time.hpp"
#include <memory>
#include <string>
#include <vector>

namespace physim {
namespace forces {

// ============================================================================
// Force Base Class
// ============================================================================

/**
 * @brief Abstract base class for force models
 *
 * Represents a physical force acting on bodies in the simulation.
 * Forces compute acceleration a = F/m for use in equations of motion.
 *
 * Examples:
 * - Point mass gravity: F = -GMm/r² * r̂
 * - J2 perturbation: Earth oblateness effects
 * - Solar radiation pressure
 * - Atmospheric drag
 * - Thrust from engines
 */
class Force {
public:
    /**
     * @brief Constructor
     * @param name Human-readable name of the force
     */
    explicit Force(std::string name) : name_(std::move(name)), enabled_(true) {}

    /**
     * @brief Virtual destructor
     */
    virtual ~Force() = default;

    /**
     * @brief Get force name
     */
    const std::string& name() const { return name_; }

    /**
     * @brief Check if force is enabled
     */
    bool is_enabled() const { return enabled_; }

    /**
     * @brief Enable or disable this force
     * @param enabled True to enable, false to disable
     */
    void set_enabled(bool enabled) { enabled_ = enabled; }

    /**
     * @brief Compute acceleration due to this force
     *
     * @param t Current time [s]
     * @param pos Position vector [m] in inertial frame
     * @param vel Velocity vector [m/s] in inertial frame
     * @param mass Mass of the body [kg]
     * @return Acceleration [m/s²]
     *
     * This is the primary interface for force computation.
     * Returns a = F/m where F is the force on the body.
     */
    virtual Vec3 acceleration(double t, const Vec3& pos, const Vec3& vel,
                             double mass) const = 0;

    /**
     * @brief Compute force (optional override)
     *
     * @param t Current time [s]
     * @param pos Position vector [m]
     * @param vel Velocity vector [m/s]
     * @param mass Mass of the body [kg]
     * @return Force [N]
     *
     * Default implementation: F = m * a
     * Some force models may override this for efficiency.
     */
    virtual Vec3 force(double t, const Vec3& pos, const Vec3& vel,
                       double mass) const {
        return mass * acceleration(t, pos, vel, mass);
    }

    /**
     * @brief Check if force depends on velocity
     * @return true if force is velocity-dependent (drag, Coriolis, etc.)
     *
     * Allows optimizations when forces are position-only.
     */
    virtual bool is_velocity_dependent() const { return false; }

    /**
     * @brief Check if force is time-dependent
     * @return true if force changes with time (thrust profile, etc.)
     *
     * Conservative forces return false.
     */
    virtual bool is_time_dependent() const { return false; }

    /**
     * @brief Compute potential energy (for conservative forces)
     *
     * @param t Current time [s]
     * @param pos Position vector [m]
     * @param mass Mass of the body [kg]
     * @return Potential energy [J]
     *
     * For non-conservative forces, returns 0.0.
     * Conservative forces should override this.
     */
    virtual double potential_energy(double t, const Vec3& pos, double mass) const {
        (void)t;
        (void)pos;
        (void)mass;
        return 0.0;
    }

protected:
    std::string name_;  ///< Force name
    bool enabled_;      ///< Whether this force is active
};

// ============================================================================
// Force Container
// ============================================================================

/**
 * @brief Container for multiple force models
 *
 * Manages a collection of forces acting on a body.
 * Computes total acceleration as sum of all enabled forces.
 */
class ForceModel {
public:
    /**
     * @brief Default constructor
     */
    ForceModel() = default;

    /**
     * @brief Add a force to the model
     * @param force Shared pointer to force
     */
    void add_force(std::shared_ptr<Force> force);

    /**
     * @brief Remove all forces
     */
    void clear() { forces_.clear(); }

    /**
     * @brief Get number of forces
     */
    size_t num_forces() const { return forces_.size(); }

    /**
     * @brief Get force by index
     * @param index Force index
     * @return Shared pointer to force
     */
    std::shared_ptr<Force> get_force(size_t index) const;

    /**
     * @brief Compute total acceleration from all enabled forces
     *
     * @param t Current time [s]
     * @param pos Position vector [m]
     * @param vel Velocity vector [m/s]
     * @param mass Mass of the body [kg]
     * @return Total acceleration [m/s²]
     *
     * Sums acceleration from all enabled forces:
     * a_total = Σ a_i
     */
    Vec3 total_acceleration(double t, const Vec3& pos, const Vec3& vel,
                           double mass) const;

    /**
     * @brief Compute total potential energy (conservative forces only)
     *
     * @param t Current time [s]
     * @param pos Position vector [m]
     * @param mass Mass of the body [kg]
     * @return Total potential energy [J]
     */
    double total_potential_energy(double t, const Vec3& pos, double mass) const;

    /**
     * @brief Enable/disable force by name
     * @param name Force name
     * @param enabled True to enable, false to disable
     * @return true if force found, false otherwise
     */
    bool set_force_enabled(const std::string& name, bool enabled);

    /**
     * @brief Get list of all force names
     */
    std::vector<std::string> get_force_names() const;

private:
    std::vector<std::shared_ptr<Force>> forces_;  ///< List of forces
};

} // namespace forces
} // namespace physim
