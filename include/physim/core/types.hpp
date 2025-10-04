/**
 * @file types.hpp
 * @brief Core type definitions for the physics simulation engine
 *
 * This file defines fundamental types used throughout the simulation,
 * including vectors, quaternions, matrices, and state representations.
 * All types are built on Eigen for high-performance SIMD operations.
 */

#pragma once

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <cstdint>
#include <string>
#include <vector>
#include <array>

namespace physim {

// ============================================================================
// Fundamental Types
// ============================================================================

/// 3D vector using double precision (position, velocity, force, etc.)
using Vec3 = Eigen::Vector3d;

/// 3x3 matrix using double precision (rotation, inertia tensor, etc.)
using Mat3 = Eigen::Matrix3d;

/// Unit quaternion for rotation representation
using Quat = Eigen::Quaterniond;

/// 6D vector for position+velocity or other 6-DOF quantities
using Vec6 = Eigen::Matrix<double, 6, 1>;

/// Dynamic-size vector
using VecX = Eigen::VectorXd;

/// Dynamic-size matrix
using MatX = Eigen::MatrixXd;

// ============================================================================
// Physical State Representation
// ============================================================================

/**
 * @brief Complete dynamic state of a rigid body
 *
 * Represents the full state vector for a rigid body in 6-DOF:
 * - Translational: position, velocity
 * - Rotational: orientation (quaternion), angular velocity
 * - Properties: mass, inertia tensor
 *
 * All quantities use SI units unless otherwise specified.
 */
struct State {
    Vec3 position;           ///< Position vector [m] in inertial frame
    Vec3 velocity;           ///< Velocity vector [m/s] in inertial frame
    Quat orientation;        ///< Orientation quaternion (unit) from body to inertial
    Vec3 angular_velocity;   ///< Angular velocity [rad/s] in body frame
    double mass;             ///< Mass [kg]
    Mat3 inertia;            ///< Inertia tensor [kg·m²] in body frame
    double time;             ///< Simulation time [s] since epoch

    /**
     * @brief Default constructor initializes to zero state
     */
    State()
        : position(Vec3::Zero()),
          velocity(Vec3::Zero()),
          orientation(Quat::Identity()),
          angular_velocity(Vec3::Zero()),
          mass(1.0),
          inertia(Mat3::Identity()),
          time(0.0) {}

    /**
     * @brief Construct state with position and velocity only
     * @param pos Position [m]
     * @param vel Velocity [m/s]
     * @param m Mass [kg]
     */
    State(const Vec3& pos, const Vec3& vel, double m = 1.0)
        : position(pos),
          velocity(vel),
          orientation(Quat::Identity()),
          angular_velocity(Vec3::Zero()),
          mass(m),
          inertia(Mat3::Identity()),
          time(0.0) {}

    /**
     * @brief Construct full 6-DOF state
     * @param pos Position [m]
     * @param vel Velocity [m/s]
     * @param quat Orientation quaternion (will be normalized)
     * @param angvel Angular velocity [rad/s]
     * @param m Mass [kg]
     * @param I Inertia tensor [kg·m²]
     * @param t Time [s]
     */
    State(const Vec3& pos, const Vec3& vel, const Quat& quat,
          const Vec3& angvel, double m, const Mat3& I, double t = 0.0)
        : position(pos),
          velocity(vel),
          orientation(quat.normalized()),
          angular_velocity(angvel),
          mass(m),
          inertia(I),
          time(t) {}
};

// ============================================================================
// Derivative State (for ODE integration)
// ============================================================================

/**
 * @brief Time derivative of State for integration
 *
 * Represents dState/dt, used by integrators to advance the state.
 * Contains velocity, acceleration, and rotational derivatives.
 */
struct StateDerivative {
    Vec3 velocity;                ///< dr/dt = v [m/s]
    Vec3 acceleration;            ///< dv/dt = a [m/s²]
    Quat orientation_derivative;  ///< dq/dt from angular velocity
    Vec3 angular_acceleration;    ///< dω/dt [rad/s²]

    StateDerivative()
        : velocity(Vec3::Zero()),
          acceleration(Vec3::Zero()),
          orientation_derivative(Quat::Identity()),
          angular_acceleration(Vec3::Zero()) {}
};

// ============================================================================
// Particle (Simple Point Mass)
// ============================================================================

/**
 * @brief Simplified representation for N-body simulations
 *
 * A particle is a point mass with position and velocity but no orientation.
 * Used for large-scale gravitational simulations where attitude is not needed.
 */
struct Particle {
    Vec3 position;     ///< Position [m]
    Vec3 velocity;     ///< Velocity [m/s]
    Vec3 acceleration; ///< Acceleration [m/s²] (computed each step)
    double mass;       ///< Mass [kg]
    uint64_t id;       ///< Unique identifier
    std::string name;  ///< Human-readable name

    Particle()
        : position(Vec3::Zero()),
          velocity(Vec3::Zero()),
          acceleration(Vec3::Zero()),
          mass(1.0),
          id(0),
          name("") {}

    Particle(const Vec3& pos, const Vec3& vel, double m,
             uint64_t particle_id = 0, const std::string& particle_name = "")
        : position(pos),
          velocity(vel),
          acceleration(Vec3::Zero()),
          mass(m),
          id(particle_id),
          name(particle_name) {}
};

// ============================================================================
// Physical Properties
// ============================================================================

/**
 * @brief Physical properties of a celestial body or spacecraft
 */
struct BodyProperties {
    double mass;              ///< Mass [kg]
    double radius;            ///< Mean radius [m]
    double gravitational_parameter; ///< GM [m³/s²]
    Mat3 inertia;             ///< Inertia tensor [kg·m²]
    double rotation_period;   ///< Rotation period [s]
    Vec3 rotation_axis;       ///< Rotation axis (unit vector)

    // Spacecraft-specific properties
    double drag_coefficient;  ///< Drag coefficient Cd (dimensionless)
    double cross_section;     ///< Cross-sectional area [m²]
    double reflectivity;      ///< Surface reflectivity [0,1]

    BodyProperties()
        : mass(1.0),
          radius(1.0),
          gravitational_parameter(1.0),
          inertia(Mat3::Identity()),
          rotation_period(86400.0),
          rotation_axis(Vec3::UnitZ()),
          drag_coefficient(2.2),
          cross_section(1.0),
          reflectivity(0.3) {}
};

// ============================================================================
// Bounding Volumes (for collision detection)
// ============================================================================

/**
 * @brief Axis-aligned bounding box
 */
struct AABB {
    Vec3 min; ///< Minimum corner
    Vec3 max; ///< Maximum corner

    AABB() : min(Vec3::Zero()), max(Vec3::Zero()) {}

    AABB(const Vec3& minimum, const Vec3& maximum)
        : min(minimum), max(maximum) {}

    /**
     * @brief Get center of bounding box
     */
    Vec3 center() const { return 0.5 * (min + max); }

    /**
     * @brief Get half-extents
     */
    Vec3 extents() const { return 0.5 * (max - min); }

    /**
     * @brief Check if point is inside AABB
     */
    bool contains(const Vec3& point) const {
        return (point.array() >= min.array()).all() &&
               (point.array() <= max.array()).all();
    }

    /**
     * @brief Check if two AABBs intersect
     */
    bool intersects(const AABB& other) const {
        return (min.array() <= other.max.array()).all() &&
               (max.array() >= other.min.array()).all();
    }
};

/**
 * @brief Sphere bounding volume
 */
struct Sphere {
    Vec3 center;  ///< Center position [m]
    double radius; ///< Radius [m]

    Sphere() : center(Vec3::Zero()), radius(1.0) {}

    Sphere(const Vec3& c, double r) : center(c), radius(r) {}

    /**
     * @brief Check if point is inside sphere
     */
    bool contains(const Vec3& point) const {
        return (point - center).squaredNorm() <= radius * radius;
    }

    /**
     * @brief Check if two spheres intersect
     */
    bool intersects(const Sphere& other) const {
        double dist_squared = (center - other.center).squaredNorm();
        double sum_radii = radius + other.radius;
        return dist_squared <= sum_radii * sum_radii;
    }
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Convert State to 13-element vector for integration
 * @param state Input state
 * @return 13D vector: [x, y, z, vx, vy, vz, qw, qx, qy, qz, ωx, ωy, ωz]
 */
inline Eigen::Matrix<double, 13, 1> state_to_vector(const State& state) {
    Eigen::Matrix<double, 13, 1> vec;
    vec.segment<3>(0) = state.position;
    vec.segment<3>(3) = state.velocity;
    vec(6) = state.orientation.w();
    vec(7) = state.orientation.x();
    vec(8) = state.orientation.y();
    vec(9) = state.orientation.z();
    vec.segment<3>(10) = state.angular_velocity;
    return vec;
}

/**
 * @brief Convert 13-element vector back to State
 * @param vec 13D vector
 * @param state Reference to state (will be modified)
 */
inline void vector_to_state(const Eigen::Matrix<double, 13, 1>& vec, State& state) {
    state.position = vec.segment<3>(0);
    state.velocity = vec.segment<3>(3);
    state.orientation = Quat(vec(6), vec(7), vec(8), vec(9)).normalized();
    state.angular_velocity = vec.segment<3>(10);
}

/**
 * @brief Compute kinetic energy of a state
 * @param state State to analyze
 * @return Kinetic energy [J]
 */
inline double kinetic_energy(const State& state) {
    double translational = 0.5 * state.mass * state.velocity.squaredNorm();
    double rotational = 0.5 * state.angular_velocity.dot(state.inertia * state.angular_velocity);
    return translational + rotational;
}

/**
 * @brief Compute linear momentum
 * @param state State to analyze
 * @return Linear momentum [kg·m/s]
 */
inline Vec3 linear_momentum(const State& state) {
    return state.mass * state.velocity;
}

/**
 * @brief Compute angular momentum in body frame
 * @param state State to analyze
 * @return Angular momentum [kg·m²/s]
 */
inline Vec3 angular_momentum(const State& state) {
    return state.inertia * state.angular_velocity;
}

} // namespace physim
