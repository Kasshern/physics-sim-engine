/**
 * @file frame.hpp
 * @brief Reference frame transformations for orbital mechanics
 *
 * Provides coordinate frame definitions and transformations between
 * common reference frames used in astrodynamics (ECI, ECEF, LVLH, etc.)
 */

#pragma once

#include "types.hpp"
#include "time.hpp"
#include "constants.hpp"
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace physim {

// ============================================================================
// Reference Frame Definitions
// ============================================================================

/**
 * @brief Reference frame enumeration
 *
 * Common coordinate frames in orbital mechanics:
 * - ICRF: International Celestial Reference Frame (inertial, J2000.0)
 * - ECI: Earth-Centered Inertial (J2000.0 equator and equinox)
 * - ECEF: Earth-Centered Earth-Fixed (rotates with Earth)
 * - LVLH: Local Vertical Local Horizontal (orbit frame)
 * - RTN: Radial-Tangential-Normal (RSW frame)
 * - NTW: Normal-Tangential-out-of-plane (satellite frame)
 */
enum class Frame {
    ICRF,   ///< International Celestial Reference Frame
    ECI,    ///< Earth-Centered Inertial (J2000.0)
    ECEF,   ///< Earth-Centered Earth-Fixed (ITRF)
    LVLH,   ///< Local Vertical Local Horizontal
    RTN,    ///< Radial-Tangential-Normal
    NTW     ///< Normal-Tangential-out-of-plane
};

// ============================================================================
// Rotation Matrix Generators
// ============================================================================

/**
 * @brief Create rotation matrix about X axis
 * @param angle Rotation angle [radians]
 * @return 3x3 rotation matrix
 */
inline Mat3 rotation_x(double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    Mat3 R;
    R << 1.0,  0.0,  0.0,
         0.0,    c,   -s,
         0.0,    s,    c;
    return R;
}

/**
 * @brief Create rotation matrix about Y axis
 * @param angle Rotation angle [radians]
 * @return 3x3 rotation matrix
 */
inline Mat3 rotation_y(double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    Mat3 R;
    R <<   c,  0.0,    s,
         0.0,  1.0,  0.0,
          -s,  0.0,    c;
    return R;
}

/**
 * @brief Create rotation matrix about Z axis
 * @param angle Rotation angle [radians]
 * @return 3x3 rotation matrix
 */
inline Mat3 rotation_z(double angle) {
    const double c = std::cos(angle);
    const double s = std::sin(angle);
    Mat3 R;
    R <<   c,   -s,  0.0,
           s,    c,  0.0,
         0.0,  0.0,  1.0;
    return R;
}

/**
 * @brief Create rotation matrix from Euler angles (Z-Y-X sequence)
 * @param yaw Rotation about Z [radians]
 * @param pitch Rotation about Y [radians]
 * @param roll Rotation about X [radians]
 * @return 3x3 rotation matrix
 */
inline Mat3 rotation_from_euler(double yaw, double pitch, double roll) {
    return rotation_z(yaw) * rotation_y(pitch) * rotation_x(roll);
}

/**
 * @brief Extract Euler angles from rotation matrix (Z-Y-X sequence)
 * @param R Rotation matrix
 * @param yaw Output yaw [radians]
 * @param pitch Output pitch [radians]
 * @param roll Output roll [radians]
 */
void euler_from_rotation(const Mat3& R, double& yaw, double& pitch, double& roll);

// ============================================================================
// ECI <-> ECEF Transformations
// ============================================================================

/**
 * @brief Transform from ECI to ECEF
 * @param r_eci Position in ECI [m]
 * @param time Current time (UT1 for Earth rotation)
 * @return Position in ECEF [m]
 *
 * Applies Earth rotation using Greenwich Mean Sidereal Time (GMST).
 * For high-precision work, also apply polar motion and nutation.
 */
Vec3 eci_to_ecef(const Vec3& r_eci, const Time& time);

/**
 * @brief Transform from ECEF to ECI
 * @param r_ecef Position in ECEF [m]
 * @param time Current time (UT1 for Earth rotation)
 * @return Position in ECI [m]
 */
Vec3 ecef_to_eci(const Vec3& r_ecef, const Time& time);

/**
 * @brief Get rotation matrix from ECI to ECEF
 * @param time Current time (UT1)
 * @return 3x3 rotation matrix
 */
Mat3 eci_to_ecef_matrix(const Time& time);

/**
 * @brief Get rotation matrix from ECEF to ECI
 * @param time Current time (UT1)
 * @return 3x3 rotation matrix
 */
inline Mat3 ecef_to_eci_matrix(const Time& time) {
    return eci_to_ecef_matrix(time).transpose();
}

// ============================================================================
// Geodetic Coordinates
// ============================================================================

/**
 * @brief Geodetic coordinates (latitude, longitude, altitude)
 */
struct GeodeticCoord {
    double latitude;   ///< Latitude [radians] (-π/2 to π/2)
    double longitude;  ///< Longitude [radians] (-π to π)
    double altitude;   ///< Altitude above reference ellipsoid [m]

    GeodeticCoord()
        : latitude(0.0), longitude(0.0), altitude(0.0) {}

    GeodeticCoord(double lat, double lon, double alt)
        : latitude(lat), longitude(lon), altitude(alt) {}
};

/**
 * @brief Convert ECEF Cartesian to geodetic coordinates
 * @param r_ecef Position in ECEF [m]
 * @return Geodetic coordinates (WGS84 ellipsoid)
 *
 * Uses iterative algorithm (converges in 3-4 iterations).
 */
GeodeticCoord ecef_to_geodetic(const Vec3& r_ecef);

/**
 * @brief Convert geodetic coordinates to ECEF Cartesian
 * @param geodetic Geodetic coordinates (WGS84)
 * @return Position in ECEF [m]
 */
Vec3 geodetic_to_ecef(const GeodeticCoord& geodetic);

// ============================================================================
// Orbit-Relative Frames
// ============================================================================

/**
 * @brief Create LVLH (Local Vertical Local Horizontal) frame
 * @param r_eci Position in ECI [m]
 * @param v_eci Velocity in ECI [m/s]
 * @return Rotation matrix from ECI to LVLH
 *
 * LVLH frame definition:
 * - X: Radial (nadir direction, -r̂)
 * - Y: Along-track (velocity direction projected to orbital plane)
 * - Z: Cross-track (orbit normal, h_hat)
 */
Mat3 lvlh_frame(const Vec3& r_eci, const Vec3& v_eci);

/**
 * @brief Create RTN (Radial-Tangential-Normal) frame
 * @param r_eci Position in ECI [m]
 * @param v_eci Velocity in ECI [m/s]
 * @return Rotation matrix from ECI to RTN
 *
 * RTN frame definition:
 * - R: Radial (r̂)
 * - T: Tangential (perpendicular to R in orbital plane)
 * - N: Normal (orbit angular momentum direction)
 */
Mat3 rtn_frame(const Vec3& r_eci, const Vec3& v_eci);

/**
 * @brief Transform vector from ECI to LVLH
 * @param vec_eci Vector in ECI
 * @param r_eci Reference position in ECI [m]
 * @param v_eci Reference velocity in ECI [m/s]
 * @return Vector in LVLH frame
 */
inline Vec3 eci_to_lvlh(const Vec3& vec_eci, const Vec3& r_eci, const Vec3& v_eci) {
    return lvlh_frame(r_eci, v_eci) * vec_eci;
}

/**
 * @brief Transform vector from LVLH to ECI
 * @param vec_lvlh Vector in LVLH
 * @param r_eci Reference position in ECI [m]
 * @param v_eci Reference velocity in ECI [m/s]
 * @return Vector in ECI frame
 */
inline Vec3 lvlh_to_eci(const Vec3& vec_lvlh, const Vec3& r_eci, const Vec3& v_eci) {
    return lvlh_frame(r_eci, v_eci).transpose() * vec_lvlh;
}

/**
 * @brief Transform vector from ECI to RTN
 * @param vec_eci Vector in ECI
 * @param r_eci Reference position in ECI [m]
 * @param v_eci Reference velocity in ECI [m/s]
 * @return Vector in RTN frame
 */
inline Vec3 eci_to_rtn(const Vec3& vec_eci, const Vec3& r_eci, const Vec3& v_eci) {
    return rtn_frame(r_eci, v_eci) * vec_eci;
}

/**
 * @brief Transform vector from RTN to ECI
 * @param vec_rtn Vector in RTN
 * @param r_eci Reference position in ECI [m]
 * @param v_eci Reference velocity in ECI [m/s]
 * @return Vector in ECI frame
 */
inline Vec3 rtn_to_eci(const Vec3& vec_rtn, const Vec3& r_eci, const Vec3& v_eci) {
    return rtn_frame(r_eci, v_eci).transpose() * vec_rtn;
}

// ============================================================================
// Topocentric Coordinates (for ground stations)
// ============================================================================

/**
 * @brief Topocentric coordinates (azimuth, elevation, range)
 */
struct TopocentricCoord {
    double azimuth;    ///< Azimuth [radians] (0 = North, π/2 = East)
    double elevation;  ///< Elevation [radians] (0 = horizon, π/2 = zenith)
    double range;      ///< Slant range [m]

    TopocentricCoord()
        : azimuth(0.0), elevation(0.0), range(0.0) {}

    TopocentricCoord(double az, double el, double r)
        : azimuth(az), elevation(el), range(r) {}
};

/**
 * @brief Convert ECI position to topocentric coordinates
 * @param r_sat_eci Satellite position in ECI [m]
 * @param r_site_ecef Ground station position in ECEF [m]
 * @param time Current time
 * @return Topocentric coordinates (azimuth, elevation, range)
 *
 * Used for satellite tracking from ground stations.
 */
TopocentricCoord eci_to_topocentric(const Vec3& r_sat_eci,
                                    const Vec3& r_site_ecef,
                                    const Time& time);

/**
 * @brief Compute look angles from ground station to satellite
 * @param observer_geodetic Observer geodetic coordinates
 * @param satellite_eci Satellite position in ECI [m]
 * @param time Current time
 * @return Topocentric coordinates
 */
TopocentricCoord compute_look_angles(const GeodeticCoord& observer_geodetic,
                                     const Vec3& satellite_eci,
                                     const Time& time);

// ============================================================================
// Frame Transformation Utilities
// ============================================================================

/**
 * @brief Generic frame transformation
 * @param vec Input vector
 * @param from_frame Source frame
 * @param to_frame Target frame
 * @param r_ref Reference position (for orbit frames)
 * @param v_ref Reference velocity (for orbit frames)
 * @param time Time (for time-dependent frames like ECEF)
 * @return Transformed vector
 */
Vec3 transform_frame(const Vec3& vec,
                     Frame from_frame,
                     Frame to_frame,
                     const Vec3& r_ref = Vec3::Zero(),
                     const Vec3& v_ref = Vec3::Zero(),
                     const Time& time = Time());

/**
 * @brief Check if frame transformation requires reference orbit
 * @param frame Frame to check
 * @return true if frame needs reference position/velocity
 */
inline bool requires_orbit_reference(Frame frame) {
    return frame == Frame::LVLH || frame == Frame::RTN || frame == Frame::NTW;
}

/**
 * @brief Check if frame transformation requires time
 * @param frame Frame to check
 * @return true if frame is time-dependent
 */
inline bool requires_time(Frame frame) {
    return frame == Frame::ECEF;
}

/**
 * @brief Compute quaternion from one frame to another
 * @param from_frame Source frame
 * @param to_frame Target frame
 * @param r_ref Reference position (for orbit frames)
 * @param v_ref Reference velocity (for orbit frames)
 * @param time Time (for time-dependent frames)
 * @return Rotation quaternion
 */
Quat frame_quaternion(Frame from_frame,
                      Frame to_frame,
                      const Vec3& r_ref = Vec3::Zero(),
                      const Vec3& v_ref = Vec3::Zero(),
                      const Time& time = Time());

} // namespace physim
