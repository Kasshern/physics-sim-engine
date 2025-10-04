/**
 * @file frame.cpp
 * @brief Implementation of reference frame transformations
 */

#include "physim/core/frame.hpp"
#include "physim/core/constants.hpp"
#include <cmath>
#include <stdexcept>

namespace physim {

// ============================================================================
// Euler Angles
// ============================================================================

void euler_from_rotation(const Mat3& R, double& yaw, double& pitch, double& roll) {
    // Extract Euler angles from rotation matrix (Z-Y-X sequence)
    pitch = std::asin(-R(2, 0));

    if (std::abs(std::cos(pitch)) > 1e-6) {
        yaw = std::atan2(R(1, 0), R(0, 0));
        roll = std::atan2(R(2, 1), R(2, 2));
    } else {
        // Gimbal lock case
        yaw = std::atan2(-R(0, 1), R(1, 1));
        roll = 0.0;
    }
}

// ============================================================================
// ECI <-> ECEF Transformations
// ============================================================================

Mat3 eci_to_ecef_matrix(const Time& time) {
    // Convert time to UT1 for Earth rotation angle
    Time ut1_time = time.to_scale(TimeScale::UT1);
    double theta = gmst(ut1_time.jd());

    // Simple rotation about Z axis (Earth's rotation)
    // For higher precision, add precession, nutation, polar motion
    return rotation_z(theta);
}

Vec3 eci_to_ecef(const Vec3& r_eci, const Time& time) {
    return eci_to_ecef_matrix(time) * r_eci;
}

Vec3 ecef_to_eci(const Vec3& r_ecef, const Time& time) {
    return ecef_to_eci_matrix(time) * r_ecef;
}

// ============================================================================
// Geodetic Coordinates
// ============================================================================

GeodeticCoord ecef_to_geodetic(const Vec3& r_ecef) {
    // Iterative algorithm for WGS84 ellipsoid
    // Based on Bowring (1976) and refined by Fukushima (2006)

    const double a = constants::earth::WGS84_A;
    const double b = constants::earth::WGS84_B;
    const double e2 = constants::earth::WGS84_E2;

    const double x = r_ecef(0);
    const double y = r_ecef(1);
    const double z = r_ecef(2);

    // Longitude is direct
    double lon = std::atan2(y, x);

    // Latitude requires iteration
    double p = std::sqrt(x * x + y * y);
    double lat = std::atan2(z, p * (1.0 - e2)); // Initial guess

    // Iterate to refine latitude
    const int max_iter = 5;
    for (int i = 0; i < max_iter; ++i) {
        double sin_lat = std::sin(lat);
        double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
        double h = p / std::cos(lat) - N;
        lat = std::atan2(z, p * (1.0 - e2 * N / (N + h)));
    }

    // Compute altitude
    double sin_lat = std::sin(lat);
    double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);
    double h = p / std::cos(lat) - N;

    return GeodeticCoord(lat, lon, h);
}

Vec3 geodetic_to_ecef(const GeodeticCoord& geodetic) {
    // WGS84 ellipsoid
    const double a = constants::earth::WGS84_A;
    const double e2 = constants::earth::WGS84_E2;

    const double lat = geodetic.latitude;
    const double lon = geodetic.longitude;
    const double h = geodetic.altitude;

    double sin_lat = std::sin(lat);
    double cos_lat = std::cos(lat);
    double sin_lon = std::sin(lon);
    double cos_lon = std::cos(lon);

    // Radius of curvature in prime vertical
    double N = a / std::sqrt(1.0 - e2 * sin_lat * sin_lat);

    Vec3 r_ecef;
    r_ecef(0) = (N + h) * cos_lat * cos_lon;
    r_ecef(1) = (N + h) * cos_lat * sin_lon;
    r_ecef(2) = (N * (1.0 - e2) + h) * sin_lat;

    return r_ecef;
}

// ============================================================================
// Orbit-Relative Frames
// ============================================================================

Mat3 lvlh_frame(const Vec3& r_eci, const Vec3& v_eci) {
    // LVLH frame:
    // X = -r̂ (nadir, towards Earth center)
    // Z = ĥ (orbit normal)
    // Y = Ẑ × X̂ (along-track, completing right-handed frame)

    Vec3 h = r_eci.cross(v_eci); // Angular momentum vector

    if (h.norm() < 1e-10) {
        throw std::runtime_error("Cannot construct LVLH frame: orbit is degenerate");
    }

    Vec3 z_hat = h.normalized();
    Vec3 x_hat = -r_eci.normalized();
    Vec3 y_hat = z_hat.cross(x_hat);

    Mat3 R_eci_to_lvlh;
    R_eci_to_lvlh.row(0) = x_hat.transpose();
    R_eci_to_lvlh.row(1) = y_hat.transpose();
    R_eci_to_lvlh.row(2) = z_hat.transpose();

    return R_eci_to_lvlh;
}

Mat3 rtn_frame(const Vec3& r_eci, const Vec3& v_eci) {
    // RTN frame:
    // R = r̂ (radial, away from Earth)
    // N = ĥ (orbit normal)
    // T = N̂ × R̂ (tangential, in-plane perpendicular to radial)

    Vec3 h = r_eci.cross(v_eci); // Angular momentum vector

    if (h.norm() < 1e-10) {
        throw std::runtime_error("Cannot construct RTN frame: orbit is degenerate");
    }

    Vec3 r_hat = r_eci.normalized();
    Vec3 n_hat = h.normalized();
    Vec3 t_hat = n_hat.cross(r_hat);

    Mat3 R_eci_to_rtn;
    R_eci_to_rtn.row(0) = r_hat.transpose();
    R_eci_to_rtn.row(1) = t_hat.transpose();
    R_eci_to_rtn.row(2) = n_hat.transpose();

    return R_eci_to_rtn;
}

// ============================================================================
// Topocentric Coordinates
// ============================================================================

TopocentricCoord eci_to_topocentric(const Vec3& r_sat_eci,
                                    const Vec3& r_site_ecef,
                                    const Time& time) {
    // Convert satellite position to ECEF
    Vec3 r_sat_ecef = eci_to_ecef(r_sat_eci, time);

    // Relative position vector
    Vec3 rho_ecef = r_sat_ecef - r_site_ecef;

    // Convert site to geodetic to get local frame
    GeodeticCoord site_geodetic = ecef_to_geodetic(r_site_ecef);

    double lat = site_geodetic.latitude;
    double lon = site_geodetic.longitude;

    // Rotation from ECEF to SEZ (South-East-Zenith) frame
    double sin_lat = std::sin(lat);
    double cos_lat = std::cos(lat);
    double sin_lon = std::sin(lon);
    double cos_lon = std::cos(lon);

    Mat3 R_ecef_to_sez;
    R_ecef_to_sez << sin_lat * cos_lon,  sin_lat * sin_lon, -cos_lat,
                     -sin_lon,            cos_lon,            0.0,
                      cos_lat * cos_lon,  cos_lat * sin_lon,  sin_lat;

    Vec3 rho_sez = R_ecef_to_sez * rho_ecef;

    // Convert SEZ to topocentric (azimuth, elevation, range)
    double range = rho_sez.norm();
    double elevation = std::asin(rho_sez(2) / range);

    // Azimuth measured from North towards East
    double azimuth = std::atan2(rho_sez(1), -rho_sez(0));

    // Normalize azimuth to [0, 2π)
    if (azimuth < 0.0) {
        azimuth += constants::TWO_PI;
    }

    return TopocentricCoord(azimuth, elevation, range);
}

TopocentricCoord compute_look_angles(const GeodeticCoord& observer_geodetic,
                                     const Vec3& satellite_eci,
                                     const Time& time) {
    Vec3 observer_ecef = geodetic_to_ecef(observer_geodetic);
    return eci_to_topocentric(satellite_eci, observer_ecef, time);
}

// ============================================================================
// Generic Frame Transformations
// ============================================================================

Vec3 transform_frame(const Vec3& vec,
                     Frame from_frame,
                     Frame to_frame,
                     const Vec3& r_ref,
                     const Vec3& v_ref,
                     const Time& time) {
    if (from_frame == to_frame) {
        return vec;
    }

    // All transformations go through ECI as intermediate frame
    Vec3 vec_eci = vec;

    // Transform from source frame to ECI
    switch (from_frame) {
        case Frame::ICRF:
        case Frame::ECI:
            // Already in ECI
            break;

        case Frame::ECEF:
            vec_eci = ecef_to_eci(vec, time);
            break;

        case Frame::LVLH:
            if (r_ref.norm() < 1e-10 || v_ref.norm() < 1e-10) {
                throw std::invalid_argument("LVLH frame requires reference position and velocity");
            }
            vec_eci = lvlh_frame(r_ref, v_ref).transpose() * vec;
            break;

        case Frame::RTN:
            if (r_ref.norm() < 1e-10 || v_ref.norm() < 1e-10) {
                throw std::invalid_argument("RTN frame requires reference position and velocity");
            }
            vec_eci = rtn_frame(r_ref, v_ref).transpose() * vec;
            break;

        case Frame::NTW:
            // NTW is similar to RTN but with different ordering
            throw std::runtime_error("NTW frame not yet implemented");
    }

    // Transform from ECI to target frame
    Vec3 result = vec_eci;

    switch (to_frame) {
        case Frame::ICRF:
        case Frame::ECI:
            // Already in ECI
            break;

        case Frame::ECEF:
            result = eci_to_ecef(vec_eci, time);
            break;

        case Frame::LVLH:
            if (r_ref.norm() < 1e-10 || v_ref.norm() < 1e-10) {
                throw std::invalid_argument("LVLH frame requires reference position and velocity");
            }
            result = lvlh_frame(r_ref, v_ref) * vec_eci;
            break;

        case Frame::RTN:
            if (r_ref.norm() < 1e-10 || v_ref.norm() < 1e-10) {
                throw std::invalid_argument("RTN frame requires reference position and velocity");
            }
            result = rtn_frame(r_ref, v_ref) * vec_eci;
            break;

        case Frame::NTW:
            throw std::runtime_error("NTW frame not yet implemented");
    }

    return result;
}

Quat frame_quaternion(Frame from_frame,
                      Frame to_frame,
                      const Vec3& r_ref,
                      const Vec3& v_ref,
                      const Time& time) {
    // Get rotation matrix
    Mat3 R = Mat3::Identity();

    // Compute transformation through ECI
    Vec3 x_axis(1, 0, 0);
    Vec3 y_axis(0, 1, 0);
    Vec3 z_axis(0, 0, 1);

    Vec3 x_transformed = transform_frame(x_axis, from_frame, to_frame, r_ref, v_ref, time);
    Vec3 y_transformed = transform_frame(y_axis, from_frame, to_frame, r_ref, v_ref, time);
    Vec3 z_transformed = transform_frame(z_axis, from_frame, to_frame, r_ref, v_ref, time);

    R.col(0) = x_transformed;
    R.col(1) = y_transformed;
    R.col(2) = z_transformed;

    return Quat(R);
}

} // namespace physim
