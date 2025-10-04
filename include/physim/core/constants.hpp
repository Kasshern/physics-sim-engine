/**
 * @file constants.hpp
 * @brief Physical and mathematical constants for astrodynamics
 *
 * All values are in SI units and sourced from IAU 2015 and CODATA 2018
 * recommended values for use in production astrodynamics software.
 */

#pragma once

#include <cmath>
#include <numbers>

namespace physim {
namespace constants {

// ============================================================================
// Mathematical Constants
// ============================================================================

/// Pi (C++20 std::numbers::pi)
constexpr double PI = std::numbers::pi;

/// 2 * Pi
constexpr double TWO_PI = 2.0 * PI;

/// Pi / 2
constexpr double HALF_PI = 0.5 * PI;

/// Degrees to radians conversion factor
constexpr double DEG_TO_RAD = PI / 180.0;

/// Radians to degrees conversion factor
constexpr double RAD_TO_DEG = 180.0 / PI;

/// Arc-seconds to radians
constexpr double ARCSEC_TO_RAD = PI / 648000.0;

// ============================================================================
// Universal Physical Constants (CODATA 2018)
// ============================================================================

/// Speed of light in vacuum [m/s]
constexpr double SPEED_OF_LIGHT = 299792458.0;

/// Gravitational constant [m³/(kg·s²)]
/// Source: CODATA 2018
constexpr double GRAVITATIONAL_CONSTANT = 6.67430e-11;

/// Astronomical unit [m]
/// IAU 2012 exact definition
constexpr double ASTRONOMICAL_UNIT = 149597870700.0;

/// Parsec [m]
constexpr double PARSEC = 3.0856775814913673e16;

/// Light-year [m]
constexpr double LIGHT_YEAR = SPEED_OF_LIGHT * 365.25 * 86400.0;

// ============================================================================
// Time Constants
// ============================================================================

/// Seconds per minute
constexpr double SECONDS_PER_MINUTE = 60.0;

/// Seconds per hour
constexpr double SECONDS_PER_HOUR = 3600.0;

/// Seconds per day
constexpr double SECONDS_PER_DAY = 86400.0;

/// Julian year [s] (exactly 365.25 days)
constexpr double JULIAN_YEAR = 365.25 * SECONDS_PER_DAY;

/// Julian century [s]
constexpr double JULIAN_CENTURY = 36525.0 * SECONDS_PER_DAY;

/// Tropical year [s] (approximate, varies slightly)
constexpr double TROPICAL_YEAR = 365.242189 * SECONDS_PER_DAY;

/// Sidereal day [s] (Earth's rotation period)
constexpr double SIDEREAL_DAY = 86164.0905;

// ============================================================================
// Solar System Bodies - Gravitational Parameters (GM) [m³/s²]
// ============================================================================
// Source: IAU 2015 / DE440 ephemeris

namespace gm {
    /// Sun GM [m³/s²]
    constexpr double SUN = 1.32712440018e20;

    /// Mercury GM [m³/s²]
    constexpr double MERCURY = 2.2032e13;

    /// Venus GM [m³/s²]
    constexpr double VENUS = 3.24859e14;

    /// Earth GM [m³/s²] (includes atmosphere)
    constexpr double EARTH = 3.986004418e14;

    /// Moon GM [m³/s²]
    constexpr double MOON = 4.9028e12;

    /// Mars GM [m³/s²]
    constexpr double MARS = 4.282837e13;

    /// Jupiter GM [m³/s²]
    constexpr double JUPITER = 1.26686534e17;

    /// Saturn GM [m³/s²]
    constexpr double SATURN = 3.7931187e16;

    /// Uranus GM [m³/s²]
    constexpr double URANUS = 5.793939e15;

    /// Neptune GM [m³/s²]
    constexpr double NEPTUNE = 6.836529e15;

    /// Pluto GM [m³/s²]
    constexpr double PLUTO = 8.71e11;
} // namespace gm

// ============================================================================
// Solar System Bodies - Radii [m]
// ============================================================================
// Source: IAU 2015 / NASA fact sheets

namespace radius {
    /// Sun mean radius [m]
    constexpr double SUN = 6.96e8;

    /// Mercury mean radius [m]
    constexpr double MERCURY = 2.4397e6;

    /// Venus mean radius [m]
    constexpr double VENUS = 6.0518e6;

    /// Earth equatorial radius [m]
    constexpr double EARTH_EQUATORIAL = 6.378137e6;

    /// Earth polar radius [m]
    constexpr double EARTH_POLAR = 6.356752e6;

    /// Earth mean radius [m]
    constexpr double EARTH = 6.371e6;

    /// Moon mean radius [m]
    constexpr double MOON = 1.7374e6;

    /// Mars mean radius [m]
    constexpr double MARS = 3.3895e6;

    /// Jupiter equatorial radius [m]
    constexpr double JUPITER = 7.1492e7;

    /// Saturn equatorial radius [m]
    constexpr double SATURN = 6.0268e7;

    /// Uranus equatorial radius [m]
    constexpr double URANUS = 2.5559e7;

    /// Neptune equatorial radius [m]
    constexpr double NEPTUNE = 2.4764e7;

    /// Pluto mean radius [m]
    constexpr double PLUTO = 1.188e6;
} // namespace radius

// ============================================================================
// Solar System Bodies - Masses [kg]
// ============================================================================
// Derived from GM / G

namespace mass {
    /// Sun mass [kg]
    constexpr double SUN = gm::SUN / GRAVITATIONAL_CONSTANT;

    /// Earth mass [kg]
    constexpr double EARTH = gm::EARTH / GRAVITATIONAL_CONSTANT;

    /// Moon mass [kg]
    constexpr double MOON = gm::MOON / GRAVITATIONAL_CONSTANT;

    /// Jupiter mass [kg]
    constexpr double JUPITER = gm::JUPITER / GRAVITATIONAL_CONSTANT;
} // namespace mass

// ============================================================================
// Earth Geophysical Constants
// ============================================================================

namespace earth {
    /// WGS84 semi-major axis (equatorial radius) [m]
    constexpr double WGS84_A = 6378137.0;

    /// WGS84 semi-minor axis (polar radius) [m]
    constexpr double WGS84_B = 6356752.314245;

    /// WGS84 flattening factor
    constexpr double WGS84_F = 1.0 / 298.257223563;

    /// WGS84 first eccentricity squared
    constexpr double WGS84_E2 = 2.0 * WGS84_F - WGS84_F * WGS84_F;

    /// Earth's angular velocity [rad/s]
    constexpr double ANGULAR_VELOCITY = 7.2921150e-5;

    /// J2 zonal harmonic coefficient (Earth oblateness)
    /// Source: EGM2008
    constexpr double J2 = 1.0826267e-3;

    /// J3 zonal harmonic coefficient
    constexpr double J3 = -2.53266e-6;

    /// J4 zonal harmonic coefficient
    constexpr double J4 = -1.61994e-6;

    /// Mean solar day on Earth [s]
    constexpr double SOLAR_DAY = SECONDS_PER_DAY;

    /// Sidereal day on Earth [s]
    constexpr double SIDEREAL_DAY_EARTH = physim::constants::SIDEREAL_DAY;
} // namespace earth

// ============================================================================
// Orbital Mechanics Constants
// ============================================================================

/// Geostationary orbit altitude [m] above Earth surface
constexpr double GEO_ALTITUDE = 35786000.0;

/// Geostationary orbit radius [m] from Earth center
constexpr double GEO_RADIUS = earth::WGS84_A + GEO_ALTITUDE;

/// Low Earth Orbit typical altitude [m]
constexpr double LEO_ALTITUDE = 400000.0;

/// ISS nominal altitude [m]
constexpr double ISS_ALTITUDE = 408000.0;

/// Escape velocity from Earth surface [m/s]
constexpr double EARTH_ESCAPE_VELOCITY = 11186.0;

// ============================================================================
// Solar Radiation and Atmospheric Constants
// ============================================================================

/// Solar constant (irradiance at 1 AU) [W/m²]
/// Source: WMO 2015
constexpr double SOLAR_CONSTANT = 1361.0;

/// Solar radiation pressure at 1 AU [N/m²] for perfect reflector
/// P = F/c where F is solar flux
constexpr double SOLAR_PRESSURE_1AU = SOLAR_CONSTANT / SPEED_OF_LIGHT;

/// Standard atmospheric pressure at sea level [Pa]
constexpr double STANDARD_PRESSURE = 101325.0;

/// Standard temperature at sea level [K]
constexpr double STANDARD_TEMPERATURE = 288.15;

/// Standard atmospheric density at sea level [kg/m³]
constexpr double STANDARD_DENSITY = 1.225;

/// Atmospheric scale height [m] (exponential model)
constexpr double SCALE_HEIGHT = 8500.0;

// ============================================================================
// Numerical Tolerance Constants
// ============================================================================

/// Machine epsilon for double precision
constexpr double EPSILON = 1e-15;

/// Default tolerance for iterative solvers
constexpr double DEFAULT_TOLERANCE = 1e-12;

/// Tolerance for quaternion normalization checks
constexpr double QUAT_NORMALIZATION_TOL = 1e-10;

/// Tolerance for energy conservation checks
constexpr double ENERGY_CONSERVATION_TOL = 1e-10;

/// Small angle threshold for series approximations [rad]
constexpr double SMALL_ANGLE_THRESHOLD = 1e-8;

// ============================================================================
// Conversion Utilities
// ============================================================================

/**
 * @brief Convert kilometers to meters
 */
constexpr double km_to_m(double km) {
    return km * 1000.0;
}

/**
 * @brief Convert meters to kilometers
 */
constexpr double m_to_km(double m) {
    return m / 1000.0;
}

/**
 * @brief Convert AU to meters
 */
constexpr double au_to_m(double au) {
    return au * ASTRONOMICAL_UNIT;
}

/**
 * @brief Convert meters to AU
 */
constexpr double m_to_au(double m) {
    return m / ASTRONOMICAL_UNIT;
}

/**
 * @brief Convert degrees to radians
 */
constexpr double deg_to_rad(double deg) {
    return deg * DEG_TO_RAD;
}

/**
 * @brief Convert radians to degrees
 */
constexpr double rad_to_deg(double rad) {
    return rad * RAD_TO_DEG;
}

/**
 * @brief Convert days to seconds
 */
constexpr double days_to_seconds(double days) {
    return days * SECONDS_PER_DAY;
}

/**
 * @brief Convert seconds to days
 */
constexpr double seconds_to_days(double seconds) {
    return seconds / SECONDS_PER_DAY;
}

/**
 * @brief Convert years to seconds (Julian year)
 */
constexpr double years_to_seconds(double years) {
    return years * JULIAN_YEAR;
}

/**
 * @brief Convert seconds to years (Julian year)
 */
constexpr double seconds_to_years(double seconds) {
    return seconds / JULIAN_YEAR;
}

} // namespace constants
} // namespace physim
