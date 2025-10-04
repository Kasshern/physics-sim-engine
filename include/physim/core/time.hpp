/**
 * @file time.hpp
 * @brief Time systems and conversions for astrodynamics
 *
 * Provides functionality for working with various time systems used in
 * orbital mechanics, including UTC, TAI, TT, TDB, and conversions between them.
 */

#pragma once

#include <chrono>
#include <string>
#include <cmath>

namespace physim {

// ============================================================================
// Time System Definitions
// ============================================================================

/**
 * @brief Time scale enumeration
 *
 * Different time systems used in astrodynamics:
 * - UTC: Coordinated Universal Time (civil time with leap seconds)
 * - TAI: International Atomic Time (continuous, no leap seconds)
 * - TT: Terrestrial Time (TAI + 32.184s)
 * - TDB: Barycentric Dynamical Time (solar system barycenter)
 * - UT1: Universal Time (based on Earth rotation)
 */
enum class TimeScale {
    UTC,  ///< Coordinated Universal Time
    TAI,  ///< International Atomic Time
    TT,   ///< Terrestrial Time
    TDB,  ///< Barycentric Dynamical Time
    UT1   ///< Universal Time
};

// ============================================================================
// Epoch Definitions
// ============================================================================

/// J2000.0 epoch (2000-01-01 12:00:00 TT) in Julian Date
constexpr double J2000_JD = 2451545.0;

/// Modified Julian Date zero point (MJD = JD - 2400000.5)
constexpr double MJD_OFFSET = 2400000.5;

/// GPS epoch (1980-01-06 00:00:00 UTC) in Julian Date
constexpr double GPS_EPOCH_JD = 2444244.5;

/// Unix epoch (1970-01-01 00:00:00 UTC) in Julian Date
constexpr double UNIX_EPOCH_JD = 2440587.5;

/// TAI - UTC offset at J2000.0 [s] (32 leap seconds)
constexpr double TAI_UTC_J2000 = 32.0;

/// TT - TAI constant offset [s]
constexpr double TT_TAI_OFFSET = 32.184;

/// TT - UTC at J2000.0 [s]
constexpr double TT_UTC_J2000 = TAI_UTC_J2000 + TT_TAI_OFFSET;

// ============================================================================
// Time Class
// ============================================================================

/**
 * @brief High-precision time representation for astrodynamics
 *
 * Stores time as Julian Date (JD) in a specified time scale.
 * Provides conversions between time systems and utilities for
 * time arithmetic and formatting.
 */
class Time {
public:
    /**
     * @brief Default constructor (J2000.0 epoch in TT)
     */
    Time() : jd_(J2000_JD), scale_(TimeScale::TT) {}

    /**
     * @brief Construct from Julian Date
     * @param julian_date Julian Date
     * @param scale Time scale
     */
    explicit Time(double julian_date, TimeScale scale = TimeScale::TT)
        : jd_(julian_date), scale_(scale) {}

    /**
     * @brief Construct from Modified Julian Date
     * @param mjd Modified Julian Date
     * @param scale Time scale
     * @return Time object
     */
    static Time from_mjd(double mjd, TimeScale scale = TimeScale::TT) {
        return Time(mjd + MJD_OFFSET, scale);
    }

    /**
     * @brief Construct from calendar date
     * @param year Year
     * @param month Month (1-12)
     * @param day Day (1-31)
     * @param hour Hour (0-23)
     * @param minute Minute (0-59)
     * @param second Second (0-59.999...)
     * @param scale Time scale
     * @return Time object
     */
    static Time from_calendar(int year, int month, int day,
                             int hour = 0, int minute = 0, double second = 0.0,
                             TimeScale scale = TimeScale::TT);

    /**
     * @brief Construct from ISO 8601 string
     * @param iso_string String in format "YYYY-MM-DDTHH:MM:SS.SSS"
     * @param scale Time scale
     * @return Time object
     */
    static Time from_iso_string(const std::string& iso_string,
                               TimeScale scale = TimeScale::TT);

    /**
     * @brief Construct from Unix timestamp
     * @param unix_time Seconds since 1970-01-01 00:00:00 UTC
     * @return Time object in UTC
     */
    static Time from_unix(double unix_time) {
        return Time(UNIX_EPOCH_JD + unix_time / 86400.0, TimeScale::UTC);
    }

    /**
     * @brief Get Julian Date
     */
    double jd() const { return jd_; }

    /**
     * @brief Get Modified Julian Date
     */
    double mjd() const { return jd_ - MJD_OFFSET; }

    /**
     * @brief Get time scale
     */
    TimeScale scale() const { return scale_; }

    /**
     * @brief Get Julian centuries since J2000.0
     * @return Centuries (1 century = 36525 days)
     */
    double centuries_since_j2000() const {
        return (jd_ - J2000_JD) / 36525.0;
    }

    /**
     * @brief Get days since J2000.0
     */
    double days_since_j2000() const {
        return jd_ - J2000_JD;
    }

    /**
     * @brief Get seconds since J2000.0
     */
    double seconds_since_j2000() const {
        return (jd_ - J2000_JD) * 86400.0;
    }

    /**
     * @brief Convert to different time scale
     * @param target_scale Target time scale
     * @return New Time object in target scale
     */
    Time to_scale(TimeScale target_scale) const;

    /**
     * @brief Convert to ISO 8601 string
     * @return String in format "YYYY-MM-DDTHH:MM:SS.SSS"
     */
    std::string to_iso_string() const;

    /**
     * @brief Convert to calendar date components
     * @param year Output year
     * @param month Output month (1-12)
     * @param day Output day (1-31)
     * @param hour Output hour (0-23)
     * @param minute Output minute (0-59)
     * @param second Output second (0-59.999...)
     */
    void to_calendar(int& year, int& month, int& day,
                    int& hour, int& minute, double& second) const;

    // Arithmetic operators
    Time operator+(double days) const { return Time(jd_ + days, scale_); }
    Time operator-(double days) const { return Time(jd_ - days, scale_); }
    double operator-(const Time& other) const { return jd_ - other.jd_; }
    Time& operator+=(double days) { jd_ += days; return *this; }
    Time& operator-=(double days) { jd_ -= days; return *this; }

    // Comparison operators
    bool operator==(const Time& other) const {
        return std::abs(jd_ - other.jd_) < 1e-9 && scale_ == other.scale_;
    }
    bool operator!=(const Time& other) const { return !(*this == other); }
    bool operator<(const Time& other) const { return jd_ < other.jd_; }
    bool operator<=(const Time& other) const { return jd_ <= other.jd_; }
    bool operator>(const Time& other) const { return jd_ > other.jd_; }
    bool operator>=(const Time& other) const { return jd_ >= other.jd_; }

private:
    double jd_;          ///< Julian Date
    TimeScale scale_;    ///< Time scale
};

// ============================================================================
// Utility Functions
// ============================================================================

/**
 * @brief Get current system time as Time object
 * @param scale Time scale to use (default UTC)
 * @return Current time
 */
Time now(TimeScale scale = TimeScale::UTC);

/**
 * @brief Convert UTC to TAI
 * @param utc_jd Julian Date in UTC
 * @return Julian Date in TAI
 *
 * Note: This is a simplified conversion assuming constant offset.
 * Production code should use a leap second table (IERS Bulletin C).
 */
double utc_to_tai(double utc_jd);

/**
 * @brief Convert TAI to UTC
 * @param tai_jd Julian Date in TAI
 * @return Julian Date in UTC
 */
double tai_to_utc(double tai_jd);

/**
 * @brief Convert TAI to TT (Terrestrial Time)
 * @param tai_jd Julian Date in TAI
 * @return Julian Date in TT
 */
inline double tai_to_tt(double tai_jd) {
    return tai_jd + TT_TAI_OFFSET / 86400.0;
}

/**
 * @brief Convert TT to TAI
 * @param tt_jd Julian Date in TT
 * @return Julian Date in TAI
 */
inline double tt_to_tai(double tt_jd) {
    return tt_jd - TT_TAI_OFFSET / 86400.0;
}

/**
 * @brief Convert TT to TDB (Barycentric Dynamical Time)
 * @param tt_jd Julian Date in TT
 * @return Julian Date in TDB
 *
 * Uses simplified Fairhead & Bretagnon 1990 algorithm.
 * For high-precision work, use full JPL ephemeris.
 */
double tt_to_tdb(double tt_jd);

/**
 * @brief Convert TDB to TT
 * @param tdb_jd Julian Date in TDB
 * @return Julian Date in TT
 */
double tdb_to_tt(double tdb_jd);

/**
 * @brief Get number of leap seconds at given date
 * @param utc_jd Julian Date in UTC
 * @return Number of leap seconds (TAI - UTC)
 *
 * Simplified version with hardcoded values.
 * Production should load from IERS data file.
 */
int get_leap_seconds(double utc_jd);

/**
 * @brief Convert calendar date to Julian Date
 * @param year Year
 * @param month Month (1-12)
 * @param day Day (1-31)
 * @param hour Hour (0-23)
 * @param minute Minute (0-59)
 * @param second Second (0-59.999...)
 * @return Julian Date
 */
double calendar_to_jd(int year, int month, int day,
                      int hour = 0, int minute = 0, double second = 0.0);

/**
 * @brief Convert Julian Date to calendar date
 * @param jd Julian Date
 * @param year Output year
 * @param month Output month (1-12)
 * @param day Output day (1-31)
 * @param hour Output hour (0-23)
 * @param minute Output minute (0-59)
 * @param second Output second (0-59.999...)
 */
void jd_to_calendar(double jd, int& year, int& month, int& day,
                   int& hour, int& minute, double& second);

/**
 * @brief Check if year is a leap year
 * @param year Year
 * @return true if leap year, false otherwise
 */
inline bool is_leap_year(int year) {
    return (year % 4 == 0 && year % 100 != 0) || (year % 400 == 0);
}

/**
 * @brief Get number of days in month
 * @param year Year
 * @param month Month (1-12)
 * @return Number of days in month
 */
int days_in_month(int year, int month);

/**
 * @brief Greenwich Mean Sidereal Time
 * @param ut1_jd Julian Date in UT1
 * @return GMST [radians]
 *
 * Computes GMST from UT1 using IAU 2006 precession model.
 */
double gmst(double ut1_jd);

/**
 * @brief Greenwich Apparent Sidereal Time
 * @param ut1_jd Julian Date in UT1
 * @return GAST [radians]
 *
 * Includes nutation correction (equation of equinoxes).
 */
double gast(double ut1_jd);

/**
 * @brief Local Mean Sidereal Time
 * @param ut1_jd Julian Date in UT1
 * @param longitude East longitude [radians]
 * @return LMST [radians]
 */
inline double lmst(double ut1_jd, double longitude) {
    return gmst(ut1_jd) + longitude;
}

} // namespace physim
