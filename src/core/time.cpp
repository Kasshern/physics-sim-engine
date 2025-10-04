/**
 * @file time.cpp
 * @brief Implementation of time system conversions
 */

#include "physim/core/time.hpp"
#include "physim/core/constants.hpp"
#include <cmath>
#include <chrono>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace physim {

// ============================================================================
// Time Class Implementation
// ============================================================================

Time Time::from_calendar(int year, int month, int day,
                         int hour, int minute, double second,
                         TimeScale scale) {
    double jd = calendar_to_jd(year, month, day, hour, minute, second);
    return Time(jd, scale);
}

Time Time::from_iso_string(const std::string& iso_string, TimeScale scale) {
    // Parse ISO 8601 format: YYYY-MM-DDTHH:MM:SS.SSS
    int year, month, day, hour, minute;
    double second;

    char sep;
    std::istringstream iss(iso_string);
    iss >> year >> sep >> month >> sep >> day >> sep
        >> hour >> sep >> minute >> sep >> second;

    if (iss.fail()) {
        throw std::invalid_argument("Invalid ISO 8601 string format");
    }

    return from_calendar(year, month, day, hour, minute, second, scale);
}

Time Time::to_scale(TimeScale target_scale) const {
    if (scale_ == target_scale) {
        return *this;
    }

    // Convert to TT as intermediate (most conversions go through TT)
    double tt_jd = jd_;

    // Source to TT
    switch (scale_) {
        case TimeScale::UTC:
            tt_jd = tai_to_tt(utc_to_tai(jd_));
            break;
        case TimeScale::TAI:
            tt_jd = tai_to_tt(jd_);
            break;
        case TimeScale::TT:
            // Already in TT
            break;
        case TimeScale::TDB:
            tt_jd = tdb_to_tt(jd_);
            break;
        case TimeScale::UT1:
            // Approximate UT1 ≈ UTC for now
            tt_jd = tai_to_tt(utc_to_tai(jd_));
            break;
    }

    // TT to target
    double target_jd = tt_jd;
    switch (target_scale) {
        case TimeScale::UTC:
            target_jd = tai_to_utc(tt_to_tai(tt_jd));
            break;
        case TimeScale::TAI:
            target_jd = tt_to_tai(tt_jd);
            break;
        case TimeScale::TT:
            // Already in TT
            break;
        case TimeScale::TDB:
            target_jd = tt_to_tdb(tt_jd);
            break;
        case TimeScale::UT1:
            // Approximate UT1 ≈ UTC for now
            target_jd = tai_to_utc(tt_to_tai(tt_jd));
            break;
    }

    return Time(target_jd, target_scale);
}

std::string Time::to_iso_string() const {
    int year, month, day, hour, minute;
    double second;
    to_calendar(year, month, day, hour, minute, second);

    std::ostringstream oss;
    oss << std::setfill('0')
        << std::setw(4) << year << "-"
        << std::setw(2) << month << "-"
        << std::setw(2) << day << "T"
        << std::setw(2) << hour << ":"
        << std::setw(2) << minute << ":"
        << std::setw(2) << static_cast<int>(second) << "."
        << std::setw(3) << static_cast<int>((second - static_cast<int>(second)) * 1000.0);

    return oss.str();
}

void Time::to_calendar(int& year, int& month, int& day,
                      int& hour, int& minute, double& second) const {
    jd_to_calendar(jd_, year, month, day, hour, minute, second);
}

// ============================================================================
// Utility Functions Implementation
// ============================================================================

Time now(TimeScale scale) {
    auto now_time = std::chrono::system_clock::now();
    auto unix_time = std::chrono::duration_cast<std::chrono::microseconds>(
        now_time.time_since_epoch()).count() / 1000000.0;

    return Time::from_unix(unix_time).to_scale(scale);
}

double utc_to_tai(double utc_jd) {
    // Get number of leap seconds at this date
    int leap_seconds = get_leap_seconds(utc_jd);
    return utc_jd + static_cast<double>(leap_seconds) / 86400.0;
}

double tai_to_utc(double tai_jd) {
    // Iterate to find correct leap second count
    double utc_approx = tai_jd - TAI_UTC_J2000 / 86400.0;
    int leap_seconds = get_leap_seconds(utc_approx);
    return tai_jd - static_cast<double>(leap_seconds) / 86400.0;
}

double tt_to_tdb(double tt_jd) {
    // Simplified Fairhead & Bretagnon 1990 algorithm
    // For high precision, use full JPL ephemeris

    double T = (tt_jd - J2000_JD) / 36525.0; // Julian centuries since J2000

    // Mean anomaly of Earth (radians)
    double M = 6.240075967 + 628.3019552 * T;

    // TDB - TT in seconds (periodic term)
    double delta_t = 0.001658 * std::sin(M) + 0.000014 * std::sin(2.0 * M);

    return tt_jd + delta_t / 86400.0;
}

double tdb_to_tt(double tdb_jd) {
    // Inverse of tt_to_tdb (iterate if needed, but one iteration is usually sufficient)
    double tt_approx = tdb_jd;

    double T = (tt_approx - J2000_JD) / 36525.0;
    double M = 6.240075967 + 628.3019552 * T;
    double delta_t = 0.001658 * std::sin(M) + 0.000014 * std::sin(2.0 * M);

    return tdb_jd - delta_t / 86400.0;
}

int get_leap_seconds(double utc_jd) {
    // Simplified leap second table
    // Production code should load from IERS Bulletin C

    // Reference: https://www.iers.org/IERS/EN/Publications/Bulletins/bulletins.html

    if (utc_jd < 2441317.5) return 10;      // Before 1972-01-01
    if (utc_jd < 2441499.5) return 10;      // 1972-01-01
    if (utc_jd < 2441683.5) return 11;      // 1972-07-01
    if (utc_jd < 2442048.5) return 12;      // 1973-01-01
    if (utc_jd < 2442413.5) return 13;      // 1974-01-01
    if (utc_jd < 2442778.5) return 14;      // 1975-01-01
    if (utc_jd < 2443144.5) return 15;      // 1976-01-01
    if (utc_jd < 2443509.5) return 16;      // 1977-01-01
    if (utc_jd < 2443874.5) return 17;      // 1978-01-01
    if (utc_jd < 2444239.5) return 18;      // 1979-01-01
    if (utc_jd < 2444786.5) return 19;      // 1980-01-01
    if (utc_jd < 2445151.5) return 20;      // 1981-07-01
    if (utc_jd < 2445516.5) return 21;      // 1982-07-01
    if (utc_jd < 2446247.5) return 22;      // 1983-07-01
    if (utc_jd < 2447161.5) return 23;      // 1985-07-01
    if (utc_jd < 2447892.5) return 24;      // 1988-01-01
    if (utc_jd < 2448257.5) return 25;      // 1990-01-01
    if (utc_jd < 2448804.5) return 26;      // 1991-01-01
    if (utc_jd < 2449169.5) return 27;      // 1992-07-01
    if (utc_jd < 2449534.5) return 28;      // 1993-07-01
    if (utc_jd < 2450083.5) return 29;      // 1994-07-01
    if (utc_jd < 2450630.5) return 30;      // 1996-01-01
    if (utc_jd < 2451179.5) return 31;      // 1997-07-01
    if (utc_jd < 2453736.5) return 32;      // 1999-01-01
    if (utc_jd < 2454832.5) return 33;      // 2006-01-01
    if (utc_jd < 2456109.5) return 34;      // 2009-01-01
    if (utc_jd < 2457204.5) return 35;      // 2012-07-01
    if (utc_jd < 2457754.5) return 36;      // 2015-07-01

    return 37;  // 2017-01-01 and after (update periodically)
}

double calendar_to_jd(int year, int month, int day,
                      int hour, int minute, double second) {
    // Algorithm from Meeus, "Astronomical Algorithms" (1998)

    // Adjust for January and February
    int a = (14 - month) / 12;
    int y = year + 4800 - a;
    int m = month + 12 * a - 3;

    // Julian day number at noon
    int jdn = day + (153 * m + 2) / 5 + 365 * y + y / 4 - y / 100 + y / 400 - 32045;

    // Add time of day
    double day_fraction = (hour - 12.0) / 24.0 + minute / 1440.0 + second / 86400.0;

    return static_cast<double>(jdn) + day_fraction;
}

void jd_to_calendar(double jd, int& year, int& month, int& day,
                   int& hour, int& minute, double& second) {
    // Algorithm from Meeus, "Astronomical Algorithms" (1998)

    // Separate integer and fractional parts
    int jdn = static_cast<int>(jd + 0.5);
    double frac = jd + 0.5 - static_cast<double>(jdn);

    int a = jdn + 32044;
    int b = (4 * a + 3) / 146097;
    int c = a - (146097 * b) / 4;
    int d = (4 * c + 3) / 1461;
    int e = c - (1461 * d) / 4;
    int m = (5 * e + 2) / 153;

    day = e - (153 * m + 2) / 5 + 1;
    month = m + 3 - 12 * (m / 10);
    year = 100 * b + d - 4800 + m / 10;

    // Extract time
    double day_seconds = frac * 86400.0;
    hour = static_cast<int>(day_seconds / 3600.0);
    day_seconds -= hour * 3600.0;
    minute = static_cast<int>(day_seconds / 60.0);
    second = day_seconds - minute * 60.0;
}

int days_in_month(int year, int month) {
    static const int days[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};

    if (month < 1 || month > 12) {
        throw std::out_of_range("Month must be between 1 and 12");
    }

    int d = days[month - 1];
    if (month == 2 && is_leap_year(year)) {
        d = 29;
    }

    return d;
}

double gmst(double ut1_jd) {
    // IAU 2006 GMST formula
    // Simplified version - for high precision use IAU SOFA library

    double T = (ut1_jd - J2000_JD) / 36525.0; // Julian centuries since J2000

    // GMST at 0h UT1 (in seconds)
    double gmst_0h = 24110.54841 + 8640184.812866 * T + 0.093104 * T * T - 6.2e-6 * T * T * T;

    // Add Earth rotation during current day
    double UT = (ut1_jd - std::floor(ut1_jd - 0.5) - 0.5) * 86400.0; // Seconds since 0h UT1
    double gmst_seconds = gmst_0h + 1.00273790935 * UT;

    // Convert to radians and normalize
    double gmst_rad = std::fmod(gmst_seconds, 86400.0) * (constants::TWO_PI / 86400.0);

    // Normalize to [0, 2π)
    while (gmst_rad < 0.0) gmst_rad += constants::TWO_PI;
    while (gmst_rad >= constants::TWO_PI) gmst_rad -= constants::TWO_PI;

    return gmst_rad;
}

double gast(double ut1_jd) {
    // For simplified implementation, GAST ≈ GMST
    // Full implementation requires nutation calculation (equation of equinoxes)

    double gmst_rad = gmst(ut1_jd);

    // Equation of equinoxes (simplified)
    double T = (ut1_jd - J2000_JD) / 36525.0;
    double omega = 125.04 - 1934.136 * T; // Longitude of ascending node of Moon
    double omega_rad = omega * constants::DEG_TO_RAD;

    double eq_equinox = 0.00264 * std::sin(omega_rad) + 0.000063 * std::sin(2.0 * omega_rad);
    double eq_equinox_rad = eq_equinox * constants::ARCSEC_TO_RAD;

    return gmst_rad + eq_equinox_rad;
}

} // namespace physim
