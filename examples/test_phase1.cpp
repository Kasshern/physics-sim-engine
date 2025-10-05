/**
 * @file test_phase1.cpp
 * @brief Simple test to verify Phase 1 foundation compiles
 */

#include "physim/core/types.hpp"
#include "physim/core/constants.hpp"
#include "physim/core/time.hpp"
#include "physim/core/frame.hpp"
#include "physim/core/logging.hpp"

#include <iostream>

int main() {
    using namespace physim;

    // Initialize logging
    logging::init(spdlog::level::info, false);
    PHYSIM_LOG_INFO("Phase 1 foundation test starting");

    // Test core types
    Vec3 position(1.0, 2.0, 3.0);
    Vec3 velocity(-0.5, 0.5, 0.0);
    State state(position, velocity, 1000.0);

    std::cout << "State created:\n";
    std::cout << "  Position: [" << state.position.transpose() << "]\n";
    std::cout << "  Velocity: [" << state.velocity.transpose() << "]\n";
    std::cout << "  Mass: " << state.mass << " kg\n";

    // Test constants
    std::cout << "\nPhysical constants:\n";
    std::cout << "  G = " << constants::GRAVITATIONAL_CONSTANT << " m^3/(kg·s^2)\n";
    std::cout << "  Earth GM = " << constants::gm::EARTH << " m^3/s^2\n";
    std::cout << "  Earth radius = " << constants::radius::EARTH << " m\n";
    std::cout << "  Speed of light = " << constants::SPEED_OF_LIGHT << " m/s\n";

    // Test time system
    Time j2000 = Time(J2000_JD, TimeScale::TT);
    std::cout << "\nTime system:\n";
    std::cout << "  J2000.0 JD = " << j2000.jd() << "\n";
    std::cout << "  J2000.0 MJD = " << j2000.mjd() << "\n";
    std::cout << "  Seconds since J2000 = " << j2000.seconds_since_j2000() << "\n";

    Time current = Time::from_calendar(2025, 10, 4, 12, 0, 0.0, TimeScale::UTC);
    std::cout << "  Current time: " << current.to_iso_string() << " UTC\n";

    // Test coordinate transformations
    Vec3 r_eci(constants::earth::WGS84_A + 400e3, 0, 0); // 400 km altitude
    Vec3 v_eci(0, 7660, 0); // Roughly circular orbit velocity

    std::cout << "\nCoordinate transformations:\n";
    std::cout << "  ECI position: [" << r_eci.transpose() << "] m\n";
    std::cout << "  ECI velocity: [" << v_eci.transpose() << "] m/s\n";

    // Test LVLH frame
    Mat3 R_lvlh = lvlh_frame(r_eci, v_eci);
    std::cout << "  LVLH frame computed successfully\n";

    // Test geodetic conversion
    Vec3 r_ecef = r_eci; // Approximate for this test
    GeodeticCoord geodetic = ecef_to_geodetic(r_ecef);
    std::cout << "  Geodetic latitude: " << constants::rad_to_deg(geodetic.latitude) << " deg\n";
    std::cout << "  Geodetic longitude: " << constants::rad_to_deg(geodetic.longitude) << " deg\n";
    std::cout << "  Geodetic altitude: " << geodetic.altitude / 1000.0 << " km\n";

    // Test utility functions
    double ke = kinetic_energy(state);
    Vec3 momentum = linear_momentum(state);
    std::cout << "\nDynamic quantities:\n";
    std::cout << "  Kinetic energy: " << ke << " J\n";
    std::cout << "  Linear momentum: [" << momentum.transpose() << "] kg·m/s\n";

    PHYSIM_LOG_INFO("Phase 1 foundation test completed successfully");
    logging::shutdown();

    std::cout << "\n✓ Phase 1 foundation compiled and executed successfully!\n";

    return 0;
}
