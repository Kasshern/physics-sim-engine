/**
 * @file test_forces.cpp
 * @brief Unit tests for force models
 *
 * Tests force accuracy, utility functions, and multi-force systems.
 */

#include <gtest/gtest.h>
#include "physim/forces/force.hpp"
#include "physim/forces/gravity.hpp"
#include "physim/forces/j2_gravity.hpp"
#include "physim/core/constants.hpp"
#include <cmath>

using namespace physim;
using namespace physim::forces;

// ============================================================================
// Test Fixtures
// ============================================================================

class ForceTest : public ::testing::Test {
protected:
    static constexpr double TOLERANCE = 1e-10;
};

// ============================================================================
// Point Mass Gravity Tests
// ============================================================================

TEST_F(ForceTest, PointMassGravity_Acceleration) {
    // Create Earth gravity
    auto earth = PointMassGravity::for_body("Earth");

    // Test at Earth's surface (r = R_earth)
    const double r_earth = 6371e3;  // m
    Vec3 pos(r_earth, 0, 0);
    Vec3 vel = Vec3::Zero();

    Vec3 accel = earth->acceleration(0.0, pos, vel, 1.0);

    // Expected: g = GM/r² pointing toward origin
    const double g_expected = constants::gm::EARTH / (r_earth * r_earth);
    const double g_magnitude = accel.norm();

    EXPECT_NEAR(g_magnitude, g_expected, TOLERANCE) << "Gravity magnitude at surface";

    // Direction should point toward origin
    Vec3 accel_dir = accel.normalized();
    Vec3 expected_dir(-1.0, 0.0, 0.0);
    EXPECT_NEAR(accel_dir.dot(expected_dir), 1.0, 1e-12) << "Gravity direction";
}

TEST_F(ForceTest, PointMassGravity_InverseSquareLaw) {
    auto earth = PointMassGravity::for_body("Earth");

    // Test at different radii
    const double r1 = 10000e3;  // 10,000 km
    const double r2 = 20000e3;  // 20,000 km

    Vec3 pos1(r1, 0, 0);
    Vec3 pos2(r2, 0, 0);
    Vec3 vel = Vec3::Zero();

    Vec3 accel1 = earth->acceleration(0.0, pos1, vel, 1.0);
    Vec3 accel2 = earth->acceleration(0.0, pos2, vel, 1.0);

    // Doubling distance should quarter the acceleration
    const double ratio = accel1.norm() / accel2.norm();
    EXPECT_NEAR(ratio, 4.0, 1e-10) << "Inverse square law";
}

TEST_F(ForceTest, PointMassGravity_CircularVelocity) {
    auto earth = PointMassGravity::for_body("Earth");

    const double r = 7000e3;  // 7,000 km (LEO)
    const double v_circ = earth->circular_velocity(r);

    // For circular orbit: v = sqrt(μ/r)
    const double v_expected = std::sqrt(constants::gm::EARTH / r);

    EXPECT_NEAR(v_circ, v_expected, TOLERANCE) << "Circular velocity";
}

TEST_F(ForceTest, PointMassGravity_EscapeVelocity) {
    auto earth = PointMassGravity::for_body("Earth");

    const double r = 7000e3;
    const double v_esc = earth->escape_velocity(r);

    // v_esc = sqrt(2μ/r) = sqrt(2) * v_circ
    const double v_circ = earth->circular_velocity(r);

    EXPECT_NEAR(v_esc, v_circ * std::sqrt(2.0), TOLERANCE) << "Escape velocity";
}

TEST_F(ForceTest, PointMassGravity_OrbitalPeriod) {
    auto earth = PointMassGravity::for_body("Earth");

    const double r = 7000e3;
    const double period = earth->orbital_period(r);

    // Kepler's third law: T = 2π * sqrt(r³/μ)
    const double period_expected = 2.0 * constants::PI *
                                   std::sqrt(r * r * r / constants::gm::EARTH);

    EXPECT_NEAR(period, period_expected, TOLERANCE) << "Orbital period";
}

TEST_F(ForceTest, PointMassGravity_PotentialEnergy) {
    auto earth = PointMassGravity::for_body("Earth");

    const double r = 10000e3;
    const double mass = 1000.0;  // 1000 kg
    Vec3 pos(r, 0, 0);

    const double U = earth->potential_energy(0.0, pos, mass);

    // U = -GMm/r
    const double U_expected = -constants::gm::EARTH * mass / r;

    EXPECT_NEAR(U, U_expected, 1e-6) << "Potential energy";
}

TEST_F(ForceTest, PointMassGravity_CenterOffset) {
    // Create gravity with non-zero center position
    Vec3 center(1000e3, 0, 0);
    auto gravity = std::make_shared<PointMassGravity>(constants::gm::EARTH, center);

    // Position relative to origin
    Vec3 pos(8000e3, 0, 0);
    Vec3 vel = Vec3::Zero();

    Vec3 accel = gravity->acceleration(0.0, pos, vel, 1.0);

    // Relative position: pos - center
    Vec3 r_rel = pos - center;
    const double r = r_rel.norm();

    // Expected acceleration magnitude
    const double g_expected = constants::gm::EARTH / (r * r);

    EXPECT_NEAR(accel.norm(), g_expected, TOLERANCE) << "Offset center magnitude";
}

TEST_F(ForceTest, PointMassGravity_NamedBodies) {
    // Test all named body constructors
    std::vector<std::string> bodies = {"Sun", "Earth", "Moon", "Mars", "Jupiter"};

    for (const auto& body : bodies) {
        auto gravity = PointMassGravity::for_body(body);
        EXPECT_NE(gravity, nullptr) << "Named body: " << body;

        // Check that GM is positive
        Vec3 pos(1e6, 0, 0);
        Vec3 vel = Vec3::Zero();
        Vec3 accel = gravity->acceleration(0.0, pos, vel, 1.0);
        EXPECT_GT(accel.norm(), 0.0) << "Non-zero gravity for " << body;
    }
}

TEST_F(ForceTest, PointMassGravity_SingularityDetection) {
    auto earth = PointMassGravity::for_body("Earth");

    // Test at origin (singularity)
    Vec3 pos_zero = Vec3::Zero();
    Vec3 vel = Vec3::Zero();

    // Should throw at r = 0
    EXPECT_THROW(earth->acceleration(0.0, pos_zero, vel, 1.0), std::runtime_error);
}

// ============================================================================
// J2 Gravity Tests
// ============================================================================

TEST_F(ForceTest, J2Gravity_PerturbationMagnitude) {
    auto j2 = std::make_shared<J2Gravity>();

    // LEO orbit
    const double r = 7000e3;
    Vec3 pos(r, 0, 0);
    Vec3 vel = Vec3::Zero();

    Vec3 accel_j2 = j2->acceleration(0.0, pos, vel, 1.0);

    // J2 acceleration should be much smaller than point mass
    auto point_mass = PointMassGravity::for_body("Earth");
    Vec3 accel_pm = point_mass->acceleration(0.0, pos, vel, 1.0);

    const double ratio = accel_j2.norm() / accel_pm.norm();

    // J2 is ~1000x smaller than point mass
    EXPECT_LT(ratio, 0.01) << "J2 perturbation is small";
    EXPECT_GT(ratio, 0.0) << "J2 perturbation is non-zero";
}

TEST_F(ForceTest, J2Gravity_NodalRegression) {
    auto j2 = std::make_shared<J2Gravity>();

    // Sun-synchronous orbit parameters
    const double a = 7178e3;  // ~800 km altitude
    const double i = 98.0 * constants::DEG_TO_RAD;

    const double omega_dot = j2->nodal_regression_rate(a, i);

    // For sun-synchronous: Ω̇ ≈ 0.9856 deg/day = 1.991e-7 rad/s
    const double omega_dot_sso = 0.9856 * constants::DEG_TO_RAD / constants::SECONDS_PER_DAY;

    EXPECT_NEAR(omega_dot, omega_dot_sso, 2e-8) << "Sun-synchronous nodal regression";
}

TEST_F(ForceTest, J2Gravity_ApsidialPrecession) {
    auto j2 = std::make_shared<J2Gravity>();

    // Circular orbit
    const double a = 7000e3;
    const double e = 0.0;
    const double i = 45.0 * constants::DEG_TO_RAD;

    const double omega_bar_dot = j2->apsidal_precession_rate(a, e, i);

    // For circular orbit, precession should be non-zero
    EXPECT_NE(omega_bar_dot, 0.0) << "Apsidal precession for circular orbit";
}

TEST_F(ForceTest, J2Gravity_SunSynchronousInclination) {
    auto j2 = std::make_shared<J2Gravity>();

    const double a = 7178e3;  // ~800 km altitude
    const double i_sso = j2->sun_synchronous_inclination(a);

    // Sun-synchronous inclination ~98 degrees
    const double i_expected = 98.0 * constants::DEG_TO_RAD;

    EXPECT_NEAR(i_sso, i_expected, 1.0 * constants::DEG_TO_RAD) << "SSO inclination";
}

TEST_F(ForceTest, J2Gravity_EquatorialVsPolar) {
    auto j2 = std::make_shared<J2Gravity>();

    const double r = 7000e3;

    // Equatorial position (z = 0)
    Vec3 pos_eq(r, 0, 0);
    Vec3 vel = Vec3::Zero();
    Vec3 accel_eq = j2->acceleration(0.0, pos_eq, vel, 1.0);

    // Polar position (x = y = 0)
    Vec3 pos_polar(0, 0, r);
    Vec3 accel_polar = j2->acceleration(0.0, pos_polar, vel, 1.0);

    // J2 effect is different at equator vs poles
    EXPECT_NE(accel_eq.norm(), accel_polar.norm()) << "Latitude dependence";
}

// ============================================================================
// ForceModel (Multi-Force) Tests
// ============================================================================

TEST_F(ForceTest, ForceModel_SingleForce) {
    ForceModel model;
    model.add_force(PointMassGravity::for_body("Earth"));

    const double r = 7000e3;
    Vec3 pos(r, 0, 0);
    Vec3 vel = Vec3::Zero();

    Vec3 accel = model.total_acceleration(0.0, pos, vel, 1.0);

    // Should match standalone force
    auto earth = PointMassGravity::for_body("Earth");
    Vec3 accel_expected = earth->acceleration(0.0, pos, vel, 1.0);

    EXPECT_NEAR((accel - accel_expected).norm(), 0.0, TOLERANCE) << "Single force";
}

TEST_F(ForceTest, ForceModel_MultipleForces) {
    ForceModel model;
    model.add_force(PointMassGravity::for_body("Earth"));
    model.add_force(std::make_shared<J2Gravity>());

    const double r = 7000e3;
    Vec3 pos(r, 0, 0);
    Vec3 vel = Vec3::Zero();

    Vec3 accel_total = model.total_acceleration(0.0, pos, vel, 1.0);

    // Should be sum of individual forces
    auto earth = PointMassGravity::for_body("Earth");
    auto j2 = std::make_shared<J2Gravity>();

    Vec3 accel_earth = earth->acceleration(0.0, pos, vel, 1.0);
    Vec3 accel_j2 = j2->acceleration(0.0, pos, vel, 1.0);
    Vec3 accel_expected = accel_earth + accel_j2;

    EXPECT_NEAR((accel_total - accel_expected).norm(), 0.0, TOLERANCE) << "Multiple forces sum";
}

TEST_F(ForceTest, ForceModel_EnableDisable) {
    ForceModel model;
    auto earth = PointMassGravity::for_body("Earth");
    auto j2 = std::make_shared<J2Gravity>();

    model.add_force(earth);
    model.add_force(j2);

    const double r = 7000e3;
    Vec3 pos(r, 0, 0);
    Vec3 vel = Vec3::Zero();

    // Both forces enabled
    Vec3 accel_both = model.total_acceleration(0.0, pos, vel, 1.0);

    // Disable J2
    j2->set_enabled(false);
    Vec3 accel_earth_only = model.total_acceleration(0.0, pos, vel, 1.0);

    // Re-enable J2
    j2->set_enabled(true);
    Vec3 accel_both_again = model.total_acceleration(0.0, pos, vel, 1.0);

    EXPECT_NEAR((accel_both - accel_both_again).norm(), 0.0, TOLERANCE) << "Re-enable force";
    EXPECT_LT(accel_earth_only.norm(), accel_both.norm()) << "Disabled force reduces total";
}

TEST_F(ForceTest, ForceModel_TotalPotentialEnergy) {
    ForceModel model;
    model.add_force(PointMassGravity::for_body("Earth"));

    const double r = 10000e3;
    const double mass = 1000.0;
    Vec3 pos(r, 0, 0);

    const double U_total = model.total_potential_energy(0.0, pos, mass);

    // Should match standalone force
    auto earth = PointMassGravity::for_body("Earth");
    const double U_expected = earth->potential_energy(0.0, pos, mass);

    EXPECT_NEAR(U_total, U_expected, 1e-6) << "Total potential energy";
}

TEST_F(ForceTest, ForceModel_VelocityDependence) {
    ForceModel model;
    auto earth = PointMassGravity::for_body("Earth");

    model.add_force(earth);

    // Gravity is not velocity dependent
    EXPECT_FALSE(earth->is_velocity_dependent()) << "Gravity not velocity dependent";
}

// ============================================================================
// Energy Conservation Tests
// ============================================================================

TEST_F(ForceTest, EnergyConservation_CircularOrbit) {
    auto earth = PointMassGravity::for_body("Earth");

    // Circular orbit at r = 7000 km
    const double r = 7000e3;
    const double v_circ = earth->circular_velocity(r);

    Vec3 pos(r, 0, 0);
    Vec3 vel(0, v_circ, 0);
    const double mass = 1.0;

    // Total energy: E = v²/2 - μ/r
    const double KE = 0.5 * mass * vel.squaredNorm();
    const double PE = earth->potential_energy(0.0, pos, mass);
    const double E_total = KE + PE;

    // For circular orbit: E = -μ/(2r)
    const double E_expected = -constants::gm::EARTH / (2.0 * r);

    EXPECT_NEAR(E_total, E_expected, 1e-6) << "Circular orbit energy";
}

TEST_F(ForceTest, EnergyConservation_EscapeTrajectory) {
    auto earth = PointMassGravity::for_body("Earth");

    const double r = 7000e3;
    const double v_esc = earth->escape_velocity(r);

    Vec3 pos(r, 0, 0);
    Vec3 vel(0, v_esc, 0);
    const double mass = 1.0;

    const double KE = 0.5 * mass * vel.squaredNorm();
    const double PE = earth->potential_energy(0.0, pos, mass);
    const double E_total = KE + PE;

    // For escape velocity: E = 0
    EXPECT_NEAR(E_total, 0.0, 1e-6) << "Escape trajectory energy";
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
