/**
 * @file test_nbody.cpp
 * @brief Unit tests for N-body system
 *
 * Tests particle and N-body system functionality.
 */

#include <gtest/gtest.h>
#include "physim/nbody/particle.hpp"
#include "physim/nbody/nbody_system.hpp"
#include "physim/integrators/rk4.hpp"
#include "physim/core/constants.hpp"
#include <cmath>

using namespace physim;
// Note: Don't use "using namespace physim::nbody" to avoid ambiguity with physim::Particle

class NBodyTest : public ::testing::Test {
protected:
    static constexpr double TOLERANCE = 1e-10;
};

// ============================================================================
// Particle Tests
// ============================================================================

TEST_F(NBodyTest, Particle_Construction) {
    Vec3 pos(1000.0, 2000.0, 3000.0);
    Vec3 vel(10.0, 20.0, 30.0);
    double mass = 5000.0;

    nbody::Particle p(pos, vel, mass, "TestParticle");

    EXPECT_EQ(p.position(), pos);
    EXPECT_EQ(p.velocity(), vel);
    EXPECT_EQ(p.mass(), mass);
    EXPECT_EQ(p.name(), "TestParticle");
}

TEST_F(NBodyTest, Particle_KineticEnergy) {
    Vec3 vel(1000.0, 0, 0);
    double mass = 2.0;
    nbody::Particle p(Vec3::Zero(), vel, mass);

    double KE = p.kinetic_energy();
    double expected = 0.5 * mass * vel.squaredNorm();

    EXPECT_NEAR(KE, expected, TOLERANCE);
}

TEST_F(NBodyTest, Particle_Momentum) {
    Vec3 vel(100.0, 200.0, 300.0);
    double mass = 5.0;
    nbody::Particle p(Vec3::Zero(), vel, mass);

    Vec3 momentum = p.momentum();
    Vec3 expected = mass * vel;

    EXPECT_NEAR((momentum - expected).norm(), 0.0, TOLERANCE);
}

TEST_F(NBodyTest, Particle_AngularMomentum) {
    Vec3 pos(1000.0, 0, 0);
    Vec3 vel(0, 500.0, 0);
    double mass = 2.0;
    nbody::Particle p(pos, vel, mass);

    Vec3 L = p.angular_momentum();
    Vec3 expected = pos.cross(mass * vel);

    EXPECT_NEAR((L - expected).norm(), 0.0, TOLERANCE);
}

TEST_F(NBodyTest, Particle_SolarSystemFactory) {
    auto earth = nbody::create_solar_system_particle("Earth", Vec3::Zero(), Vec3::Zero());

    EXPECT_NE(earth, nullptr);
    EXPECT_EQ(earth->name(), "Earth");
    EXPECT_GT(earth->mass(), 0.0);
}

// ============================================================================
// NBodySystem Tests
// ============================================================================

TEST_F(NBodyTest, NBodySystem_AddRemoveParticles) {
    nbody::NBodySystem system;

    EXPECT_TRUE(system.empty());
    EXPECT_EQ(system.size(), 0);

    nbody::Particle p1(Vec3(1000, 0, 0), Vec3::Zero(), 1000.0);
    system.add_particle(p1);

    EXPECT_FALSE(system.empty());
    EXPECT_EQ(system.size(), 1);

    system.clear();
    EXPECT_TRUE(system.empty());
}

TEST_F(NBodyTest, NBodySystem_TwoBodyForces) {
    nbody::NBodySystem system;

    const double m1 = 1e24;  // kg
    const double m2 = 1e20;  // kg
    const double r = 1e8;    // m

    nbody::Particle p1(Vec3(-r/2, 0, 0), Vec3::Zero(), m1);
    nbody::Particle p2(Vec3(r/2, 0, 0), Vec3::Zero(), m2);

    system.add_particle(p1);
    system.add_particle(p2);

    system.compute_accelerations();

    const auto& p1_f = system.particle(0);
    const auto& p2_f = system.particle(1);

    // Check acceleration directions (should point toward each other)
    EXPECT_GT(p1_f.acceleration()(0), 0.0) << "p1 should accelerate toward p2 (+x)";
    EXPECT_LT(p2_f.acceleration()(0), 0.0) << "p2 should accelerate toward p1 (-x)";

    // Check Newton's third law: F1 = -F2 => a1/a2 = -m2/m1
    const double ratio = p1_f.acceleration().norm() / p2_f.acceleration().norm();
    const double expected_ratio = m2 / m1;

    EXPECT_NEAR(ratio, expected_ratio, 1e-6) << "Newton's third law";
}

TEST_F(NBodyTest, NBodySystem_CenterOfMass) {
    nbody::NBodySystem system;

    // Two equal masses at symmetric positions
    const double mass = 1000.0;
    nbody::Particle p1(Vec3(1000, 0, 0), Vec3::Zero(), mass);
    nbody::Particle p2(Vec3(-1000, 0, 0), Vec3::Zero(), mass);

    system.add_particle(p1);
    system.add_particle(p2);

    Vec3 com = system.center_of_mass();

    // COM should be at origin
    EXPECT_NEAR(com.norm(), 0.0, TOLERANCE) << "COM of symmetric system";
}

TEST_F(NBodyTest, NBodySystem_TotalMomentum) {
    nbody::NBodySystem system;

    nbody::Particle p1(Vec3::Zero(), Vec3(100, 0, 0), 1000.0);
    nbody::Particle p2(Vec3::Zero(), Vec3(-50, 0, 0), 2000.0);

    system.add_particle(p1);
    system.add_particle(p2);

    Vec3 total_p = system.total_momentum();
    Vec3 expected = p1.momentum() + p2.momentum();

    EXPECT_NEAR((total_p - expected).norm(), 0.0, TOLERANCE);
}

TEST_F(NBodyTest, NBodySystem_EnergyTracking) {
    nbody::NBodySystem system;

    const double m = 1e24;
    const double r = 1e8;
    const double v = 1e4;

    nbody::Particle p1(Vec3(-r, 0, 0), Vec3(0, v, 0), m);
    nbody::Particle p2(Vec3(r, 0, 0), Vec3(0, -v, 0), m);

    system.add_particle(p1);
    system.add_particle(p2);

    system.reset_energy_tracking();

    auto stats = system.compute_stats();

    EXPECT_EQ(stats.num_particles, 2);
    EXPECT_GT(stats.kinetic_energy, 0.0);
    EXPECT_LT(stats.potential_energy, 0.0);
    EXPECT_NEAR(stats.energy_error, 0.0, TOLERANCE) << "Just reset tracking";
}

TEST_F(NBodyTest, NBodySystem_CircularOrbitIntegration) {
    nbody::NBodySystem system;

    // Central mass
    const double M = constants::gm::SUN / constants::GRAVITATIONAL_CONSTANT;
    nbody::Particle sun(Vec3::Zero(), Vec3::Zero(), M, "Sun");
    system.add_particle(sun);

    // Orbiting body in circular orbit
    const double r = 1e11;  // 100 million km
    const double v_circ = std::sqrt(constants::gm::SUN / r);
    nbody::Particle planet(Vec3(r, 0, 0), Vec3(0, v_circ, 0), 1e20, "Planet");
    system.add_particle(planet);

    system.reset_energy_tracking();

    // Integrate for a short time
    integrators::RK4 rk4;
    const double dt = 3600.0;   // 1 hour
    const double t_final = 86400.0;  // 1 day

    system.propagate(t_final, dt, rk4);

    auto stats_f = system.compute_stats();

    // Energy should be conserved
    EXPECT_LT(stats_f.energy_error, 1e-6) << "Energy conservation in circular orbit";

    // Orbit radius should be approximately constant
    const auto& planet_f = system.particle(1);
    const double r_final = planet_f.position().norm();
    const double r_error = std::abs(r_final - r) / r;

    EXPECT_LT(r_error, 0.01) << "Orbit radius variation < 1%";
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
