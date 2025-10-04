/**
 * @file test_integrators.cpp
 * @brief Unit tests for ODE integrators
 *
 * Tests integrator convergence, order of accuracy, and error control.
 */

#include <gtest/gtest.h>
#include "physim/integrators/rk4.hpp"
#include "physim/integrators/rk45.hpp"
#include <cmath>

using namespace physim;
using namespace physim::integrators;

// ============================================================================
// Test Fixtures
// ============================================================================

class IntegratorTest : public ::testing::Test {
protected:
    // Simple harmonic oscillator: y'' = -ω²y
    // Analytical solution: y(t) = A*cos(ωt) + B*sin(ωt)
    static DerivativeFunction harmonic_oscillator(double omega) {
        return [omega](double t, const VecX& y) -> VecX {
            (void)t;
            VecX dydt(2);
            dydt(0) = y(1);           // dy/dt = v
            dydt(1) = -omega * omega * y(0);  // dv/dt = -ω²y
            return dydt;
        };
    }

    // Exponential decay: y' = -λy
    // Analytical solution: y(t) = y0 * exp(-λt)
    static DerivativeFunction exponential_decay(double lambda) {
        return [lambda](double t, const VecX& y) -> VecX {
            (void)t;
            VecX dydt(1);
            dydt(0) = -lambda * y(0);
            return dydt;
        };
    }
};

// ============================================================================
// RK4 Tests
// ============================================================================

TEST_F(IntegratorTest, RK4_HarmonicOscillator) {
    RK4 integrator;

    const double omega = 1.0;
    auto f = harmonic_oscillator(omega);

    // Initial conditions: y(0) = 1, v(0) = 0
    // Analytical: y(t) = cos(ωt)
    VecX y0(2);
    y0 << 1.0, 0.0;

    const double t0 = 0.0;
    const double tf = 2.0 * M_PI;  // One period
    const double dt = 0.01;

    VecX yf = integrator.integrate(t0, y0, tf, dt, f);

    // After one period, should return to initial state
    const double error_y = std::abs(yf(0) - 1.0);
    const double error_v = std::abs(yf(1) - 0.0);

    EXPECT_LT(error_y, 1e-6) << "Position error after one period";
    EXPECT_LT(error_v, 1e-6) << "Velocity error after one period";

    // Energy should be conserved (E = v²/2 + ω²y²/2)
    const double E0 = 0.5 * (y0(1) * y0(1) + omega * omega * y0(0) * y0(0));
    const double Ef = 0.5 * (yf(1) * yf(1) + omega * omega * yf(0) * yf(0));
    const double energy_error = std::abs((Ef - E0) / E0);

    EXPECT_LT(energy_error, 1e-6) << "Energy conservation";
}

TEST_F(IntegratorTest, RK4_ExponentialDecay) {
    RK4 integrator;

    const double lambda = 0.5;
    auto f = exponential_decay(lambda);

    VecX y0(1);
    y0 << 1.0;

    const double t0 = 0.0;
    const double tf = 2.0;
    const double dt = 0.01;

    VecX yf = integrator.integrate(t0, y0, tf, dt, f);

    // Analytical solution: y(2) = exp(-0.5 * 2) = exp(-1)
    const double y_analytical = std::exp(-lambda * tf);
    const double error = std::abs(yf(0) - y_analytical);

    EXPECT_LT(error, 1e-6) << "Exponential decay accuracy";
}

TEST_F(IntegratorTest, RK4_ConvergenceOrder) {
    // Test that RK4 has 4th order convergence: error ∝ h⁴
    RK4 integrator;

    const double omega = 1.0;
    auto f = harmonic_oscillator(omega);

    VecX y0(2);
    y0 << 1.0, 0.0;

    const double t0 = 0.0;
    const double tf = 1.0;

    // Analytical solution at tf: y(1) = cos(1)
    const double y_exact = std::cos(omega * tf);

    // Test with different step sizes
    std::vector<double> step_sizes = {0.1, 0.05, 0.025};
    std::vector<double> errors;

    for (double dt : step_sizes) {
        VecX yf = integrator.integrate(t0, y0, tf, dt, f);
        errors.push_back(std::abs(yf(0) - y_exact));
    }

    // Check convergence order: error(h/2) / error(h) ≈ 2⁴ = 16
    for (size_t i = 0; i < errors.size() - 1; ++i) {
        double ratio = errors[i] / errors[i + 1];
        EXPECT_NEAR(ratio, 16.0, 2.0) << "4th order convergence at step " << i;
    }
}

// ============================================================================
// RK45 Tests
// ============================================================================

TEST_F(IntegratorTest, RK45_AdaptiveStep) {
    RK45 integrator;

    const double omega = 1.0;
    auto f = harmonic_oscillator(omega);

    VecX y0(2);
    y0 << 1.0, 0.0;

    const double t0 = 0.0;
    const double tf = 2.0 * M_PI;
    const double dt = 0.1;
    const double tol = 1e-12;

    VecX yf = integrator.integrate(t0, y0, tf, dt, f, true, tol);

    const auto& stats = integrator.stats();

    // Check that adaptive stepping worked
    EXPECT_GT(stats.num_steps, 0);
    EXPECT_LT(stats.num_steps, tf / 0.01) << "Should use larger steps than fixed dt=0.01";

    // Check accuracy
    const double error_y = std::abs(yf(0) - 1.0);
    const double error_v = std::abs(yf(1) - 0.0);

    EXPECT_LT(error_y, 1e-8) << "RK45 adaptive position accuracy";
    EXPECT_LT(error_v, 1e-8) << "RK45 adaptive velocity accuracy";

    // Max error should be within tolerance
    EXPECT_LT(stats.max_error, 10.0) << "Max error within bounds";
}

TEST_F(IntegratorTest, RK45_ErrorControl) {
    RK45 integrator;

    const double lambda = 1.0;
    auto f = exponential_decay(lambda);

    VecX y0(1);
    y0 << 1.0;

    const double t0 = 0.0;
    const double tf = 2.0;
    const double dt = 1.0;

    // Test with different tolerances
    std::vector<double> tolerances = {1e-6, 1e-9, 1e-12};
    std::vector<double> errors;

    for (double tol : tolerances) {
        integrator.reset_stats();
        VecX yf = integrator.integrate(t0, y0, tf, dt, f, true, tol);

        const double y_exact = std::exp(-lambda * tf);
        errors.push_back(std::abs(yf(0) - y_exact));
    }

    // Tighter tolerance should give smaller errors
    for (size_t i = 0; i < errors.size() - 1; ++i) {
        EXPECT_LT(errors[i + 1], errors[i]) <<
            "Tighter tolerance should improve accuracy";
    }
}

TEST_F(IntegratorTest, RK45_FixedStepFallback) {
    // Test that RK45 works in fixed-step mode (adaptive=false)
    RK45 integrator;

    const double omega = 1.0;
    auto f = harmonic_oscillator(omega);

    VecX y0(2);
    y0 << 1.0, 0.0;

    const double t0 = 0.0;
    const double tf = 1.0;
    const double dt = 0.01;

    VecX yf = integrator.integrate(t0, y0, tf, dt, f, false);  // Fixed step

    const double y_exact = std::cos(omega * tf);
    const double error = std::abs(yf(0) - y_exact);

    EXPECT_LT(error, 1e-6) << "RK45 fixed-step mode accuracy";
}

// ============================================================================
// Edge Cases
// ============================================================================

TEST_F(IntegratorTest, ZeroStepSize) {
    RK4 integrator;

    VecX y0(1);
    y0 << 1.0;

    auto f = [](double t, const VecX& y) -> VecX {
        (void)t;
        return y;
    };

    // Should throw for dt <= 0
    EXPECT_THROW(integrator.integrate(0.0, y0, 1.0, 0.0, f), std::invalid_argument);
    EXPECT_THROW(integrator.integrate(0.0, y0, 1.0, -0.1, f), std::invalid_argument);
}

TEST_F(IntegratorTest, InvalidTimeRange) {
    RK4 integrator;

    VecX y0(1);
    y0 << 1.0;

    auto f = [](double t, const VecX& y) -> VecX {
        (void)t;
        return y;
    };

    // Should throw for tf <= t0
    EXPECT_THROW(integrator.integrate(1.0, y0, 0.0, 0.1, f), std::invalid_argument);
    EXPECT_THROW(integrator.integrate(1.0, y0, 1.0, 0.1, f), std::invalid_argument);
}

// ============================================================================
// Main
// ============================================================================

int main(int argc, char** argv) {
    ::testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}
