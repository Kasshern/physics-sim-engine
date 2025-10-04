# Unit tests for physics simulation engine

# Test: Integrators
add_executable(test_integrators test_integrators.cpp)
target_link_libraries(test_integrators PRIVATE physim GTest::gtest GTest::gtest_main)
gtest_discover_tests(test_integrators)

# Future tests:
# add_executable(test_forces test_forces.cpp)
# target_link_libraries(test_forces PRIVATE physim GTest::gtest GTest::gtest_main)
# gtest_discover_tests(test_forces)
