# Examples for physics simulation engine

# Phase 1 foundation test
add_executable(test_phase1 test_phase1.cpp)
target_link_libraries(test_phase1 PRIVATE physim)

# Phase 2: Earth-Moon two-body problem
add_executable(earth_moon_orbit earth_moon_orbit.cpp)
target_link_libraries(earth_moon_orbit PRIVATE physim)

# Future examples:
# add_executable(01_solar_system 01_solar_system.cpp)
# target_link_libraries(01_solar_system PRIVATE physim)

message(STATUS "Examples configuration complete")
