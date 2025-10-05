# Examples for physics simulation engine

# Phase 1 foundation test
add_executable(test_phase1 test_phase1.cpp)
target_link_libraries(test_phase1 PRIVATE physim)

# Phase 2: Earth-Moon two-body problem
add_executable(earth_moon_orbit earth_moon_orbit.cpp)
target_link_libraries(earth_moon_orbit PRIVATE physim)

# Phase 3: N-body examples
add_executable(two_body_nbody two_body_nbody.cpp)
target_link_libraries(two_body_nbody PRIVATE physim)

add_executable(inner_solar_system inner_solar_system.cpp)
target_link_libraries(inner_solar_system PRIVATE physim)

message(STATUS "Examples configuration complete")
