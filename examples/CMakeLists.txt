# Examples for physics simulation engine
# Phase 1: Foundation test
# Phase 2+: Full examples will be added

# Phase 1 foundation test
add_executable(test_phase1 test_phase1.cpp)
target_link_libraries(test_phase1 PRIVATE physim)

# Example executables will be added as:
# add_executable(01_solar_system 01_solar_system.cpp)
# target_link_libraries(01_solar_system PRIVATE physim)

message(STATUS "Examples configuration complete")
