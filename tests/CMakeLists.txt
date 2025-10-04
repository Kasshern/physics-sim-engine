# Test suite for physics simulation engine

find_package(GTest REQUIRED)

enable_testing()
include(GoogleTest)

# Unit tests subdirectory
add_subdirectory(unit)

# Integration tests subdirectory
# add_subdirectory(integration)

# Validation tests subdirectory
# add_subdirectory(validation)

message(STATUS "Test suite configuration complete")
