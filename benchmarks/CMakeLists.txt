# Benchmarks for physics simulation engine
# Phase 1: Placeholder - benchmarks will be added in Phase 2+

find_package(benchmark REQUIRED)

# Benchmark executables will be added as:
# add_executable(integrator_perf integrator_perf.cpp)
# target_link_libraries(integrator_perf PRIVATE physim benchmark::benchmark)

message(STATUS "Benchmarks configuration complete (benchmarks will be added in Phase 2)")
