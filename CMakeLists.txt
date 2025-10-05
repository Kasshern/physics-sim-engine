cmake_minimum_required(VERSION 3.25)
project(PhysicsSimEngine VERSION 0.1.0 LANGUAGES CXX)

# Set C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

# Project options
option(PHYSIM_BUILD_TESTS "Build tests" ON)
option(PHYSIM_BUILD_BENCHMARKS "Build benchmarks" OFF)
option(PHYSIM_BUILD_EXAMPLES "Build examples" ON)
option(PHYSIM_BUILD_PYTHON "Build Python bindings" OFF)
option(PHYSIM_ENABLE_CUDA "Enable CUDA GPU acceleration" OFF)
option(PHYSIM_ENABLE_VISUALIZATION "Build visualization support" OFF)
option(PHYSIM_USE_ASAN "Enable AddressSanitizer" OFF)
option(PHYSIM_USE_TSAN "Enable ThreadSanitizer" OFF)
option(PHYSIM_USE_UBSAN "Enable UndefinedBehaviorSanitizer" OFF)

# Enable CUDA if requested (will be added in Phase 5)
if(PHYSIM_ENABLE_CUDA)
    enable_language(CUDA)
    set(CMAKE_CUDA_STANDARD 17)
    set(CMAKE_CUDA_STANDARD_REQUIRED ON)
endif()

# Set default build type
if(NOT CMAKE_BUILD_TYPE AND NOT CMAKE_CONFIGURATION_TYPES)
    set(CMAKE_BUILD_TYPE "Release" CACHE STRING "Choose the type of build." FORCE)
    set_property(CACHE CMAKE_BUILD_TYPE PROPERTY STRINGS "Debug" "Release" "MinSizeRel" "RelWithDebInfo")
endif()

# Compiler-specific flags
if(CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang")
    add_compile_options(
        -Wall
        -Wextra
        -Wpedantic
        -Wconversion
        -Wshadow
        -Wnon-virtual-dtor
        -Wold-style-cast
        -Wcast-align
        -Wunused
        -Woverloaded-virtual
    )

    if(CMAKE_BUILD_TYPE MATCHES "Release")
        add_compile_options(-O3 -march=native -ffast-math -funroll-loops)
    endif()

    # Sanitizer support
    if(PHYSIM_USE_ASAN)
        add_compile_options(-fsanitize=address -fno-omit-frame-pointer -g)
        add_link_options(-fsanitize=address)
    endif()

    if(PHYSIM_USE_TSAN)
        add_compile_options(-fsanitize=thread -g)
        add_link_options(-fsanitize=thread)
    endif()

    if(PHYSIM_USE_UBSAN)
        add_compile_options(-fsanitize=undefined -g)
        add_link_options(-fsanitize=undefined)
    endif()
elseif(MSVC)
    add_compile_options(/W4 /permissive-)
    if(CMAKE_BUILD_TYPE MATCHES "Release")
        add_compile_options(/O2 /arch:AVX2)
    endif()
endif()

# Find required packages
find_package(Eigen3 3.4 QUIET NO_MODULE)
if(NOT Eigen3_FOUND)
    message(WARNING "Eigen3 not found. Please install via: brew install eigen or vcpkg install eigen3")
endif()

find_package(Boost 1.70 QUIET)
if(NOT Boost_FOUND)
    message(WARNING "Boost not found. Please install via: brew install boost or vcpkg install boost")
endif()

find_package(spdlog QUIET)
if(NOT spdlog_FOUND)
    message(WARNING "spdlog not found. Please install via: brew install spdlog or vcpkg install spdlog")
endif()

find_package(fmt QUIET)
if(NOT fmt_FOUND)
    message(WARNING "fmt not found. Please install via: brew install fmt or vcpkg install fmt")
endif()

# Optional dependencies based on build options
if(PHYSIM_ENABLE_CUDA)
    find_package(CUDAToolkit REQUIRED)
endif()

if(PHYSIM_ENABLE_VISUALIZATION)
    find_package(glfw3 REQUIRED)
    find_package(imgui REQUIRED)
endif()

# Create main library target
add_library(physim SHARED)

target_include_directories(physim
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
)

target_link_libraries(physim
    PUBLIC
        $<$<TARGET_EXISTS:Eigen3::Eigen>:Eigen3::Eigen>
        $<$<TARGET_EXISTS:Boost::boost>:Boost::boost>
        $<$<TARGET_EXISTS:spdlog::spdlog>:spdlog::spdlog>
        $<$<TARGET_EXISTS:fmt::fmt>:fmt::fmt>
)

# Compiler feature requirements
target_compile_features(physim PUBLIC cxx_std_20)

# Add subdirectories
add_subdirectory(src)

if(PHYSIM_BUILD_TESTS)
    enable_testing()
    add_subdirectory(tests)
endif()

if(PHYSIM_BUILD_BENCHMARKS)
    add_subdirectory(benchmarks)
endif()

if(PHYSIM_BUILD_EXAMPLES)
    add_subdirectory(examples)
endif()

if(PHYSIM_BUILD_PYTHON)
    add_subdirectory(python)
endif()

# Installation rules
include(GNUInstallDirs)

install(TARGETS physim
    EXPORT PhysicsSimEngineTargets
    LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(DIRECTORY include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(EXPORT PhysicsSimEngineTargets
    FILE PhysicsSimEngineTargets.cmake
    NAMESPACE PhysicsSimEngine::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PhysicsSimEngine
)

# Generate version file
include(CMakePackageConfigHelpers)
write_basic_package_version_file(
    "${CMAKE_CURRENT_BINARY_DIR}/PhysicsSimEngineConfigVersion.cmake"
    VERSION ${PROJECT_VERSION}
    COMPATIBILITY AnyNewerVersion
)

# Generate config file
configure_package_config_file(
    "${CMAKE_CURRENT_SOURCE_DIR}/cmake/PhysicsSimEngineConfig.cmake.in"
    "${CMAKE_CURRENT_BINARY_DIR}/PhysicsSimEngineConfig.cmake"
    INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PhysicsSimEngine
)

install(FILES
    "${CMAKE_CURRENT_BINARY_DIR}/PhysicsSimEngineConfig.cmake"
    "${CMAKE_CURRENT_BINARY_DIR}/PhysicsSimEngineConfigVersion.cmake"
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/PhysicsSimEngine
)

# Print configuration summary
message(STATUS "")
message(STATUS "Physics Simulation Engine v${PROJECT_VERSION}")
message(STATUS "================================")
message(STATUS "Build type:              ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ compiler:            ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "C++ standard:            C++${CMAKE_CXX_STANDARD}")
message(STATUS "Build tests:             ${PHYSIM_BUILD_TESTS}")
message(STATUS "Build benchmarks:        ${PHYSIM_BUILD_BENCHMARKS}")
message(STATUS "Build examples:          ${PHYSIM_BUILD_EXAMPLES}")
message(STATUS "Build Python bindings:   ${PHYSIM_BUILD_PYTHON}")
message(STATUS "Enable CUDA:             ${PHYSIM_ENABLE_CUDA}")
message(STATUS "Enable visualization:    ${PHYSIM_ENABLE_VISUALIZATION}")
message(STATUS "Enable AddressSanitizer: ${PHYSIM_USE_ASAN}")
message(STATUS "================================")
message(STATUS "")
