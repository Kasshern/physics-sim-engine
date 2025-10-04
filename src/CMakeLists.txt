# Source files for the physics simulation library

target_sources(physim PRIVATE
    # Core module
    core/types.cpp
    core/time.cpp
    core/frame.cpp
    core/logging.cpp

    # Integrators module (Phase 2)
    integrators/integrator.cpp
    integrators/rk4.cpp
    integrators/rk45.cpp

    # Forces module (Phase 2)
    forces/force.cpp
    forces/gravity.cpp
    forces/j2_gravity.cpp

    # Future modules
    # integrators/dopri.cpp
    # nbody/particle.cpp
    # nbody/system.cpp
    # nbody/propagator.cpp
    # attitude/quaternion.cpp
    # attitude/dynamics.cpp
    # io/checkpoint.cpp
    # io/telemetry.cpp
)

# Add CUDA sources if enabled
if(PHYSIM_ENABLE_CUDA)
    target_sources(physim PRIVATE
        # gpu/kernels.cu
        # gpu/device_vector.cpp
        # gpu/reduction.cu
    )

    set_source_files_properties(
        # gpu/kernels.cu
        # gpu/reduction.cu
        PROPERTIES LANGUAGE CUDA
    )

    target_link_libraries(physim PUBLIC CUDA::cudart)
endif()

# Add visualization sources if enabled
if(PHYSIM_ENABLE_VISUALIZATION)
    target_sources(physim PRIVATE
        # viz/renderer.cpp
        # viz/camera.cpp
        # viz/scene.cpp
        # viz/gui.cpp
    )

    target_link_libraries(physim PUBLIC
        glfw
        imgui::imgui
    )
endif()
