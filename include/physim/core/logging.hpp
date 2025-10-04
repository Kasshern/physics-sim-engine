/**
 * @file logging.hpp
 * @brief Logging infrastructure using spdlog
 *
 * Provides centralized logging for the physics simulation engine.
 * Uses spdlog for high-performance, asynchronous logging.
 */

#pragma once

#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <memory>
#include <string>

namespace physim {
namespace logging {

/**
 * @brief Initialize logging system
 * @param level Log level (trace, debug, info, warn, error, critical)
 * @param log_to_file Enable file logging
 * @param log_file Path to log file
 *
 * Sets up console and optionally file logging with colored output.
 */
void init(spdlog::level::level_enum level = spdlog::level::info,
          bool log_to_file = false,
          const std::string& log_file = "physim.log");

/**
 * @brief Get the default logger
 * @return Shared pointer to logger
 */
std::shared_ptr<spdlog::logger> get_logger();

/**
 * @brief Set logging level
 * @param level New log level
 */
void set_level(spdlog::level::level_enum level);

/**
 * @brief Flush all loggers
 *
 * Useful before program termination or checkpointing.
 */
void flush();

/**
 * @brief Shutdown logging system
 *
 * Call before program exit to ensure all logs are written.
 */
void shutdown();

} // namespace logging
} // namespace physim

// Convenience macros for logging
#define PHYSIM_LOG_TRACE(...)    physim::logging::get_logger()->trace(__VA_ARGS__)
#define PHYSIM_LOG_DEBUG(...)    physim::logging::get_logger()->debug(__VA_ARGS__)
#define PHYSIM_LOG_INFO(...)     physim::logging::get_logger()->info(__VA_ARGS__)
#define PHYSIM_LOG_WARN(...)     physim::logging::get_logger()->warn(__VA_ARGS__)
#define PHYSIM_LOG_ERROR(...)    physim::logging::get_logger()->error(__VA_ARGS__)
#define PHYSIM_LOG_CRITICAL(...) physim::logging::get_logger()->critical(__VA_ARGS__)
