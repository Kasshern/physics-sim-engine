/**
 * @file logging.cpp
 * @brief Implementation of logging infrastructure
 */

#include "physim/core/logging.hpp"
#include <spdlog/spdlog.h>
#include <spdlog/sinks/stdout_color_sinks.h>
#include <spdlog/sinks/basic_file_sink.h>
#include <spdlog/sinks/rotating_file_sink.h>
#include <iostream>
#include <vector>

namespace physim {
namespace logging {

namespace {
    std::shared_ptr<spdlog::logger> g_logger = nullptr;
}

void init(spdlog::level::level_enum level, bool log_to_file, const std::string& log_file) {
    if (g_logger) {
        // Already initialized
        return;
    }

    try {
        std::vector<spdlog::sink_ptr> sinks;

        // Console sink with colors
        auto console_sink = std::make_shared<spdlog::sinks::stdout_color_sink_mt>();
        console_sink->set_level(level);
        console_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%^%l%$] [%n] %v");
        sinks.push_back(console_sink);

        // File sink if requested
        if (log_to_file) {
            // Rotating file sink: 10MB max size, 3 rotated files
            auto file_sink = std::make_shared<spdlog::sinks::rotating_file_sink_mt>(
                log_file, 1024 * 1024 * 10, 3);
            file_sink->set_level(spdlog::level::trace); // Log everything to file
            file_sink->set_pattern("[%Y-%m-%d %H:%M:%S.%e] [%l] [%n] %v");
            sinks.push_back(file_sink);
        }

        // Create logger with both sinks
        g_logger = std::make_shared<spdlog::logger>("physim", sinks.begin(), sinks.end());
        g_logger->set_level(level);
        g_logger->flush_on(spdlog::level::warn); // Auto-flush on warnings and errors

        // Register as default logger
        spdlog::set_default_logger(g_logger);

        PHYSIM_LOG_INFO("Physics simulation engine logging initialized");
        PHYSIM_LOG_INFO("Log level: {}", spdlog::level::to_string_view(level));
        if (log_to_file) {
            PHYSIM_LOG_INFO("Logging to file: {}", log_file);
        }
    }
    catch (const spdlog::spdlog_ex& ex) {
        std::cerr << "Log initialization failed: " << ex.what() << std::endl;
    }
}

std::shared_ptr<spdlog::logger> get_logger() {
    if (!g_logger) {
        // Initialize with defaults if not already done
        init();
    }
    return g_logger;
}

void set_level(spdlog::level::level_enum level) {
    if (g_logger) {
        g_logger->set_level(level);
        PHYSIM_LOG_INFO("Log level changed to: {}", spdlog::level::to_string_view(level));
    }
}

void flush() {
    if (g_logger) {
        g_logger->flush();
    }
}

void shutdown() {
    if (g_logger) {
        PHYSIM_LOG_INFO("Shutting down logging system");
        g_logger->flush();
        spdlog::shutdown();
        g_logger.reset();
    }
}

} // namespace logging
} // namespace physim
