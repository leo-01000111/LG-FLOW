#include "utils/Logger.hpp"

#include <iostream>

// ── Singleton ────────────────────────────────────────────────────────────────

Logger& Logger::get()
{
    static Logger instance;
    return instance;
}

Logger::Logger()
    : m_minLevel(LogLevel::INFO)
{
}

// ── Public interface ─────────────────────────────────────────────────────────

void Logger::setLevel(LogLevel level)
{
    m_minLevel = level;
}

void Logger::log(LogLevel level, std::string_view message)
{
    // Filter: messages below the configured minimum level are silently dropped.
    // LogLevel is ordered DEBUG(0) < INFO(1) < WARN(2) < ERROR(3); the enum
    // underlying-int comparison makes this work without explicit numeric casts.
    if (level < m_minLevel)
        return;

    // Output format (stable for tests):
    //   [-----] [LEVEL] message\n
    //
    // The "[-----]" timestamp placeholder will be replaced with a real
    // HH:MM:SS timestamp in Milestone 2. The format of the LEVEL tag is
    // fixed-width (5 chars) so log lines are easy to parse/grep.
    std::string_view tag;
    switch (level)
    {
        case LogLevel::DEBUG: tag = "DEBUG"; break;
        case LogLevel::INFO:  tag = "INFO "; break;
        case LogLevel::WARN:  tag = "WARN "; break;
        case LogLevel::ERROR: tag = "ERROR"; break;
        default:              tag = "?????"; break;
    }

    std::cout << "[-----] [" << tag << "] " << message << '\n';
}

void Logger::debug(std::string_view message) { log(LogLevel::DEBUG, message); }
void Logger::info (std::string_view message) { log(LogLevel::INFO,  message); }
void Logger::warn (std::string_view message) { log(LogLevel::WARN,  message); }
void Logger::error(std::string_view message) { log(LogLevel::ERROR, message); }
