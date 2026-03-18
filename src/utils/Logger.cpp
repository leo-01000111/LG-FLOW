#include "utils/Logger.hpp"

#include <chrono>
#include <ctime>
#include <iostream>
#include <string_view>

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
    if (level < m_minLevel)
        return;

    // STUB: timestamp and formatted output will be implemented in Milestone 2
    std::string_view tag;
    switch (level)
    {
        case LogLevel::DEBUG: tag = "DEBUG"; break;
        case LogLevel::INFO:  tag = "INFO "; break;
        case LogLevel::WARN:  tag = "WARN "; break;
        case LogLevel::ERROR: tag = "ERROR"; break;
    }

    std::cout << "[-----] [" << tag << "] " << message << '\n';
}

void Logger::debug(std::string_view message) { log(LogLevel::DEBUG, message); }
void Logger::info (std::string_view message) { log(LogLevel::INFO,  message); }
void Logger::warn (std::string_view message) { log(LogLevel::WARN,  message); }
void Logger::error(std::string_view message) { log(LogLevel::ERROR, message); }
