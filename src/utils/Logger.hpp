#pragma once

#include <string>
#include <string_view>

/**
 * @brief Severity levels for log messages.
 */
enum class LogLevel
{
    DEBUG,
    INFO,
    WARN,
    ERROR
};

/**
 * @brief Simple timestamped logger with severity levels.
 *
 * All output goes to stdout in the format:
 *   [HH:MM:SS] [LEVEL] message
 *
 * Access the singleton via Logger::get().
 */
class Logger
{
public:
    /**
     * @brief Returns the global Logger singleton.
     * @return Reference to the single Logger instance.
     */
    static Logger& get();

    /**
     * @brief Logs a message at the given severity level.
     * @param level   Severity of the message.
     * @param message Human-readable log text.
     */
    void log(LogLevel level, std::string_view message);

    /** @brief Convenience wrapper: logs at DEBUG level. */
    void debug(std::string_view message);

    /** @brief Convenience wrapper: logs at INFO level. */
    void info(std::string_view message);

    /** @brief Convenience wrapper: logs at WARN level. */
    void warn(std::string_view message);

    /** @brief Convenience wrapper: logs at ERROR level. */
    void error(std::string_view message);

    /**
     * @brief Sets the minimum level that will be printed.
     *
     * Messages below this level are silently discarded.
     * Default: LogLevel::INFO.
     *
     * @param level Minimum severity to display.
     */
    void setLevel(LogLevel level);

    // Non-copyable singleton
    Logger(const Logger&)            = delete;
    Logger& operator=(const Logger&) = delete;

private:
    Logger();

    LogLevel m_minLevel;
};
