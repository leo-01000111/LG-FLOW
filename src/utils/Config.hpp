#pragma once

#include <optional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>

/**
 * @brief Parses and stores key=value configuration files.
 *
 * Supports a simple line-oriented format:
 *   - Lines starting with '#' are comments.
 *   - Blank lines are ignored.
 *   - All other lines must have the form: key = value
 *
 * Example:
 * @code
 *   Config cfg;
 *   cfg.load("case.cfg");
 *   double dt = cfg.get<double>("solver.dt");
 *   int nx    = cfg.get<int>("mesh.Nx", 64);   // default 64
 * @endcode
 */
class Config
{
public:
    Config() = default;

    /**
     * @brief Loads key=value pairs from a file on disk.
     *
     * @param filepath Path to the configuration file.
     * @throws std::runtime_error if the file cannot be opened or a line is malformed.
     */
    void load(const std::string& filepath);

    /**
     * @brief Returns the value for a key, converting it to type T.
     *
     * @tparam T Target type (must be parseable from a string via std::istringstream).
     * @param key The configuration key to look up.
     * @return The parsed value.
     * @throws std::runtime_error if the key does not exist or conversion fails.
     */
    template <typename T>
    [[nodiscard]] T get(const std::string& key) const;

    /**
     * @brief Returns the value for a key, or a default if the key is absent.
     *
     * @tparam T Target type.
     * @param key          The configuration key to look up.
     * @param defaultValue Value returned when the key is missing.
     * @return The parsed value, or defaultValue.
     */
    template <typename T>
    [[nodiscard]] T get(const std::string& key, const T& defaultValue) const;

    /**
     * @brief Checks whether a key exists in the loaded configuration.
     * @param key The key to query.
     * @return true if the key is present.
     */
    [[nodiscard]] bool has(const std::string& key) const;

private:
    std::unordered_map<std::string, std::string> m_entries;
};

// ── Template implementations ─────────────────────────────────────────────────

template <typename T>
T Config::get(const std::string& key) const
{
    auto it = m_entries.find(key);
    if (it == m_entries.end())
        throw std::runtime_error("Config: missing required key '" + key + "'");

    std::istringstream ss(it->second);
    T value{};
    if (!(ss >> value))
        throw std::runtime_error("Config: cannot parse value for key '" + key + "'");

    return value;
}

template <typename T>
T Config::get(const std::string& key, const T& defaultValue) const
{
    if (!has(key))
        return defaultValue;

    return get<T>(key);
}
