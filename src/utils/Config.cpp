#include "utils/Config.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

// ── Duplicate-key policy ──────────────────────────────────────────────────────
// If the same key appears more than once in a config file, the LAST value wins.
// This is intentional: it allows a base config to be overridden by appending
// lines, which is convenient for derived test cases.
// Rationale for last-wins over error: simpler composition, no silent data loss
// (the final value is always the most recent one seen).

void Config::load(const std::string& filepath)
{
    // Clear any previously loaded state so reload does not retain stale keys.
    m_entries.clear();

    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Config: cannot open file '" + filepath + "'");

    std::string line;
    int lineNumber = 0;
    while (std::getline(file, line))
    {
        ++lineNumber;

        // Strip leading whitespace before anything else
        auto start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos)
            continue;  // blank line
        line = line.substr(start);

        if (line[0] == '#')
            continue;  // comment line

        auto sep = line.find('=');
        if (sep == std::string::npos)
            throw std::runtime_error(
                "Config: malformed line " + std::to_string(lineNumber) +
                " in '" + filepath + "': no '=' separator found");

        std::string key   = line.substr(0, sep);
        std::string value = line.substr(sep + 1);

        // Trim both key and value of surrounding whitespace
        auto trim = [](std::string& s) {
            auto b = s.find_first_not_of(" \t");
            auto e = s.find_last_not_of(" \t\r\n");
            s = (b == std::string::npos) ? "" : s.substr(b, e - b + 1);
        };
        trim(key);
        trim(value);

        if (key.empty())
            throw std::runtime_error(
                "Config: malformed line " + std::to_string(lineNumber) +
                " in '" + filepath + "': key is empty");

        // last-wins: overwrite any previous entry for this key
        m_entries[key] = value;
    }
}

bool Config::has(const std::string& key) const
{
    return m_entries.count(key) > 0;
}
