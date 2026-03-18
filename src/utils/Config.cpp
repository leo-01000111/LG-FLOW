#include "utils/Config.hpp"

#include <fstream>
#include <sstream>
#include <stdexcept>

void Config::load(const std::string& filepath)
{
    std::ifstream file(filepath);
    if (!file.is_open())
        throw std::runtime_error("Config: cannot open file '" + filepath + "'");

    std::string line;
    int lineNumber = 0;
    while (std::getline(file, line))
    {
        ++lineNumber;

        // Strip leading/trailing whitespace
        auto start = line.find_first_not_of(" \t\r\n");
        if (start == std::string::npos)
            continue;  // blank line
        line = line.substr(start);

        if (line[0] == '#')
            continue;  // comment

        auto sep = line.find('=');
        if (sep == std::string::npos)
            throw std::runtime_error(
                "Config: malformed line " + std::to_string(lineNumber) +
                " in '" + filepath + "': no '=' found");

        std::string key   = line.substr(0, sep);
        std::string value = line.substr(sep + 1);

        // Trim key and value
        auto trim = [](std::string& s) {
            auto b = s.find_first_not_of(" \t");
            auto e = s.find_last_not_of(" \t\r\n");
            s = (b == std::string::npos) ? "" : s.substr(b, e - b + 1);
        };
        trim(key);
        trim(value);

        m_entries[key] = value;
    }
}

bool Config::has(const std::string& key) const
{
    return m_entries.count(key) > 0;
}
