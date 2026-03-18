#include "utils/Logger.hpp"

#include <gtest/gtest.h>
#include <sstream>
#include <string>

// ── Test fixture: captures stdout and restores Logger state ───────────────────
//
// Logger is a singleton so level changes persist across tests. The fixture
// saves the current level (via setLevel calls) and resets it in TearDown so
// each test starts from a known state.
//
// Output is captured by swapping std::cout's stream buffer with a local
// ostringstream. The original buffer is always restored in TearDown even if
// the test throws, preventing a corrupted std::cout for the rest of the suite.

class LoggerTest : public ::testing::Test
{
protected:
    void SetUp() override
    {
        // Redirect std::cout to m_captured
        m_prevBuf = std::cout.rdbuf(m_captured.rdbuf());
        // Ensure a deterministic starting level for every test
        Logger::get().setLevel(LogLevel::INFO);
    }

    void TearDown() override
    {
        // Always restore cout before any EXPECT macros write to it
        std::cout.rdbuf(m_prevBuf);
        // Reset Logger level to default for subsequent tests
        Logger::get().setLevel(LogLevel::INFO);
    }

    /// Returns everything written to std::cout during this test.
    std::string captured() const { return m_captured.str(); }

    /// Clears the capture buffer (useful for multi-step tests).
    void clearCapture() { m_captured.str(""); m_captured.clear(); }

private:
    std::ostringstream m_captured;
    std::streambuf*    m_prevBuf{nullptr};
};

// ── Level filtering ───────────────────────────────────────────────────────────

TEST_F(LoggerTest, Info_PassesThroughAtDefaultLevel)
{
    Logger::get().info("marker-info");
    EXPECT_NE(captured().find("marker-info"), std::string::npos);
}

TEST_F(LoggerTest, Warn_PassesThroughAtDefaultLevel)
{
    Logger::get().warn("marker-warn");
    EXPECT_NE(captured().find("marker-warn"), std::string::npos);
}

TEST_F(LoggerTest, Error_PassesThroughAtDefaultLevel)
{
    Logger::get().error("marker-error");
    EXPECT_NE(captured().find("marker-error"), std::string::npos);
}

TEST_F(LoggerTest, Debug_FilteredAtDefaultInfoLevel)
{
    Logger::get().debug("marker-debug");
    EXPECT_EQ(captured().find("marker-debug"), std::string::npos);
}

TEST_F(LoggerTest, Debug_PassesThroughWhenLevelSetToDebug)
{
    Logger::get().setLevel(LogLevel::DEBUG);
    Logger::get().debug("marker-debug-enabled");
    EXPECT_NE(captured().find("marker-debug-enabled"), std::string::npos);
}

TEST_F(LoggerTest, Info_FilteredWhenLevelSetToError)
{
    Logger::get().setLevel(LogLevel::ERROR);
    Logger::get().info("marker-info-filtered");
    EXPECT_EQ(captured().find("marker-info-filtered"), std::string::npos);
}

TEST_F(LoggerTest, Warn_FilteredWhenLevelSetToError)
{
    Logger::get().setLevel(LogLevel::ERROR);
    Logger::get().warn("marker-warn-filtered");
    EXPECT_EQ(captured().find("marker-warn-filtered"), std::string::npos);
}

TEST_F(LoggerTest, Error_PassesThroughWhenLevelSetToError)
{
    Logger::get().setLevel(LogLevel::ERROR);
    Logger::get().error("marker-error-passes");
    EXPECT_NE(captured().find("marker-error-passes"), std::string::npos);
}

// ── Output format ─────────────────────────────────────────────────────────────

TEST_F(LoggerTest, OutputFormat_ContainsLevelTag)
{
    Logger::get().info("format-check");
    // Stable format: "[-----] [INFO ] message"
    EXPECT_NE(captured().find("[INFO ]"), std::string::npos);
    EXPECT_NE(captured().find("format-check"), std::string::npos);
}

TEST_F(LoggerTest, OutputFormat_WarnTag)
{
    Logger::get().warn("warn-format");
    EXPECT_NE(captured().find("[WARN ]"), std::string::npos);
}

TEST_F(LoggerTest, OutputFormat_ErrorTag)
{
    Logger::get().error("error-format");
    EXPECT_NE(captured().find("[ERROR]"), std::string::npos);
}

// ── Singleton identity ────────────────────────────────────────────────────────

TEST_F(LoggerTest, Singleton_TwoGetCalls_ReturnSameInstance)
{
    EXPECT_EQ(&Logger::get(), &Logger::get());
}
