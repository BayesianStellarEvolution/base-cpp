// Copyright (c) Embry-Riddle Aeronautical University. All rights reserved.

#include <gtest/gtest.h>

#include <cmath>
#include <string>

#include "yaml-cpp/yaml.h"

class ScatterConfigTest : public ::testing::Test {
protected:
    std::string confDir;

    void SetUp() override {
        confDir = std::string(TEST_CONF_DIR) + "/";
    }
};

TEST_F(ScatterConfigTest, LoadUncrowdedConfig) {
    std::string configPath = confDir + "scatter_uncrowded.yaml";
    YAML::Node config;

    ASSERT_NO_THROW(config = YAML::LoadFile(configPath))
        << "Failed to load " << configPath;

    EXPECT_TRUE(config["brightLimit"]) << "Missing brightLimit";
    EXPECT_TRUE(config["faintLimit"]) << "Missing faintLimit";
    EXPECT_TRUE(config["relevantFilt"]) << "Missing relevantFilt";
    EXPECT_TRUE(config["limitS2N"]) << "Missing limitS2N";
    EXPECT_TRUE(config["minSigma"]) << "Missing minSigma";
    EXPECT_TRUE(config["exposures"]) << "Missing exposures";
    EXPECT_TRUE(config["s2nCoeffs"]) << "Missing s2nCoeffs";

    EXPECT_DOUBLE_EQ(config["minSigma"].as<double>(), 0.01);

    auto exposures = config["exposures"];
    EXPECT_TRUE(exposures["U"]) << "Missing U in exposures";
    EXPECT_TRUE(exposures["V"]) << "Missing V in exposures";
    EXPECT_TRUE(exposures["G"]) << "Missing G (Gaia) in exposures";

    auto s2n = config["s2nCoeffs"];
    EXPECT_TRUE(s2n["U"]) << "Missing U in s2nCoeffs";
    EXPECT_TRUE(s2n["G"]) << "Missing G (Gaia) in s2nCoeffs";

    auto uCoeffs = s2n["U"].as<std::vector<double>>();
    ASSERT_EQ(uCoeffs.size(), 2u);
    EXPECT_NEAR(uCoeffs[0], 9.33989, 1e-4);
    EXPECT_NEAR(uCoeffs[1], 0.3375778, 1e-6);
}

TEST_F(ScatterConfigTest, LoadCrowdedConfig) {
    std::string configPath = confDir + "scatter_crowded.yaml";
    YAML::Node config;

    ASSERT_NO_THROW(config = YAML::LoadFile(configPath))
        << "Failed to load " << configPath;

    EXPECT_TRUE(config["brightLimit"]) << "Missing brightLimit";
    EXPECT_TRUE(config["faintLimit"]) << "Missing faintLimit";
    EXPECT_TRUE(config["relevantFilt"]) << "Missing relevantFilt";
    EXPECT_TRUE(config["limitS2N"]) << "Missing limitS2N";
    EXPECT_TRUE(config["minSigma"]) << "Missing minSigma";
    EXPECT_TRUE(config["exposures"]) << "Missing exposures";
    EXPECT_TRUE(config["s2nCoeffs"]) << "Missing s2nCoeffs";

    EXPECT_DOUBLE_EQ(config["minSigma"].as<double>(), 0.04);

    auto s2n = config["s2nCoeffs"];
    auto uvCoeffs = s2n["UVf275w"].as<std::vector<double>>();
    ASSERT_EQ(uvCoeffs.size(), 2u);
    EXPECT_NEAR(uvCoeffs[0], 6.7302, 1e-4);
    EXPECT_NEAR(uvCoeffs[1], 0.2090, 1e-4);
}

TEST_F(ScatterConfigTest, S2nCoeffsHaveTwoElements) {
    for (const auto& configName : {"scatter_uncrowded.yaml", "scatter_crowded.yaml"}) {
        std::string configPath = confDir + configName;
        YAML::Node config = YAML::LoadFile(configPath);

        auto s2n = config["s2nCoeffs"];
        for (auto it = s2n.begin(); it != s2n.end(); ++it) {
            std::string filter = it->first.as<std::string>();
            auto coeffs = it->second.as<std::vector<double>>();
            EXPECT_EQ(coeffs.size(), 2u)
                << "Filter " << filter << " in " << configName
                << " should have exactly 2 S/N coefficients";
        }
    }
}

TEST_F(ScatterConfigTest, ExposuresArePositive) {
    for (const auto& configName : {"scatter_uncrowded.yaml", "scatter_crowded.yaml"}) {
        std::string configPath = confDir + configName;
        YAML::Node config = YAML::LoadFile(configPath);

        auto exposures = config["exposures"];
        for (auto it = exposures.begin(); it != exposures.end(); ++it) {
            std::string filter = it->first.as<std::string>();
            double expTime = it->second.as<double>();
            EXPECT_GT(expTime, 0.0)
                << "Filter " << filter << " in " << configName
                << " should have positive exposure time";
        }
    }
}

TEST_F(ScatterConfigTest, ConfigsHaveMatchingS2nForExposures) {
    for (const auto& configName : {"scatter_uncrowded.yaml", "scatter_crowded.yaml"}) {
        std::string configPath = confDir + configName;
        YAML::Node config = YAML::LoadFile(configPath);

        auto exposures = config["exposures"];
        auto s2n = config["s2nCoeffs"];

        for (auto it = exposures.begin(); it != exposures.end(); ++it) {
            std::string filter = it->first.as<std::string>();
            EXPECT_TRUE(s2n[filter])
                << "Filter " << filter << " has exposure but no S/N coefficients in "
                << configName;
        }
    }
}
