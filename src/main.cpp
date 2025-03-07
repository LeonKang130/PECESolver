//
// Created by Leon Kang on 2025/3/6.
//

#include <fstream>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include "ode.h"
#include "solver.h"

[[nodiscard]] auto linspace(float start, float end, int n) {
    std::vector<float> x;
    x.resize(n);
    std::ranges::generate(x, [start, end, n, i = 0]() mutable {
        return start + static_cast<float>(i++) * (end - start) / static_cast<float>(n - 1);
    });
    return x;
}

TEST(Part1, Problem1) {
    using ODESystem = pece::ODESystem<2>;
    auto f1 = [](const float t, const Eigen::Vector2f& y) noexcept {
        return Eigen::Vector2f{
            -y(0),
            -10.0f * (y(1) - t * t) + 2.0f * t
        };
    };
    const ODESystem ode{f1};
    pece::Solver<2> solver;
    constexpr auto num_points = 64;
    const auto ts = linspace(0.0f, 1.0f, num_points);
    const auto ys = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{1.0f, 2.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    });
    nlohmann::json json;
    json["ts"] = ts;
    json["ys"] = ys;
    json["title"] = "Part 1 - Fixed Step Size";
    std::ofstream file("part-1.json");
    file << json.dump(4);
    file.close();
}

TEST(Part2, Problem1LowPrecision) {
    const auto ode = pece::PredatorPrey();
    pece::Solver<2> solver;
    constexpr auto num_points = 128;
    const auto ts = linspace(0.0f, 100.0f, num_points);
    const auto ys = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{10.0f, 10.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-3f);
    nlohmann::json json;
    json["ts"] = ts;
    json["ys"] = ys;
    json["title"] = "Part 2 - Predator-Prey Problem (Low Precision)";
    std::ofstream file("part-2-1-1.json");
    file << json.dump(4);
    file.close();
}

TEST(Part2, Problem1HighPrecision) {
    const auto ode = pece::PredatorPrey();
    pece::Solver<2> solver;
    constexpr auto num_points = 128;
    const auto ts = linspace(0.0f, 100.0f, num_points);
    const auto ys = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{10.0f, 10.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-6f);
    nlohmann::json json;
    json["ts"] = ts;
    json["ys"] = ys;
    json["title"] = "Part 2 - Predator-Prey Problem (High Precision)";
    std::ofstream file("part-2-1-2.json");
    file << json.dump(4);
    file.close();
}

TEST(Part2, Problem2LowPrecision) {
    const auto ode = pece::VanDerPol();
    pece::Solver<2> solver;
    constexpr auto num_points = 128;
    const auto ts = linspace(0.0f, 11.0f, num_points);
    const auto ys = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{2.0f, 0.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-3f);
    nlohmann::json json;
    json["ts"] = ts;
    json["ys"] = ys;
    json["title"] = "Part 2 - Van der Pol's equation (Low Precision)";
    std::ofstream file("part-2-2-1.json");
    file << json.dump(4);
    file.close();
}

TEST(Part2, Problem2HighPrecision) {
    const auto ode = pece::VanDerPol();
    pece::Solver<2> solver;
    constexpr auto num_points = 128;
    const auto ts = linspace(0.0f, 11.0f, num_points);
    const auto ys = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{2.0f, 0.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-6f);
    nlohmann::json json;
    json["ts"] = ts;
    json["ys"] = ys;
    json["title"] = "Part 2 - Van der Pol's equation (High Precision)";
    std::ofstream file("part-2-2-2.json");
    file << json.dump(4);
    file.close();
}

TEST(Part2, Problem3LowPrecision) {
    constexpr auto num_points = 128;
    using ODESystem = pece::ODESystem<num_points>;
    std::vector<float> ts{0.0f, 0.25f, 0.5f, 0.6f, 0.8f, 1.0f};
    pece::Solver<num_points> solver;

}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}