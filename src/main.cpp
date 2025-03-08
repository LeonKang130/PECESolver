//
// Created by Leon Kang on 2025/3/6.
//

#include <fstream>
#include <iostream>
#include <gtest/gtest.h>
#include <nlohmann/json.hpp>
#include "ode.h"
#include "solver.h"
#include "plot.h"

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
    pece::Figure figure{"Part 1 - Fixed Step Size"};
    figure.AddPlot("y1", ts, ys[0]);
    figure.AddPlot("y2", ts, ys[1]);
    std::ofstream file("data/part-1.json");
    file << figure.ToString();
    file.close();
}

TEST(Part2, Problem1LowAccuracy) {
    const auto ode = pece::PredatorPrey();
    pece::Solver<2> solver;
    constexpr auto num_points = 256;
    const auto ts = linspace(0.0f, 100.0f, num_points);
    const auto [ys, hs] = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{10.0f, 10.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-3f);
    pece::Figure figure1{"Part 2 - Predator Prey (Low Accuracy, Solution)"};
    figure1.AddPlot("y1", ts, ys[0]);
    figure1.AddPlot("y2", ts, ys[1]);
    std::ofstream file1("data/part-2-1-1.json");
    file1 << figure1.ToString();
    file1.close();
    pece::Figure figure2{"Part 2 - Predator Prey (Low Accuracy, Step Size)"};
    figure2.AddPlot("h", ts, hs);
    std::ofstream file2("data/part-2-1-2.json");
    file2 << figure2.ToString();
    file2.close();
    pece::Figure figure3{"Part 2 - Predator Prey (Low Accuracy, Relative)"};
    figure3.AddPlot("y1", ys[1], ys[0]);
    std::ofstream file3("data/part-2-1-3.json");
    file3 << figure3.ToString();
    file3.close();
}

TEST(Part2, Problem1HighAccuracy) {
    const auto ode = pece::PredatorPrey();
    pece::Solver<2> solver;
    constexpr auto num_points = 256;
    const auto ts = linspace(0.0f, 100.0f, num_points);
    const auto [ys, hs] = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{10.0f, 10.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-6f);
    pece::Figure figure1{"Part 2 - Predator Prey (High Accuracy, Solution)"};
    figure1.AddPlot("y1", ts, ys[0]);
    figure1.AddPlot("y2", ts, ys[1]);
    std::ofstream file1("data/part-2-1-4.json");
    file1 << figure1.ToString();
    file1.close();
    pece::Figure figure2{"Part 2 - Predator Prey (High Accuracy, Step Size)"};
    figure2.AddPlot("h", ts, hs);
    std::ofstream file2("data/part-2-1-5.json");
    file2 << figure2.ToString();
    file2.close();
    pece::Figure figure3{"Part 2 - Predator Prey (High Accuracy, Relative)"};
    figure3.AddPlot("y1", ys[1], ys[0]);
    std::ofstream file3("data/part-2-1-6.json");
    file3 << figure3.ToString();
    file3.close();
}

TEST(Part2, Problem2LowAccuracy) {
    const auto ode = pece::VanDerPol();
    pece::Solver<2> solver;
    constexpr auto num_points = 256;
    const auto ts = linspace(0.0f, 11.0f, num_points);
    const auto [ys, hs] = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{2.0f, 0.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-3f);
    pece::Figure figure1{"Part 2 - Van der Pol's (Low Accuracy, Solution)"};
    figure1.AddPlot("y1", ts, ys[0]);
    figure1.AddPlot("y2", ts, ys[1]);
    std::ofstream file1("data/part-2-2-1.json");
    file1 << figure1.ToString();
    file1.close();
    pece::Figure figure2{"Part 2 - Van der Pol's (Low Accuracy, Step Size)"};
    figure2.AddPlot("h", ts, hs);
    std::ofstream file2("data/part-2-2-2.json");
    file2 << figure2.ToString();
    file2.close();
}

TEST(Part2, Problem2HighAccuracy) {
    const auto ode = pece::VanDerPol();
    pece::Solver<2> solver;
    constexpr auto num_points = 256;
    const auto ts = linspace(0.0f, 11.0f, num_points);
    const auto [ys, hs] = solver.Solve(ode, ts, pece::Solver<2>::StartingValues{
        .y0 = Eigen::Vector2f{2.0f, 0.0f},
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-6f);
    pece::Figure figure1{"Part 2 - Van der Pol's (High Accuracy, Solution)"};
    figure1.AddPlot("y1", ts, ys[0]);
    figure1.AddPlot("y2", ts, ys[1]);
    std::ofstream file1("data/part-2-2-3.json");
    file1 << figure1.ToString();
    file1.close();
    pece::Figure figure2{"Part 2 - Van der Pol's (High Accuracy, Step Size)"};
    figure2.AddPlot("h", ts, hs);
    std::ofstream file2("data/part-2-2-4.json");
    file2 << figure2.ToString();
    file2.close();
}

TEST(Part2, Problem3LowAccuracy) {
    constexpr auto num_points = 256;
    const auto ode = pece::MethodOfLines<num_points>();
    pece::Solver<num_points> solver;
    std::vector ts{0.0f, 0.25f, 0.5f, 0.6f, 0.8f, 1.0f};
    const auto xs = linspace(0.0f, 1.0f, num_points);
    Eigen::Vector<float, num_points> y0;
    for (auto i = 0; i < num_points; ++i) {
        y0(i) = std::expf(-10.0f * xs[i]);
    }
    const auto [ys, hs] = solver.Solve(ode, ts, pece::Solver<num_points>::StartingValues{
        .y0 = y0,
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-3f);
    pece::Figure figure1{"Part 2 - Method Of Lines (Low Accuracy, Solution)"};
    std::vector<std::array<float, num_points>> us;
    us.resize(ts.size());
    for (auto i = 0; i < ts.size(); ++i) {
        for (auto j = 0; j < num_points; ++j) {
            us[i][j] = ys[j][i];
        }
        figure1.AddPlot(std::format("t = {}", ts[i]), xs, us[i]);
    }
    std::ofstream file1("data/part-2-3-1.json");
    file1 << figure1.ToString();
    file1.close();
    pece::Figure figure2{"Part 2 - Method Of Lines (Low Accuracy, Step Size)"};
    figure2.AddPlot("h", ts, hs);
    std::ofstream file2("data/part-2-3-2.json");
    file2 << figure2.ToString();
    file2.close();
}

TEST(Part2, Problem3HighAccuracy) {
    constexpr auto num_points = 256;
    const auto ode = pece::MethodOfLines<num_points>();
    pece::Solver<num_points> solver;
    std::vector ts{0.0f, 0.25f, 0.5f, 0.6f, 0.8f, 1.0f};
    const auto xs = linspace(0.0f, 1.0f, num_points);
    Eigen::Vector<float, num_points> y0;
    for (auto i = 0; i < num_points; ++i) {
        y0(i) = std::expf(-10.0f * xs[i]);
    }
    const auto [ys, hs] = solver.Solve(ode, ts, pece::Solver<num_points>::StartingValues{
        .y0 = y0,
        .t0 = 0.0f,
        .h0 = 0.01f
    }, 1e-6f);
    pece::Figure figure1{"Part 2 - Method Of Lines (High Accuracy, Solution)"};
    std::vector<std::array<float, num_points>> us;
    us.resize(ts.size());
    for (auto i = 0; i < ts.size(); ++i) {
        for (auto j = 0; j < num_points; ++j) {
            us[i][j] = ys[j][i];
        }
        figure1.AddPlot(std::format("t = {}", ts[i]), xs, us[i]);
    }
    std::ofstream file1("data/part-2-3-3.json");
    file1 << figure1.ToString();
    file1.close();
    pece::Figure figure2{"Part 2 - Method Of Lines (High Accuracy, Step Size)"};
    figure2.AddPlot("h", ts, hs);
    std::ofstream file2("data/part-2-3-4.json");
    file2 << figure2.ToString();
    file2.close();
}

int main(int argc, char **argv) {
    testing::InitGoogleTest(&argc, argv);
    return RUN_ALL_TESTS();
}