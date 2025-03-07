//
// Created by Leon Kang on 2025/3/6.
//

#ifndef SOLVER_H
#define SOLVER_H

#include <vector>
#include <span>
#include "ode.h"

namespace pece {
    template <int dim>
    class Solver {
    public:
        struct StartingValues {
            Eigen::Vector<float, dim> y0;
            float t0, h0; // initial time and step size
        };
        struct State {
            Eigen::Vector<float, dim> y_prev, y_curr, f_prev, f_curr;
            float t_prev, t_curr;
        };
        // No adaptivity, does not require second and third derivatives
        [[nodiscard]] auto Solve(const ODESystem<dim>& ode, const std::span<const float> ts, const StartingValues& sv) {
            std::vector<Eigen::Vector<float, dim>> ys;
            ys.reserve(ts.size());
            auto [y_curr, t_curr, h] = sv;
            auto f_curr = ode.F1(t_curr, y_curr);
            // Use midpoint method to obtain starting values not provided
            auto y_midpoint = y_curr - 0.5f * h * ode.F1(t_curr, y_curr);
            auto y_prev = y_curr - h * ode.F1(t_curr - 0.5f * h, y_midpoint);
            auto f_prev = ode.F1(t_curr - h, y_prev);
            auto iter = ts.begin();
            while (iter != ts.end()) {
                auto t_next = t_curr + h;
                // Prediction using 2nd order Adams-Bashforth method
                auto y_prediction = y_curr + 0.5f * h * (3.0f * f_curr - f_prev);
                auto f_prediction = ode.F1(t_next, y_prediction);
                // Correction using 2nd order Adams-Moulton method
                auto y_next = y_curr + 0.5f * h * (f_curr + f_prediction);
                auto f_next = ode.F1(t_next, y_next);
                while (iter != ts.end() && *iter <= t_next) {
                    ys.emplace_back(y_curr + 0.5f * (*iter++ - t_curr) * (3.0f * f_next - f_curr));
                }
                f_prev = f_curr;
                y_curr = y_next;
                f_curr = f_next;
                t_curr = t_next;
            }
            return _flatten(ys);
        }
        // With adaptive step size based on error estimation and error tolerance
        [[nodiscard]] auto Solve(const ODESystem<dim>& ode, const std::span<const float> ts, const StartingValues& sv, float tol) {
            std::vector<Eigen::Vector<float, dim>> ys;
            ys.reserve(ts.size());
            std::vector<float> hs;
            hs.reserve(ts.size());
            auto [y_curr, t_curr, h_curr] = sv;
            auto f_curr = ode.F1(t_curr, y_curr);
            // Use forward Euler method to obtain starting values not provided, since the problem is non-stiff
            auto h_prev = std::min(h_curr, std::sqrtf(1.8f * tol / (ode.F2(t_curr, y_curr).norm() + 1e-6f)));
            auto y_prev = y_curr - h_prev * f_curr;
            auto f_prev = ode.F1(t_curr - h_prev, y_prev);
            auto iter = ts.begin();
            while (iter != ts.end()) {
                while (true) {
                    auto t_next = t_curr + h_curr;
                    // Prediction using 2nd order Adams-Bashforth method
                    auto y_prediction = y_curr + h_curr * f_curr + 0.5f * h_curr * h_curr / h_prev * (f_curr - f_prev);
                    auto f_prediction = ode.F1(t_next, y_prediction);
                    // Correction using 2nd order Adams-Moulton method
                    auto y_next = y_curr + 0.5f * h_curr * (f_curr + f_prediction);
                    auto f_next = ode.F1(t_next, y_next);
                    // Error estimation using 2nd order Adams-Moulton method
                    auto error = (h_curr * h_curr * h_curr / 12.0f * ode.F3(t_next, y_next)).norm();
                    if (error < tol) {
                        while (iter != ts.end() && *iter <= t_next) {
                            auto h_star = *iter - t_curr;
                            ys.emplace_back(y_curr + h_star * f_next + 0.5f * h_star * h_star / h_curr * (f_next - f_curr));
                            hs.emplace_back(h_curr);
                            ++iter;
                        }
                        h_prev = h_curr;
                        h_curr *= error == 0.0f ? 1.0f : std::powf(0.9f * tol / error, 1.0f / 3.0f);
                        f_prev = f_curr;
                        y_curr = y_next;
                        f_curr = f_next;
                        t_curr = t_next;
                        break;
                    }
                    h_curr *= error == 0.0f ? 1.0f : std::powf(0.9f * tol / error, 1.0f / 3.0f);
                }
            }
            return std::make_pair(_flatten(ys), hs);
        }
    private:
        [[nodiscard]] static auto _flatten(std::vector<Eigen::Vector<float, dim>> ys) {
            std::array<std::vector<float>, dim> flattened_ys;
            for (auto i = 0; i < dim; ++i) {
                flattened_ys[i].resize(ys.size());
                std::ranges::transform(ys, flattened_ys[i].begin(), [i](const Eigen::Vector<float, dim>& y) {
                    return y(i);
                });
            }
            return flattened_ys;
        }
    };
} // namespace pece

#endif //SOLVER_H
