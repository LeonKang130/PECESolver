//
// Created by Leon Kang on 2025/3/6.
//

#ifndef ODE_H
#define ODE_H

#include <functional>
#include <stdexcept>
#include <Eigen/Dense>
#include <spdlog/spdlog.h>

namespace pece {
    template <int dim>
    class ODESystem {
        struct DefaultDerivative {
            auto operator()(float, const Eigen::Vector<float, dim>&) -> Eigen::Vector<float, dim> {
                spdlog::critical("Derivative function not provided for ODE system");
                throw std::runtime_error("Not implemented");
            }
        };
        using Derivative =
            std::function<Eigen::Vector<float, dim>(float, const Eigen::Vector<float, dim>&)>;
        Derivative _f1, _f2, _f3;
    public:
        explicit ODESystem(Derivative&& f1, Derivative&& f2 = DefaultDerivative{}, Derivative&& f3 = DefaultDerivative{})
            : _f1{std::move(f1)}, _f2 {std::move(f2)}, _f3 {std::move(f3)} {}
        [[nodiscard]] static auto Dimension() noexcept { return dim; }
        [[nodiscard]] auto F1(float t, const Eigen::Vector<float, dim>& y) const { return _f1(t, y); }
        [[nodiscard]] auto F2(float t, const Eigen::Vector<float, dim>& y) const { return _f2(t, y); }
        [[nodiscard]] auto F3(float t, const Eigen::Vector<float, dim>& y) const { return _f3(t, y); }
    };

    static auto PredatorPrey() noexcept {
        auto f1 = [](const float t, const Eigen::Vector2f& y) noexcept {
            return Eigen::Vector2f{
                0.25f * y(0) - 0.01f * y(0) * y(1),
                -y(1) + 0.01f * y(0) * y(1)
            };
        };
        auto f2 = [](const float t, const Eigen::Vector2f& y) noexcept {
            auto f1 = Eigen::Vector2f{
                0.25f * y(0) - 0.01f * y(0) * y(1),
                -y(1) + 0.01f * y(0) * y(1)
            };
            return Eigen::Vector2f{
                (0.25f - 0.01f * y(1)) * f1(0) - 0.01f * y(0) * f1(1),
                (0.01f * y(0) - 1.0f) * f1(1) + 0.01f * y(1) * f1(0)
            };
        };
        auto f3 = [](const float t, const Eigen::Vector2f& y) noexcept {
            auto f1 = Eigen::Vector2f{
                0.25f * y(0) - 0.01f * y(0) * y(1),
                -y(1) + 0.01f * y(0) * y(1)
            };
            auto f2 = Eigen::Vector2f{
                (0.25f - 0.01f * y(1)) * f1(0) - 0.01f * y(0) * f1(1),
                (0.01f * y(0) - 1.0f) * f1(1) + 0.01f * y(1) * f1(0)
            };
            return Eigen::Vector2f{
                (0.25f - 0.01f * y(1)) * f2(0) - 0.01f * y(0) * f2(1) - 0.02f * f1(0) * f1(1),
                (0.01f * y(0) - 1.0f) * f2(1) + 0.01f * y(1) * f2(0) + 0.02f * f1(0) * f1(1)
            };
        };
        return ODESystem<2>{f1, f2, f3};
    }

    static auto VanDerPol() noexcept {
        auto f1 = [](const float t, const Eigen::Vector2f& y) noexcept {
            return Eigen::Vector2f{
                y(1),
                2.0f * ((1.0f - y(0) * y(0)) * y(1) - y(0))
            };
        };
        auto f2 = [](const float t, const Eigen::Vector2f& y) noexcept {
            const auto f1 = Eigen::Vector2f{
                y(1),
                2.0f * ((1.0f - y(0) * y(0)) * y(1) - y(0))
            };
            return Eigen::Vector2f{
                f1(1),
                2.0f * ((1.0f - y(0) * y(0)) * f1(1) - (2.0f * y(0) * y(1) + 1.0f) * f1(0))
            };
        };
        auto f3 = [](const float t, const Eigen::Vector2f& y) noexcept {
            const auto f1 = Eigen::Vector2f{
                y(1),
                2.0f * ((1.0f - y(0) * y(0)) * y(1) - y(0))
            };
            const auto f2 = Eigen::Vector2f{
                f1(1),
                2.0f * ((1.0f - y(0) * y(0)) * f1(1) - (2.0f * y(0) * y(1) + 1.0f) * f1(0))
            };
            return Eigen::Vector2f{
                f2(1),
                -(2.0f * y(0) * y(1) + 1.0f) * f2(0) - 4.0f * y(0) * f1(0) * f1(1) - 2.0f * f1(0) * f1(0) * y(1) + (1.0f - y(0) * y(0)) * f2(1)
            };
        };
        return ODESystem<2>{f1, f2, f3};
    }

    template <int dim>
    static auto MethodOfLines() noexcept {
        auto f1 = [](const float t, const Eigen::Vector<float, dim>& y) noexcept {

        };
    }
} // namespace pece

#endif //ODE_H
