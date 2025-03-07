//
// Created by Leon Kang on 2025/3/7.
//

#ifndef PLOT_H
#define PLOT_H

#include <nlohmann/json.hpp>
#include <string>

namespace pece {
    struct Plot {
        std::string label;
        std::span<const float> xs;
        std::span<const float> ys;
    };

    class Figure {
        std::string _title;
        std::vector<Plot> _plots;
    public:
        explicit Figure(const std::string_view title) : _title(title) {}
        void AddPlot(const std::string_view label, const std::span<const float> xs, const std::span<const float> ys) {
            _plots.emplace_back(Plot {std::string{label}, xs, ys});
        }
        [[nodiscard]] auto ToString() const noexcept {
            nlohmann::json figure, plots;
            for (const auto&[label, xs, ys] : _plots) {
                plots.emplace_back(nlohmann::json{
                    {"label", label},
                    {"xs", xs},
                    {"ys", ys}
                });
            }
            figure["title"] = _title;
            figure["plots"] = plots;
            return figure.dump(4);
        }
    };
} // namespace pece

#endif //PLOT_H
