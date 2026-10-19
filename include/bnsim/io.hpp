// Reproducible, explicitly labeled outputs, GPL-2.0.
#pragma once
#include "simulation.hpp"
#include <filesystem>

namespace BNSim {
class RunWriter {
public:
    RunWriter(const std::filesystem::path& parent, const Universe& initial,
              double requested_duration, double requested_timestep, const std::string& scenario);
    const std::filesystem::path& directory() const { return directory_; }
    void snapshot(const Universe& universe) const;
private:
    std::filesystem::path directory_;
    void metadata(const Universe&, double duration, double timestep, const std::string& scenario) const;
};
} // namespace BNSim
