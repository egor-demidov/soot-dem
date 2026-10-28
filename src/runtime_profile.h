#pragma once
#include <array>
#include <chrono>
#include <cstdlib>
#include <iomanip>
#include <ostream>

// Opt-in wall-clock counters. All scope boundaries are on the calling thread;
// threaded kernels finish before their scope exits. Nested categories overlap.
class CapillaryRuntimeProfile {
public:
    enum Stage { Momentum, ContactLine, Projection, DetectCells, Remap, WallVelocity, WallLoad, Count };
    using Clock = std::chrono::steady_clock;
    bool enabled = std::getenv("CAPILLARY_PROFILE") != nullptr;
    std::array<double, Count> seconds{};
    std::array<unsigned long long, Count> calls{};
    struct Scope {
        CapillaryRuntimeProfile& profile;
        Stage stage;
        Clock::time_point start;
        Scope(CapillaryRuntimeProfile& p, Stage s) : profile(p), stage(s) {
            if (profile.enabled) start = Clock::now();
        }
        ~Scope() {
            if (profile.enabled) {
                profile.seconds[stage] += std::chrono::duration<double>(Clock::now()-start).count();
                ++profile.calls[stage];
            }
        }
    };
    void reset() { seconds.fill(0); calls.fill(0); }
    void print(std::ostream& out, double stepSeconds, unsigned long long steps) const {
        if (!enabled) return;
        static constexpr const char* names[] = {"momentum", "contact_line", "projection", "detect_cells", "remap_inclusive", "wall_velocity", "wall_load"};
        const auto oldFlags=out.flags();
        const auto oldPrecision=out.precision();
        out << std::setprecision(9) << "Runtime profile: steps=" << steps << " seconds=" << stepSeconds << " seconds_per_step=" << stepSeconds/steps << '\n';
        for (int i=0;i<Count;++i) out << "Runtime stage: " << names[i] << " seconds=" << seconds[i] << " calls=" << calls[i] << '\n';
        out.flags(oldFlags);out.precision(oldPrecision);
    }
};
