#include "../src/time_window.h"
#include <cassert>
#include <limits>
int main() {
    // Non-binary durations and a shortened final window exercise rounding.
    for (double end : {0.01, 0.3, 20., 40.}) {
        double time = 0;
        int n = 0;
        while (time < end) {
            const double endpoint = advection_endpoint(time, .037, end);
            while (time < endpoint) {
                double dt = checked_acoustic_step(time, endpoint, .0031);
                time = std::min(endpoint, time + dt);
                assert(time <= end);
                assert(++n < 100000);
            }
        }
        assert(time == end);
    }
    for (double bad : {0., -1., std::numeric_limits<double>::quiet_NaN()}) {
        bool rejected = false;
        try { checked_acoustic_step(0., 1., bad); }
        catch (const std::runtime_error &) { rejected = true; }
        assert(rejected);
    }
    bool rejected = false;
    try { checked_acoustic_step(1., 2., 1e-30); }
    catch (const std::runtime_error &) { rejected = true; }
    assert(rejected);
}
