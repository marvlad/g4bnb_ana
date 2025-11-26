#ifndef TIME_BINNING_H
#define TIME_BINNING_H

#include <vector>
#include <utility>
#include <stdexcept>
#include <iostream>

/// DElta_t bin definitions
static const std::vector<std::pair<double,double>> dt_bins = {
    {-0.3, -0.1},
    {-0.1,  0.1},
    { 0.1,  0.3},
    { 0.3,  0.5},
    { 0.5,  0.7},
    { 0.7,  0.9},
    { 0.9,  1.1}
};

/// Return total number of Δt bins
inline size_t deltaT_bin_count() {
    return dt_bins.size();
}

/// Convert Delta_t to a bin index in [0, N-1], or -1 if outside all bins.
inline int deltaT_to_bin(double dT) {
    for (size_t i = 0; i < dt_bins.size(); ++i) {
        if (dT >= dt_bins[i].first && dT <= dt_bins[i].second)
            return static_cast<int>(i);
    }
    return -1;  // No valid bin
}

/// Validate that user histogram vector matches required size.
/// Will throw (or warn) the first time it is called incorrectly.
inline void validate_dt_hist_count(size_t provided) {
    const size_t required = dt_bins.size();
    if (provided != required) {
        std::cerr 
            << "[TimeBinning ERROR] Wrong number of Delta_t histograms!  Provided "
            << provided << ", but expected exactly " << required << ".\n"
            << "This analysis assumes " << required 
            << " time bins defined in TimeBinning.h.\n";
        throw std::runtime_error("Incorrect Δt histogram count");
    }
}
#endif

