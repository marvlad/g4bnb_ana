#ifndef HISTOGRAM_CONFIG_H
#define HISTOGRAM_CONFIG_H

#include <string>

struct AxisConfig {
    int nbins;
    double min;
    double max;
};

struct TH1Config {
    std::string name;
    std::string title;
    AxisConfig x;
};

struct TH2Config {
    std::string name;
    std::string title;
    AxisConfig x;
    AxisConfig y;
};

#endif
