#ifndef HISTOGRAM_MANAGER_H
#define HISTOGRAM_MANAGER_H

#include <string>
#include <TH2F.h>
#include "HistGroup.h"

class HistogramManager {
public:

    struct HistGroupConfig {
        std::string prefix;
        int nGroups;
        std::string title;
        int nx; double x1, x2;
    };

    struct TH2Config {
        std::string name;
        std::string title;
        int nx; double x1, x2;
        int ny; double y1, y2;
    };

    HistogramManager(const HistGroupConfig& hgConf,
                     const TH2Config& th2AConf,
                     const TH2Config& th2BConf);

    ~HistogramManager();

    HistGroup<TH1F>* HistoGroup;   // one grouped 1D histogram
    TH2F* H2A;                     // first standalone 2D histogram
    TH2F* H2B;                     // second standalone 2D histogram

    void WriteAll();

private:
    TH2F* MakeTH2F(const TH2Config& cfg);
};

#endif
