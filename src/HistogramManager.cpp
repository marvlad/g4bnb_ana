#include "HistogramManager.h"

HistogramManager::HistogramManager(const HistGroupConfig& hgConf,
                                   const TH2Config& th2AConf,
                                   const TH2Config& th2BConf)
{
    // Create 1 grouped histogram set
    HistoGroup = new HistGroup<TH1F>(
        hgConf.prefix,
        hgConf.nGroups,
        hgConf.title,
        hgConf.nx, hgConf.x1, hgConf.x2
    );

    // Create the standalone TH2 histograms
    H2A = MakeTH2F(th2AConf);
    H2B = MakeTH2F(th2BConf);
}

HistogramManager::~HistogramManager()
{
    delete HistoGroup;
    delete H2A;
    delete H2B;
}

TH2F* HistogramManager::MakeTH2F(const TH2Config& cfg)
{
    return new TH2F(
        cfg.name.c_str(),
        cfg.title.c_str(),
        cfg.nx, cfg.x1, cfg.x2,
        cfg.ny, cfg.y1, cfg.y2
    );
}

void HistogramManager::WriteAll()
{
    HistoGroup->Write();
    H2A->Write();
    H2B->Write();
}
