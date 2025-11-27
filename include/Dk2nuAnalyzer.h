#ifndef DK2NU_ANALYZER_H
#define DK2NU_ANALYZER_H

#include <vector>
#include <string>
#include <TChain.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>
#include <TMath.h> 
#include <TVector3.h>

#include "dk2nu/tree/dk2nu.h"
#include "HistGroup.h"
#include "TimeBinning.h"
#include "HistogramManager.h"

class Dk2nuAnalyzer {
public:
    explicit Dk2nuAnalyzer(const std::vector<std::string>& files);
    void SetVerbosity(int v) { verbosity_ = v; }
    int  GetVerbosity() const { return verbosity_; }

    void Run();
    void Save(const std::string& outname);
    void PrintInfo() const;

private:
    int verbosity_ = 1;  // 0 = silent, 1 = normal, 2 = verbose, 3 = debug
    void SetupChain();
    void BookHistograms();
    double ApplyTimeShift(double timeres, TRandom3* r, TH1D* prof, bool smear=true);
    bool IsGoodFile(const std::string& fname) const;
    void PrintDk2nu(int eventid);


    int totalFiles_ = 0;
    int goodFiles_ = 0;

    std::vector<std::string> files_;
    int verbosity;

    TChain* chain_;
    bsim::Dk2Nu* dk2nu_;

    TRandom3* rand1_;
    TRandom3* rand2_;

    // Example histogram (you will add ALL yours)
    TH2F* smeared_vs_Energy;
    HistogramManager* hEnumuPion;
    HistogramManager* hEnumuKaon;
    HistogramManager* hEnuePion;
    HistogramManager* hEnueKaon;
    HistogramManager* hEnumuAll;
    HistogramManager* hEnueAll;
};
#endif
