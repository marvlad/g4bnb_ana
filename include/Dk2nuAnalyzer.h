#ifndef DK2NU_ANALYZER_H
#define DK2NU_ANALYZER_H

#include <vector>
#include <string>
#include <TChain.h>
#include <TRandom3.h>
#include <TH1.h>
#include <TH2.h>

#include "dk2nu/tree/dk2nu.h"

class Dk2nuAnalyzer {
public:
    explicit Dk2nuAnalyzer(const std::vector<std::string>& files);

    void Run();
    void Save(const std::string& outname);

private:
    void SetupChain();
    void BookHistograms();
    double ApplyTimeShift(double timeres, TRandom3* r, TH1D* prof, bool smear=true);

    std::vector<std::string> files_;

    TChain* chain_;
    bsim::Dk2Nu* dk2nu_;

    TRandom3* rand1_;
    TRandom3* rand2_;

    // Example histogram (you will add ALL yours)
    TH2F* smeared_vs_Energy;
};

#endif
