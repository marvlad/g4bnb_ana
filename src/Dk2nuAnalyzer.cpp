#include "Dk2nuAnalyzer.h"
#include <iostream>
#include <TFile.h>
#include <cmath>

Dk2nuAnalyzer::Dk2nuAnalyzer(const std::vector<std::string>& files)
    : files_(files), chain_(nullptr), dk2nu_(nullptr)
{
    rand1_ = new TRandom3(40411199);
    rand2_ = new TRandom3(40411200);

    BookHistograms();
    SetupChain();
}

void Dk2nuAnalyzer::BookHistograms() {

    smeared_vs_Energy =
        new TH2F("smeared_vs_Energy",
                 "Smeared vs Energy;#Delta t [ns];E_{#nu} [GeV]",
                 100, -3., 20.,
                 100, 0., 6.);
}

void Dk2nuAnalyzer::SetupChain() {

    chain_ = new TChain("dk2nuTree");

    for (auto& f : files_) {
        std::cout << "Adding input file: " << f << std::endl;
        chain_->Add(f.c_str());
    }

    dk2nu_ = new bsim::Dk2Nu();

    chain_->SetBranchAddress("dk2nu", &dk2nu_);
}

double Dk2nuAnalyzer::ApplyTimeShift(double timeres, TRandom3* r, TH1D* prof, bool smear) {

    double det = 0, profs = 0;

    if (smear) {
        det = r->Gaus(0., timeres);
        if (prof) profs = prof->GetRandom() * 1e9;
    }

    return det + profs;
}

void Dk2nuAnalyzer::Run() {

    Long64_t N = chain_->GetEntries();
    std::cout << "Total entries: " << N << std::endl;

    for (Long64_t i = 0; i < N; i++) {

        chain_->GetEntry(i);

        int Nance = dk2nu_->ancestor.size();
        if (Nance < 2) continue;

        double Enu = dk2nu_->nuray[1].E;

        double startz = dk2nu_->ancestor[Nance-1].startz;
        double startT = dk2nu_->ancestor[Nance-1].startt;

        double dT = startT - ((startz + 360.) / 29.98);

        double smear = ApplyTimeShift(0.5, rand1_, nullptr);
        double dTr = dT + smear;

        smeared_vs_Energy->Fill(dTr, Enu);
    }
}

void Dk2nuAnalyzer::Save(const std::string& outname) {

    TFile* f = TFile::Open(outname.c_str(), "RECREATE");

    smeared_vs_Energy->Write();

    f->Close();

    std::cout << "Saved output --->  " << outname << std::endl;
}
