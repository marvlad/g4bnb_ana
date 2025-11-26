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

bool Dk2nuAnalyzer::IsGoodFile(const std::string& fname) const
{
    TFile f(fname.c_str(), "READ");

    if (f.IsZombie()) {
        std::cerr << "[BAD FILE] Cannot open: " << fname << std::endl;
        return false;
    }

    if (f.TestBit(TFile::kRecovered)) {
        std::cerr << "[BAD FILE] CORRUPTED (Recovered) : " << fname << std::endl;
        return false;
    }

    TTree* t = (TTree*) f.Get("dk2nuTree");
    if (!t) {
        std::cerr << "[BAD FILE] Missing dk2nuTree: " << fname << std::endl;
        return false;
    }

    if (t->GetEntries() == 0) {
        std::cerr << "[BAD FILE] Empty dk2nuTree: " << fname << std::endl;
        return false;
    }

    return true;
}

void Dk2nuAnalyzer::BookHistograms() {

    smeared_vs_Energy = new TH2F("smeared_vs_Energy",
                 "Smeared vs Energy;#Delta t [ns];E_{#nu} [GeV]",
                 100, -3., 20.,
                 100, 0., 6.);

    hEnuPion =  new HistGroup<TH1F>("h_Enu_pion", 7, ";Neutrino Energy (GeV);Counts", 40, 0, 8);
}

void Dk2nuAnalyzer::SetupChain()
{
    chain_ = new TChain("dk2nuTree");

    totalFiles_ = files_.size();
    goodFiles_ = 0;

    for (auto& f : files_) {
        if (IsGoodFile(f)) {
            std::cout << "[GOOD] Adding: " << f << std::endl;
            chain_->Add(f.c_str());
            goodFiles_++;
        }
        else {
            std::cout << "[SKIPPED] " << f << std::endl;
        }
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
        double wgt = dk2nu_->nuray[1].wgt*dk2nu_->decay.nimpwt;
        if (Nance < 2) continue;

        double Enu = dk2nu_->nuray[1].E;

        double startz = dk2nu_->ancestor[Nance-1].startz;
        double startT = dk2nu_->ancestor[Nance-1].startt;

        double dT = startT - ((startz + 360.) / 29.98);

        double smear = ApplyTimeShift(0.5, rand1_, nullptr);
        double dTr = dT + smear;

        smeared_vs_Energy->Fill(dTr, Enu);

        int idx = deltaT_to_bin(dT);

        if (idx >= 0) {
           validate_dt_hist_count(hEnuPion->size()); // optional after init
           hEnuPion->Fill(idx, Enu, 0, 0, wgt);
        }

        hEnuPion->FillAll(Enu, 0, 0, wgt);

    }
}

void Dk2nuAnalyzer::Save(const std::string& outname) {

    TFile* f = TFile::Open(outname.c_str(), "RECREATE");

    smeared_vs_Energy->Write();
    hEnuPion->Write();

    f->Close();

    std::cout << "Saved output --->  " << outname << std::endl;
    delete smeared_vs_Energy;
    delete hEnuPion;
}

void Dk2nuAnalyzer::PrintInfo() const
{
    double eff = (totalFiles_ > 0)
        ? 100.0 * goodFiles_ / totalFiles_
        : 0.0;

    std::cout << "\n========== File Quality Report ==========\n";
    std::cout << "  Total input files:  " << totalFiles_ << "\n";
    std::cout << "  Good files:         " << goodFiles_ << "\n";
    std::cout << "  Bad files:          " << totalFiles_ - goodFiles_ << "\n";
    std::cout << "  Efficiency:         " << eff << " %\n";
    std::cout << "=========================================\n\n";
}
