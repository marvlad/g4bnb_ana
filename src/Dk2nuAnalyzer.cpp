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
                 "Smeared vs Energy;E_{#nu} [GeV];#Delta t t[ns]",
                 100, -3., 20.,
                 100, 0., 1.);

    std::string histLabel = ";Neutrino Energy (GeV);Counts";
    std::string hist2DLabel = "Energy [GeV];#Delta t [ns]";

    // Numu + Pion
    HistogramManager::HistGroupConfig hEnumuPionCfg = {"h_Enumu_pion", 7, histLabel, 40, 0, 8};
    HistogramManager::TH2Config hEnumuPionh2Acfg = {"h_deltat_vs_energy_numu_pion", ";Neutrino " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    HistogramManager::TH2Config hEnumuPionh2Bcfg = { "h_deltat_vs_energy_parent_pion", ";Pion " + hist2DLabel,  100, 0, 10, 100, 0, 1};
    hEnumuPion = new HistogramManager(hEnumuPionCfg, hEnumuPionh2Acfg, hEnumuPionh2Bcfg);

    // Nue + Pion
    HistogramManager::HistGroupConfig hEnuePionCfg = {"h_Enue_pion", 7, histLabel, 40, 0, 8};
    HistogramManager::TH2Config hEnuePionh2Acfg = {"h_deltat_vs_energy_nue_pion", ";Neutrino " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    HistogramManager::TH2Config hEnuePionh2Bcfg = { "h_deltat_vs_energy_parent_pion_e", ";Pion " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    hEnuePion = new HistogramManager(hEnuePionCfg, hEnuePionh2Acfg, hEnuePionh2Bcfg);

    // Numu + Kaon
    HistogramManager::HistGroupConfig hEnumuKaonCfg = {"h_Enumu_kaon", 7, histLabel, 40, 0, 8};
    HistogramManager::TH2Config hEnumuKaonh2Acfg = {"h_deltat_vs_energy_numu_kaon", ";Neutrino " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    HistogramManager::TH2Config hEnumuKaonh2Bcfg = { "h_deltat_vs_energy_parent_kaon", ";Kaon " + hist2DLabel,  100, 0, 10, 100, 0, 1};
    hEnumuKaon = new HistogramManager(hEnumuKaonCfg, hEnumuKaonh2Acfg, hEnumuKaonh2Bcfg);

    // Numu + Kaon_e
    HistogramManager::HistGroupConfig hEnueKaonCfg = {"h_Enue_kaon", 7, histLabel, 40, 0, 8};
    HistogramManager::TH2Config hEnueKaonh2Acfg = {"h_deltat_vs_energy_nue_kaon", ";Neutrino " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    HistogramManager::TH2Config hEnueKaonh2Bcfg = { "h_deltat_vs_energy_parent_kaon_e", ";Kaon " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    hEnueKaon = new HistogramManager(hEnueKaonCfg, hEnueKaonh2Acfg, hEnueKaonh2Bcfg);

    // Numu all
    HistogramManager::HistGroupConfig hEnumuAllCfg = {"h_Enumu_all", 7, histLabel, 40, 0, 8};
    HistogramManager::TH2Config hEnumuAllh2Acfg = {"h_deltat_vs_energy_numu_all", ";Neutrino " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    HistogramManager::TH2Config hEnumuAllh2Bcfg = { "h_deltat_vs_energy_parent_numu_all", ";All " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    hEnumuAll = new HistogramManager(hEnumuAllCfg, hEnumuAllh2Acfg, hEnumuAllh2Bcfg);

    // Nue all
    HistogramManager::HistGroupConfig hEnueAllCfg = {"h_Enue_all", 7, histLabel, 40, 0, 8};
    HistogramManager::TH2Config hEnueAllh2Acfg = {"h_deltat_vs_energy_nue_all", ";Neutrino " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    HistogramManager::TH2Config hEnueAllh2Bcfg = { "h_deltat_vs_energy_parent_nue_all", ";All " + hist2DLabel, 100, 0, 10, 100, 0, 1};
    hEnueAll = new HistogramManager(hEnueAllCfg, hEnueAllh2Acfg, hEnueAllh2Bcfg);

}

void Dk2nuAnalyzer::SetupChain()
{
    chain_ = new TChain("dk2nuTree");

    totalFiles_ = files_.size();
    goodFiles_ = 0;

    for (auto& f : files_) {
        if (IsGoodFile(f)) {
            if(verbosity==2) std::cout << "[GOOD] Adding: " << f << std::endl;
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

void Dk2nuAnalyzer::PrintDk2nu(int eventid){
    std::cout << "\n =========================> Start of Event ID: " << eventid << std::endl;
    int Nance = dk2nu_->ancestor.size();
    std::cout << "Total ancestors: " << Nance << std::endl;
    for (int i=0; i<Nance; i++){
        std::cout << "Index: " << i << ", PDG: " << dk2nu_->ancestor[i].pdg << std::endl;
        std::cout << "\t Start x,y,z,t: " << dk2nu_->ancestor[i].startx<<","<<dk2nu_->ancestor[i].starty<<","<<dk2nu_->ancestor[i].startz<<","<<dk2nu_->ancestor[i].startt << std::endl;
        //cout << "Stop x,y,z: "<<dk2nu_->ancestor[i].stopx<<","<<dk2nu_->ancestor[i].stopy<<","<<dk2nu_->ancestor[i].stopz<<std::endl;
        std::cout << "\t Momentum start px,py,pz: "<<dk2nu_->ancestor[i].startpx<<", "<<dk2nu_->ancestor[i].startpy<<", "<<dk2nu_->ancestor[i].startpz << std::endl;
        std::cout << "\t Momentum stop px,py,pz: "<<dk2nu_->ancestor[i].stoppx<<", "<<dk2nu_->ancestor[i].stoppy<<", "<<dk2nu_->ancestor[i].stoppz<<std::endl;
        std::cout << "\t Proc: "<<dk2nu_->ancestor[i].proc<<std::endl;
        std::cout << "\t imat: "<<dk2nu_->ancestor[i].imat<<std::endl;
        std::cout << "\t ivol: "<<dk2nu_->ancestor[i].ivol<<std::endl;
        std::cout << "\t nucleus: "<<dk2nu_->ancestor[i].nucleus<<std::endl;
        std::cout << "\t Polarization: x,y,z "<<dk2nu_->ancestor[i].polx<<", "<<dk2nu_->ancestor[i].poly<<", "<<dk2nu_->ancestor[i].polz<<std::endl;
        //std::cout << "\t Ray P: x,y,z "<<dk2nu_->nuray[i+2].px<<", "<<dk2nu_->nuray[i].py<<", "<<dk2nu_->nuray[i].pz<<std::endl;
    }
    std::cout << "decay.xyz "<<dk2nu_->decay.vx<<", "<<dk2nu_->decay.vy<<", "<<dk2nu_->decay.vz<<std::endl; 
}

void Dk2nuAnalyzer::Run() {

    int nparticles,nutype,nuparent,mpotnum;
    double nustartx,nustarty,nustartz,nustartT;
    double nuE, bwt0, bwt1, bwt2, decay_wt;
    double nupx,nupy,nupz,nupperp,nuang;
    double parentpx,parentpy,parentpz,parentpperp;
    double parentE,parentang;
    double dT,dTr;
 
    bwt0=0.; bwt1=0.; bwt2=0.; decay_wt=0.;

    double nuendx = 0.; 
    double nuendy = 0.; 
    double nuendz = 0.; 

    // -----
    Long64_t N = chain_->GetEntries();
    std::cout << "Total entries: " << N << std::endl;

    int three_body_decay = 0;

    for (Long64_t i = 0; i < N; i++) {
        chain_->GetEntry(i);
        int Nance = dk2nu_->ancestor.size();
  
        //mpotnum = Npot;
        nparticles= Nance;
  
        nutype = dk2nu_->ancestor[Nance-1].pdg;
        nuparent = dk2nu_->ancestor[Nance-2].pdg;
    
        //Access to the origin (x,y,z,t)
        nustartx = dk2nu_->ancestor[Nance-1].startx;
        nustarty = dk2nu_->ancestor[Nance-1].starty;
        nustartz = dk2nu_->ancestor[Nance-1].startz;
        nustartT = dk2nu_->ancestor[Nance-1].startt;

        nupx = dk2nu_->nuray[1].px;
        nupy = dk2nu_->nuray[1].py;
        nupz = dk2nu_->nuray[1].pz;
    
        nupperp = sqrt(nupx*nupx+nupy*nupy);
        nuang = TMath::ATan(nupperp/nupz);
    
        parentpx = dk2nu_->ancestor[Nance-2].startpx;
        parentpy = dk2nu_->ancestor[Nance-2].startpy;
        parentpz = dk2nu_->ancestor[Nance-2].startpz;
        parentpperp = sqrt(parentpx*parentpx+parentpy*parentpy);
        parentE = dk2nu_->decay.ppenergy;
        parentang = TMath::ATan(parentpperp/parentpz);
    
        nuE = dk2nu_->nuray[1].E;
    
        decay_wt = dk2nu_->decay.nimpwt;
    
        bwt0 = dk2nu_->nuray[0].wgt*dk2nu_->decay.nimpwt;
        bwt1 = dk2nu_->nuray[1].wgt*dk2nu_->decay.nimpwt;
        bwt2 = dk2nu_->nuray[2].wgt*dk2nu_->decay.nimpwt;
    
        dT = nustartT - ((nustartz+360.)/29.98);
    
        bool smearing = true;  // Enable or disable smearing as needed
        double timeresolution = 0.5; // Adjust as needed
    
        // Apply time shift with smearing
        TH1D* profileShift = NULL;
        double smearedShift = ApplyTimeShift(0.5, rand1_, profileShift); // Adjust timeresolution as needed
        dTr = dT;// + smearedShift;
    
        //smeared_vs_Energy->Fill(nuE, dTr, bwt1);
    
        //int Nance = dk2nu_->ancestor.size();
        //if(Nance == 3)
        //continue;
    
        double startx1 = dk2nu_->ancestor[1].startx;
        double starty1 = dk2nu_->ancestor[1].starty;
        double startz1 = dk2nu_->ancestor[1].startz;
     
        double startxf = dk2nu_->ancestor[Nance - 1].startx;
        double startyf = dk2nu_->ancestor[Nance - 1].starty;
        double startzf = dk2nu_->ancestor[Nance - 1].startz;
    
	TVector3 startParent(startx1,starty1,startz1);
	TVector3 startGrandparent(startxf,startyf,startzf);
        
	double tracklength = dk2nu_->ancestor[1].tracklength/10.;
     
        double decayt = dk2nu_->ancestor[Nance - 1].startt - dk2nu_->ancestor[1].startt;
    
        double E_parent = dk2nu_->decay.ppenergy;
        double Enu = dk2nu_->nuray[1].E;
        double wgt = dk2nu_->nuray[1].wgt*dk2nu_->decay.nimpwt;
        double c = 29.9792; //cm/ns 
     
        if (verbosity==2) PrintDk2nu(i);
        smeared_vs_Energy->Fill(dTr, Enu);
        
        if(Nance == 3) three_body_decay++;

        double deltaz1 = startzf - startz1;
	
        bool fullg4tracklenght = true;    
        double length3d;
        if(fullg4tracklenght) length3d = tracklength;
        else length3d = (startParent - startGrandparent).Mag();

        bool length1d = false;
        double deltat;
	if(length1d) deltat = decayt - (deltaz1/c);
	else deltat = decayt - (length3d/c);

        //h_deltat_vs_energy_numu_pion->Fill(Enu,deltat,wgt);
        //h_deltat_vs_energy_parent_pion->Fill(E_parent,deltat,wgt);
        //h_Enu_pion_all_dis->Fill((old_deltat?(startzf - startz1):distance(startx1,starty1,startz1,startxf,startyf,startzf)) - tracklength);
        //h_Enu_pion_all->Fill(Enu,wgt);

        smeared_vs_Energy->Fill(nuE, deltat, bwt1);

        // Pion + nu_mu 
        // =========================
        if(dk2nu_->ancestor[1].pdg == 211 && dk2nu_->decay.ntype == 14 && Nance == 3){
          if(verbosity==2) std::cout << "PDG = " << dk2nu_->ancestor[1].pdg << ", corresponding to pions " << std::endl;
          int idx = deltaT_to_bin(deltat);
          if (idx >= 0) {
             validate_dt_hist_count(hEnumuPion->HistoGroup->size()); // optional after init
             hEnumuPion->HistoGroup->Fill(idx, Enu, 0, 0, wgt);
          }
          hEnumuPion->HistoGroup->FillAll(Enu, 0, 0, wgt);
          hEnumuPion->H2A->Fill(Enu,deltat,wgt);
          hEnumuPion->H2B->Fill(E_parent,deltat,wgt);
        }

        // Pion + nu_e 
        // =========================
        if(dk2nu_->ancestor[1].pdg == 211 && dk2nu_->decay.ntype == 12 && Nance == 3){
          if(verbosity==2) std::cout << "PDG = " << dk2nu_->ancestor[1].pdg << ", corresponding to pions from nue " << std::endl;
          int idx = deltaT_to_bin(deltat);
          if (idx >= 0) {
             validate_dt_hist_count(hEnuePion->HistoGroup->size()); // optional after init
             hEnuePion->HistoGroup->Fill(idx, Enu, 0, 0, wgt);
          }
          hEnuePion->HistoGroup->FillAll(Enu, 0, 0, wgt);
          hEnuePion->H2A->Fill(Enu,deltat,wgt);
          hEnuePion->H2B->Fill(E_parent,deltat,wgt);
        }

        // Kaon + nu_mu 
        // =========================
        if(dk2nu_->ancestor[1].pdg == 321 && dk2nu_->decay.ntype == 14 && Nance == 3){
          if(verbosity==2) std::cout << "PDG = " << dk2nu_->ancestor[1].pdg << ", corresponding to kaons " << std::endl;
          int idx = deltaT_to_bin(deltat);
          if (idx >= 0) {
             validate_dt_hist_count(hEnumuKaon->HistoGroup->size()); // optional after init
             hEnumuKaon->HistoGroup->Fill(idx, Enu, 0, 0, wgt);
          }
          hEnumuKaon->HistoGroup->FillAll(Enu, 0, 0, wgt);
          hEnumuKaon->H2A->Fill(Enu,deltat,wgt);
          hEnumuKaon->H2B->Fill(E_parent,deltat,wgt);
        }

        // Kaon + nu_e
        // =========================
        if(dk2nu_->ancestor[1].pdg == 321 && dk2nu_->decay.ntype == 12 && Nance == 3){
          if(verbosity==2) std::cout << "PDG = " << dk2nu_->ancestor[1].pdg << ", corresponding to kaons " << std::endl;
          int idx = deltaT_to_bin(deltat);
          if (idx >= 0) {
             validate_dt_hist_count(hEnueKaon->HistoGroup->size()); // optional after init
             hEnueKaon->HistoGroup->Fill(idx, Enu, 0, 0, wgt);
          }
          hEnueKaon->HistoGroup->FillAll(Enu, 0, 0, wgt);
          hEnueKaon->H2A->Fill(Enu,deltat,wgt);
          hEnueKaon->H2B->Fill(E_parent,deltat,wgt);
        }

        // Numu all 
        // =========================
        if(dk2nu_->decay.ntype == 14 && Nance == 3){
          if(verbosity==2) std::cout << "PDG = " << dk2nu_->ancestor[1].pdg << ", corresponding numu only " << std::endl;
          int idx = deltaT_to_bin(deltat);
          if (idx >= 0) {
             validate_dt_hist_count(hEnumuAll->HistoGroup->size()); // optional after init
             hEnumuAll->HistoGroup->Fill(idx, Enu, 0, 0, wgt);
          }
          hEnumuAll->HistoGroup->FillAll(Enu, 0, 0, wgt);
          hEnumuAll->H2A->Fill(Enu,deltat,wgt);
          hEnumuAll->H2B->Fill(E_parent,deltat,wgt);
        }

        // Nue all 
        // =========================
        if(dk2nu_->decay.ntype == 12 && Nance == 3){
          if(verbosity==2) std::cout << "PDG = " << dk2nu_->ancestor[1].pdg << ", corresponding to nue only " << std::endl;
          int idx = deltaT_to_bin(deltat);
          if (idx >= 0) {
             validate_dt_hist_count(hEnueAll->HistoGroup->size()); // optional after init
             hEnueAll->HistoGroup->Fill(idx, Enu, 0, 0, wgt);
          }
          hEnueAll->HistoGroup->FillAll(Enu, 0, 0, wgt);
          hEnueAll->H2A->Fill(Enu,deltat,wgt);
          hEnueAll->H2B->Fill(E_parent,deltat,wgt);
        }

    }

    std::cout << "Three body decay: " << three_body_decay << std::endl;
}

void Dk2nuAnalyzer::Save(const std::string& outname) {

    TFile* f = TFile::Open(outname.c_str(), "RECREATE");

    smeared_vs_Energy->Write();
    hEnumuPion->WriteAll();
    hEnumuKaon->WriteAll();
    hEnuePion->WriteAll();
    hEnueKaon->WriteAll();
    hEnumuAll->WriteAll();
    hEnueAll->WriteAll();

    f->Close();

    std::cout << "Saved output --->  " << outname << std::endl;
    delete smeared_vs_Energy;
    delete hEnumuPion;
    delete hEnumuKaon;
    delete hEnuePion;
    delete hEnueKaon;
    delete hEnumuAll;
    delete hEnueAll;
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

