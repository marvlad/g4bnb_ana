#ifndef HISTGROUP_H
#define HISTGROUP_H

#include <vector>
#include <string>
#include <stdexcept>

#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>

template <typename H>
class HistGroup {
public:
    HistGroup(const std::string& prefix,
              int nGroups,
              const std::string& title,
              int nx, double x1, double x2,
              int ny = 1, double y1 = 0, double y2 = 1,
              int nz = 1, double z1 = 0, double z2 = 1);

    void Fill(int groupIndex, double x, double y = 0, double z = 0, double w = 1.0);
    void FillAll(double x, double y = 0, double z = 0, double w = 1.0);

    H* Get(int i)   { return bins_.at(i); }
    H* GetAll()     { return all_; }

    /// NEW: number of groups
    int size() const { return nGroups_; }

    void Write();

private:
    std::string prefix_;
    std::string title_;

    std::vector<H*> bins_;
    H* all_;

    int nGroups_;          ///< NEW — store number of histograms

    int nx_; double x1_, x2_;
    int ny_; double y1_, y2_;
    int nz_; double z1_, z2_;
};

// -----------------------------------------------------------------------------
// Implementation
// -----------------------------------------------------------------------------

template <typename H>
HistGroup<H>::HistGroup(const std::string& prefix,
                        int nGroups,
                        const std::string& title,
                        int nx, double x1, double x2,
                        int ny, double y1, double y2,
                        int nz, double z1, double z2)
    : prefix_(prefix), title_(title),
      nGroups_(nGroups),            // <--- store it!
      nx_(nx), x1_(x1), x2_(x2),
      ny_(ny), y1_(y1), y2_(y2),
      nz_(nz), z1_(z1), z2_(z2)
{
    // Create N group histograms
    for (int i = 0; i < nGroups_; i++) {
        std::string name = prefix_ + "_bin" + std::to_string(i);

        if constexpr (std::is_same<H, TH1F>::value)
            bins_.push_back(new TH1F(name.c_str(), title_.c_str(), nx_, x1_, x2_));

        else if constexpr (std::is_same<H, TH2F>::value)
            bins_.push_back(new TH2F(name.c_str(), title_.c_str(),
                                     nx_, x1_, x2_, ny_, y1_, y2_));

        else if constexpr (std::is_same<H, TH3F>::value)
            bins_.push_back(new TH3F(name.c_str(), title_.c_str(),
                                     nx_, x1_, x2_, ny_, y1_, y2_,
                                     nz_, z1_, z2_));
    }

    // Create "all" histogram
    std::string allname = prefix_ + "_all";

    if constexpr (std::is_same<H, TH1F>::value)
        all_ = new TH1F(allname.c_str(), title_.c_str(), nx_, x1_, x2_);

    else if constexpr (std::is_same<H, TH2F>::value)
        all_ = new TH2F(allname.c_str(), title_.c_str(),
                        nx_, x1_, x2_, ny_, y1_, y2_);

    else if constexpr (std::is_same<H, TH3F>::value)
        all_ = new TH3F(allname.c_str(), title_.c_str(),
                        nx_, x1_, x2_, ny_, y1_, y2_, nz_, z1_, z2_);
}

template <typename H>
void HistGroup<H>::Fill(int idx, double x, double y, double z, double w) {
    if (idx < 0 || idx >= nGroups_)   // <--- now using nGroups_
        throw std::out_of_range("Invalid histogram subgroup index");

    if constexpr (std::is_same<H, TH1F>::value)
        bins_[idx]->Fill(x, w);
    else if constexpr (std::is_same<H, TH2F>::value)
        bins_[idx]->Fill(x, y, w);
    else if constexpr (std::is_same<H, TH3F>::value)
        bins_[idx]->Fill(x, y, z, w);
}

template <typename H>
void HistGroup<H>::FillAll(double x, double y, double z, double w) {
    if constexpr (std::is_same<H, TH1F>::value)
        all_->Fill(x, w);
    else if constexpr (std::is_same<H, TH2F>::value)
        all_->Fill(x, y, w);
    else if constexpr (std::is_same<H, TH3F>::value)
        all_->Fill(x, y, z, w);
}

template <typename H>
void HistGroup<H>::Write() {
    for (auto h : bins_) h->Write();
    all_->Write();
}

#endif

