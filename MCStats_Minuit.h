#ifndef __MCSTATS_MINUIT__
#define __MCSTATS_MINUIT__

#include <algorithm>
#include <iostream>
#include <vector>
#include <string>

#include "TTree.h"
#include "TRandom3.h"
#include "TFile.h"
#include "TF1.h"
#include "TH1D.h"
#include "TMath.h"

#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
 
enum test_type {
  // Use only statistical Poisson likelihood
  LLHStat_StatOnly, 

  // Use Pearson statistic, i.e. assuming Gaussian variation with Poisson error on data only
  LLHStat_PearsonNoPen, 

  // Use Pearson statistic, i.e. assuming Gaussian variation with Poisson error on both data and MC 
  LLHStat_PearsonPen, 

  // Use statistical Poisson likelihood by fluctuating the MC around the MC central value and using that for the MC in the likelihood
  //LLHStat_StatWithFluc, 

  // Use statistical Poisson likelihood and introduce a normalisation penalty relative the uncertainty in the bin
  //LLHStat_StatWithNorm, 

  // Barlow-Beeston with Gaussian
  LLH_BarlowBeestonGauss,

  // Barlow-Beeston with Poisson
  //LLH_BarlowBeestonPoisson,

  // Tianlu LLH
  LLH_Tianlu,

  nLLH
};

// A class to hold information about the fit
class storer {
  public:
    storer(int);
    ~storer();
    void Setup();
    void ReWeight();
    double GetLikelihood();
    double DoEval(const double *x) ;

    TH1D *GetData() { return hData; };
    TH1D *GetMC() { return hMCnorm; };
    TH1D *GetUnscaledMC() { return hMC; };

    bool CheckMe();

    void SetMean(double meanp) { mean = meanp; }
    void SetSigma(double sigmap) { sigma = sigmap; }
    void SetNorm(double normp) { norm = normp; }
    void SetMCfact(double mcfactorp) { MCfactor = mcfactorp; }

    void SetTestType(test_type type) { testtype = type; }

    double GetRefLLH() { return RefLLH; }

  private:
    TH1D* hData;
    TH1D* hMC;
    TH1D* hMCnorm;
    TF1* bkg;
    TF1* sig;
    test_type testtype;

    double mean;
    double sigma;
    double norm;
    double MCfactor;

    // Proposed
    double meanp;
    double sigmap;

    // Reference LLH
    double RefLLH;

    // The raw vector of MC
    std::vector<double> MC;
    std::vector<double> binnorm;
    std::vector<double> binuncert;

    TRandom3 *rnd;
};
#endif
