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

enum test_type {
  // Use only statistical Poisson likelihood
  LLHStat_StatOnly, 

  // Use only statistical Poisson likelihood with a MC statistics penalty
  LLHStat_StatOnlyPen, 

  // Use Pearsons statistic, i.e. assuming Gaussian variation with Poisson error on data only
  LLHStat_PearsonNoPen, 

  // Use Pearson statistic, i.e. assuming Gaussian variation with Poisson error on both data and MC 
  LLHStat_PearsonPen, 

  // Use statistical Poisson likelihood by fluctuating the MC around the MC central value and using that for the MC in the likelihood
  LLHStat_StatWithFluc, 

  // Use statistical Poisson likelihood and introduce a normalisation penalty relative the uncertainty in the bin
  LLHStat_StatWithNorm, 

  // Barlow-Beeston with Gaussian
  LLH_BarlowBeestonGauss,

  // Barlow-Beeston with Poisson
  LLH_BarlowBeestonPoisson
};

// Standard LLH calculator
// Add in MC statistical penalty
double calcLLH(TH1D* hdata, TH1D* hmc, TH1D *herrors = NULL, bool penalty = false) {
  double LLH=0;
  for(int i=1; i<=hdata->GetNbinsX(); i++) {
    double dat=hdata->GetBinContent(i);
    double mc=hmc->GetBinContent(i);
    if(dat>0) LLH+=mc-dat+dat*TMath::Log(dat/mc);
    else LLH+=mc;
    
    if (penalty && herrors->GetBinContent(i) > 0) LLH += 1./herrors->GetBinError(i);
  }
  return LLH;
}

// Simple LLH calculator with a random variation on the bin content
// herrors is the nominal MC histogram
double calcLLHwithvar(TH1D* hdata, TH1D* hmc, TH1D* herrors, TRandom3* rnd) {
  double LLH=0;
  for(int i=1; i<=hdata->GetNbinsX(); i++) {
    double dat=hdata->GetBinContent(i);
    // Throw the MC within the fractional MC statistical error
    double mc;
    if (herrors->GetBinContent(i) > 0) mc = rnd->Gaus(hmc->GetBinContent(i), herrors->GetBinError(i)/herrors->GetBinContent(i)*hmc->GetBinContent(i));
    else mc = hmc->GetBinContent(i);

    if (dat > 0 && mc > 0) LLH += mc-dat+dat*TMath::Log(dat/mc);
    else LLH += mc;
  }
  return LLH;
}


double calcLLHwithnorm(TH1D* hdata, TH1D* hmc,std::vector<double> norm, std::vector<double> uncert) {
  double LLH=0;
  for(int i=1; i<=hdata->GetNbinsX(); i++) {
    double dat=hdata->GetBinContent(i);
    double mc=hmc->GetBinContent(i);
    if(dat>0) LLH+=mc-dat+dat*TMath::Log(dat/mc);
    else LLH+=mc;
    if (uncert[i-1] > 0) LLH+=0.5*pow((norm[i-1]-1)/uncert[i-1],2);
  }
  return LLH;
}

// Likelihood a la Pearson, which assumes a Gaussian distributed statistic with an error of sqrt(N) --- i.e. its error is Poissonian
// withPen = true penalises on both data and MC statistics
// withPen = false penalises on data only
double calcLLHwithPen(TH1D* hdata, TH1D* hmc, TH1D *herrors, bool withPen = true) {
  double LLH=0;
  for(int i=1; i<=hdata->GetNbinsX(); i++) {
    double dat=hdata->GetBinContent(i);
    double mc=hmc->GetBinContent(i);
    double mcraw = herrors->GetBinContent(i+1);
    if (withPen && dat+mc>0) LLH += (dat-mc)*(dat-mc)/(dat+mcraw);
    else if (!withPen) {
      if (dat != 0 && dat > 0) LLH += (dat-mc)*(dat-mc)/(dat);
      else if (mc > 0) LLH += (dat-mc)*(dat-mc)/(mcraw);
    }
  }
  return LLH;
}

void reweight(TH1D* hmcpred, std::vector<double> MC, TF1* nom, TF1* pred) {
  // Reset the Monte Carlo prediction
  hmcpred->Reset();
  // Reweight each MC event
  for (int i = 0; i < MC.size(); i++) {
    double preded = pred->Eval(MC[i]);
    double nomed = nom->Eval(MC[i]);
    double fill = preded/nomed;
    // Cap at 10
    //if (fill > 100) fill = 100;
    hmcpred->Fill(MC[i], fill);
  }
}

// Reweight by changing the bin normalisation and doing a reweight
void reweightwithnorm(TH1D* hmcpred, std::vector<double> MC, TF1* nom, TF1* pred, std::vector<double> binnorm) {
  hmcpred->Reset();
  for (int i = 0; i < MC.size(); i++) {
    int bin = hmcpred->FindBin(MC[i]);
    double fill = pred->Eval(MC[i])/nom->Eval(MC[i]);
    //if (fill > 100) fill = 100;
    hmcpred->AddBinContent(bin, binnorm[bin-1] * fill);
  }
}


// Calculate the normalisation probabilities using Newton Rhapson
double calcp(double data, double mc, double mcraw) {
  // Just follows 1993 Barlow Beeston paper
  // Observed data events
  double di = data;
  // Predicted scaled MC events
  double fi = mc;

  // Parameter of convenience
  double ti = 1-di/fi;

  double ai = mcraw;

  double pi = fi/(ai-fi-di);

  double Ai = fi/pi;

  // fi is the generated number of statistics in this bin
  // pi is the normalisation of the events
  // Ai is the true number of MC statistics in this bin
  // fi = pi*Ai
  // Want to find the pi (the normalisation parameters) and the Ai we don't really care about
  // First guess of pi
  //double pi = -1.0/ti;
  
  // Needs to be positive, puts constraint on pi
  //double Ai = ai/(1+pi*ti);
  // Solve eq 15 for the ti

  // sum over i ti*Ai = 0 (eq 16)

  // The likelihood we're trying to maximise
  double statllh = 0;
  if (di>0) statllh += fi-di+di*TMath::Log(di/fi);
  else statllh += fi;
  double bbllh = Ai-ai+ai*TMath::Log(ai/Ai);
  //double bbllh = ai*TMath::Log(Ai)-Ai;

  return statllh+bbllh;
}

//double calcp_lite(double dat, double mc, double mcraw) {
double calcp_lite(double dat, double mc, double uncert) {
  // fractional error
  if (uncert == 0) return 0;
  double fracerror = uncert;

  double temp = (mc*fracerror*fracerror-1);
  double temp2 = temp*temp + 4*dat*fracerror*fracerror;
  if (temp2 < 0) {
    std::cerr << "eurgh no root for b^2-4ac=" << temp2 << std::endl;
    throw;
  }

  double beta = (-1*temp + sqrt(temp2))/2.;

  double newmc = mc*beta;

  double stat = 0;
  if (dat == 0) stat = newmc;
  else if (newmc > 0) stat = newmc-dat+dat*TMath::Log(dat/newmc);
  //else return 9999999999999;
  double pen = (beta-1)*(beta-1)/(2*fracerror*fracerror);
  //double llh = -1*dat*TMath::Log(newmc)+newmc+(beta-1)*(beta-1)/(2*fracerror*fracerror);
  //std::cout << "dat, mc, beta, newmc, fracerror  = " << dat << " " << mc << " " << beta << " " << newmc << " " << fracerror << std::endl;
  //std::cout << "stat, pen = " << stat << " " << pen << " " << beta << std::endl;
  //std::cout << "stat/pen = " << stat/pen << std::endl;
  return stat+pen;
}

//BarlowBeeston() {
double LLHBarlowBeeston(TH1D *hData, TH1D *hMC, std::vector<double> binuncert) {
  double llht = 0.0;

  for (int i = 0; i < hData->GetXaxis()->GetNbins()+1; ++i) {
    // Calculate the pi for each i bin
    //double llh = calcp(hData->GetBinContent(i+1), hMC->GetBinContent(i+1), hMCraw->GetBinContent(i+1));
    double llh = calcp_lite(hData->GetBinContent(i+1), hMC->GetBinContent(i+1), binuncert[i]);
    //std::cout << i << " = " << llh << std::endl;
    // Add to the total likelihood
    llht += llh;
  }

  return llht;
}


void MCStats(int fittype, double normIn=200., double MCfactorIn = 10) {
  if (fittype > 7 || fittype < 0) {
    std::cerr << "I take 0, 1, 2, 3, 4, 5, 6 as arguments" << std::endl;
    std::cerr << "  0 = LLHStat_StatOnly" << std::endl;
    std::cerr << "  1 = LLHStat_StatOnlyPen" << std::endl;
    std::cerr << "  2 = LLHStat_PearsonNoPen" << std::endl;
    std::cerr << "  3 = LLHStat_PearsonPen" << std::endl;
    std::cerr << "  4 = LLHStat_StatWithFluc" << std::endl;
    std::cerr << "  5 = LLHStat_StatWithNorm" << std::endl;
    std::cerr << "  6 = LLH_BarlowBeestonGaus" << std::endl;
    std::cerr << "  7 = LLH_BarlowBeestonPoiss" << std::endl;
    exit(-1);
  }
  test_type testtype = static_cast<test_type>(fittype);
  // random number
  TRandom3* rnd = new TRandom3(10);
  // The mean of the parameter we want to measure
  double mean_prod = -2.5;
  // The sigma of the parameter we want to measure
  double sigma_prod = 3;

  // The number of data events
  double norm = normIn;

  // The factor of more MC that we have over data
  double MCfactor = MCfactorIn;

  std::vector<double> MC;

  // Vector of the binedges
  std::vector<double> binedges;
  int minx = -4;
  int maxx = 4;
  int currentindex=0;
  // At least 20 MC events per bin
  int minnum = 0;
  // Minimum width of each bin
  double mindist = 0.2;

  // Make some new MC following the Gaussian
  for (int i = 0; i < int(norm*MCfactor); i++) {
    double num = 0.7*rnd->Gaus(mean_prod, sigma_prod);
    if (num < minx || num > maxx) {
      i--;
      continue;
    }
    MC.push_back(num);
  }
  // Sort the MC
  sort(MC.begin(), MC.end()); 

  // Set some uniform binning
  for (double it = minx; it <= maxx; it += mindist) binedges.push_back(it);
  /*
  // Set the binning starting at -5 to +5
  // With at least 20 events per MC bin
  // And a bin width of 0.2
  for (int i = 0; i < MC.size()-minnum; i++) {
    // Check distance between points is larger than the minimum distance
    if (MC[i] - MC[currentindex] < mindist) continue;
    // 
    else if (i - currentindex < minnum) continue;
    else {
      binedges.push_back(MC[i]);
      currentindex=i;
    }
  }
  */
  binedges.push_back(maxx);

  //for (int i = 0; i < binedges.size(); ++i) std::cout << binedges[i] << std::endl;

  // The parameter values
  TH1D* hMC = new TH1D("hMC", ";p_{#mu} Fake", binedges.size()-1, &binedges[0]);

  // Fill the histogram with our MC
  for (int i=0; i<MC.size(); i++) {
    // Throw away rubbish
    if (MC[i] < binedges.front() || MC[i] > binedges.back()) continue;
    hMC->Fill(MC[i]);
  }
  hMC->Sumw2();
  hMC->Draw();

  // The bin normalisation, and uncerainty
  std::vector<double> binnorm, binnormp, binuncert;
  for(int i=1; i<=hMC->GetNbinsX(); i++) {
    binnorm.push_back(1);
    binnormp.push_back(1);
    if (hMC->GetBinContent(i) > 0) binuncert.push_back(hMC->GetBinError(i)/hMC->GetBinContent(i));
    else binuncert.push_back(0);
  }

  // The scaled MC histogram
  TH1D* hMCnorm = (TH1D*)hMC->Clone("hMCnorm");
  hMCnorm->Scale(1./MCfactor);

  // The data histogram
  TH1D* hData = (TH1D*)hMCnorm->Clone("hData");

  // Number of data points as a Gaussian around the N with error sqrt(N)
  //int ndatapts = rnd->Gaus(norm, sqrt(norm));
  int ndatapts = norm;
  hData->Reset();
  //for(int i=0; i<ndatapts; i++) hData->Fill(rnd->Gaus(mean, sigma));

  // Try making the data a scaled version of the MC
  for (int i = 0; i < hData->GetXaxis()->GetNbins()+1; ++i) {
    hData->SetBinContent(i+1, hMCnorm->GetBinContent(i+1));
    hData->SetBinError(i+1, sqrt(hData->GetBinContent(i+1)));
  }

  double mean = 0;
  double sigma = 1;

  TF1* gausnominal = new TF1("gausnom","gaus", MC.front(), MC.back());
  gausnominal->SetParameters(1, mean, sigma);
  TF1* gauspred = new TF1("gauspred","gaus", MC.front(), MC.back());
  gauspred->SetParameters(1, mean, sigma);

  // The proposed parameters
  double meanp = rnd->Gaus(mean, 0.5);
  double sigmap = rnd->Gaus(sigma, 0.5);
  //double sigmap = sigma;
  double normp = norm;

  gauspred->SetParameters(1, meanp, sigmap);
  // Do the first reweight
  if (testtype != LLHStat_StatWithNorm) reweight(hMCnorm, MC, gausnominal, gauspred);
  else reweightwithnorm(hMCnorm, MC, gausnominal, gauspred, binnormp);
  hMCnorm->Scale(norm/double(MC.size()));
  double LLH;

  // Calculate the likelihood
  if      (testtype == LLHStat_StatOnly) LLH = calcLLH(hData, hMCnorm, hMC, false);
  else if (testtype == LLHStat_StatOnlyPen) LLH = calcLLH(hData, hMCnorm, hMC, true);
  else if (testtype == LLHStat_PearsonNoPen) LLH = calcLLHwithPen(hData, hMCnorm, hMC, false);
  else if (testtype == LLHStat_PearsonPen) LLH = calcLLHwithPen(hData, hMCnorm, hMC, true);
  else if (testtype == LLHStat_StatWithFluc) LLH = calcLLHwithvar(hData, hMCnorm, hMC, rnd);
  else if (testtype == LLHStat_StatWithNorm) LLH = calcLLHwithnorm(hData, hMCnorm, binnormp, binuncert);
  else if (testtype == LLH_BarlowBeestonGauss) LLH = LLHBarlowBeeston(hData, hMCnorm, binuncert);

  // Proposed likelihood
  double LLHp = LLH;

  // Number of MCMC steps
  const int nsteps = 1E4;

  TString filename = "MCStats_fittype_";
  filename += testtype;
  filename += Form("_nData%2.2f_MCscale%1.1f_nsteps%i_muSigmaOnly_dataclone", normIn, MCfactorIn, nsteps);
  filename += ".root";

  double accProb;
  TFile* file = new TFile(filename,"RECREATE");
  TTree* tree = new TTree("mcstats","");
  tree->Branch("mean",&mean,"mean/D");
  tree->Branch("sigma",&sigma,"sigma/D");
  tree->Branch("norm",&norm,"norm/D");
  tree->Branch("llh",&LLH,"llh/D");
  tree->Branch("accprob",&accProb,"accProb/D");
  if (testtype == LLHStat_StatWithNorm) {
    for(int i=0; i<binnorm.size(); i++) {
      TString name="bin";
      name+=i;
      TString des = name;
      des+="/D";
      tree->Branch(name,&binnorm[i],des);
    }
  }

  //std::cout << meanp << " " << sigmap << " " << normp << std::endl;
  // Update the starting step
  mean = meanp;
  sigma = sigmap;

  // Run the MCMC
  int naccept = 0;
  for(int i=0; i<nsteps; i++) {

    if (i % (nsteps/20) == 0) {
      std::cout << "On step " << i << "/" << nsteps << " (" << int(double(i)/double(nsteps)*100.0) << "%)" << std::endl;
    }

    // The proposed parameters in the MCMC
    // Gaussian proposal around previously accepted step
    meanp = rnd->Gaus(mean, 0.02);
    sigmap = rnd->Gaus(sigma, 0.02);
    //normp = rnd->Gaus(norm, sqrt(norm)*0.5);
    //sigmap = sigma;
    normp = norm;

    // For method 2 throw the random bin
    if (testtype == LLHStat_StatWithNorm) {
      for (int j = 0; j < binnormp.size(); j++) {
        binnormp[j] = rnd->Gaus(binnorm[j], binuncert[j]*0.1);
      }
    }

    gauspred->SetParameters(1, meanp, sigmap);
    // reweight the monte carlo prediction
    if (testtype != LLHStat_StatWithNorm) reweight(hMCnorm, MC, gausnominal, gauspred);
    else reweightwithnorm(hMCnorm, MC, gausnominal, gauspred, binnormp);
    // Scale the MC back to data "POT"
    hMCnorm->Scale(norm/double(MC.size()));

    // Full Barlow-Beeston, which solves the normalisation in each bin
    if      (testtype == LLHStat_StatOnly) LLHp = calcLLH(hData, hMCnorm);
    else if (testtype == LLHStat_StatOnlyPen) LLHp = calcLLH(hData, hMCnorm, hMC, true);
    else if (testtype == LLHStat_PearsonNoPen) LLHp = calcLLHwithPen(hData, hMCnorm, hMC, false);
    else if (testtype == LLHStat_PearsonPen) LLHp = calcLLHwithPen(hData, hMCnorm, hMC, true);
    else if (testtype == LLHStat_StatWithFluc) LLHp = calcLLHwithvar(hData, hMCnorm, hMC, rnd);
    else if (testtype == LLHStat_StatWithNorm) LLHp = calcLLHwithnorm(hData, hMCnorm, binnormp, binuncert);
    else if (testtype == LLH_BarlowBeestonGauss) LLHp = LLHBarlowBeeston(hData, hMCnorm, binuncert);

    // Do the simple MCMC accept-reject
    accProb = TMath::Min(1., TMath::Exp(LLH-LLHp));
    double fRandom = rnd->Rndm();
    // Accept
    if (fRandom <= accProb) {
      mean = meanp;
      sigma = sigmap;
      norm = normp;
      LLH = LLHp;
      naccept++;
      // Update the bin normalisation parameter if accepted
      if (testtype == LLHStat_StatWithNorm) {
        for (int j = 0; j < binnorm.size(); j++) binnorm[j] = binnormp[j];
      }
    }
    tree->Fill();
  }
  std::cout << "fittype " << testtype << " had " << double(naccept)/double(nsteps)*100.0 << "% acceptance" << std::endl;
  tree->Write();
  hData->Write();
  hMC->Write();
  hMCnorm->Write();
  file->Close();
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Need three arguments: fittype, n_data, MC_scaling" << std::endl;
    return -1;
  }
  MCStats(std::atoi(argv[1]), std::atof(argv[2]), std::atof(argv[3]));
  return 0;
}
