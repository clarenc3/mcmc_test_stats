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

// Standard LLH calculator
double calcLLH(TH1D* hdata, TH1D* hmc) {
  double LLH=0;
  for(int i=1; i<=hdata->GetNbinsX(); i++) {
    double dat=hdata->GetBinContent(i);
    double mc=hmc->GetBinContent(i);
    if(dat>0) LLH+=mc-dat+dat*TMath::Log(dat/mc);
    else LLH+=mc;
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
    /*
    std::cout << "***START" << std::endl;
    std::cout << i << " " << LLH << std::endl;
    std::cout << dat << " " << mc << std::endl;
    std::cout << LLH << std::endl;
    std::cout << "***END" << std::endl;
    */
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

double calcLLHwithPen(TH1D* hdata, TH1D* hmc, bool withPen = true) {
  double LLH=0;
  for(int i=1; i<=hdata->GetNbinsX(); i++) {
    double dat=hdata->GetBinContent(i);
    double mc=hmc->GetBinContent(i);
    if (withPen && dat+mc>0) LLH += (dat-mc)*(dat-mc)/(dat+mc);
    else {
      if (dat != 0 && dat > 0) LLH += (dat-mc)*(dat-mc)/(dat);
      else if (mc > 0) LLH += (dat-mc)*(dat-mc)/(mc);
    }
  }
  return LLH;
}

//void reweight(TH1D* hmcpred, std::vector<double> MC, double mean, double sigma, double norm, TF1* nom, TF1* pred) {
void reweight(TH1D* hmcpred, std::vector<double> MC, TF1* nom, TF1* pred) {
  // Reset the Monte Carlo prediction
  hmcpred->Reset();
  //pred->SetParameters(1, mean, sigma);
  // Reweight each MC event
  for (int i = 0; i < MC.size(); i++) {
    double fill = pred->Eval(MC[i])/nom->Eval(MC[i]);
    //std::cout << nom->Eval(MC[i]) << std::endl;
    hmcpred->Fill(MC[i], fill);
  }
  // Scale back to data events
  //hmcpred->Scale(norm/double(MC.size()));
}

//void reweightwithnorm(TH1D* hmcpred, std::vector<double> MC, double mean, double sigma, double norm, TF1* nom, TF1* pred, std::vector<double> binnorm) {
void reweightwithnorm(TH1D* hmcpred, std::vector<double> MC, TF1* nom, TF1* pred, std::vector<double> binnorm) {
  hmcpred->Reset();
  //pred->SetParameters(1, mean, sigma);
  for (int i = 0; i < MC.size(); i++) {
    int bin = hmcpred->FindBin(MC[i]);
    hmcpred->AddBinContent(bin, binnorm[bin-1] * pred->Eval(MC[i])/nom->Eval(MC[i]));
  }
  //hmcpred->Scale(norm/double(MC.size()));
}


// Calculate the normalisation probabilities using Newton Rhapson
double calcp(double data, double mc, double MCscaling) {
  double minllh = 99999;


  // Just follows 1993 Barlow Beeston paper
  // Observed data events
  double di = data;
  // Predicted scaled MC events
  double fi = mc;

  // Parameter of convenience
  double ti = 1-di/fi;

  double ai = fi * MCscaling;

  // fi is the generated number of statistics in this bin
  // pi is the normalisation of the events
  // Ai is the true number of MC statistics in this bin
  // fi = pi*Ai
  // Want to find the pi (the normalisation parameters) and the Ai we don't really care about
  // First guess of pi
  double pi = -1.0/ti;
  
  // Needs to be positive, puts constraint on pi
  double Ai = ai/(1+pi*ti);
  // Solve eq 15 for the ti

  // sum over i ti*Ai = 0 (eq 16)

  // The likelihood we're trying to maximise
  //double llh = di*TMath::Log(fi)-fi + ai*TMath::Log(Ai)-Ai;

  return minllh;
}

//BarlowBeeston() {
double LLHBarlowBeeston(TH1D *hData, TH1D *hMC, double MCscaling) {
  double llht = 0.0;

  for (int i = 0; i < hData->GetXaxis()->GetNbins()+1; ++i) {
    // Calculate the pi for each i bin
    double llh = calcp(hData->GetBinContent(i+1), hMC->GetBinContent(i+1), MCscaling);
    // Add to the total likelihood
    llht += llh;
  }

  return llht;
}


void MCStats(int fittype, double normIn=200., double MCfactorIn = 10) {
  if (fittype > 4 || fittype < 0) {
    std::cerr << "I take 0, 1, 2, 3, 4 as arguments" << std::endl;
    exit(-1);
  }
  // random number
  TRandom3* rnd = new TRandom3(10);
  // The mean of the parameter we want to measure
  double mean = 0;
  // The sigma of the parameter we want to measure
  double sigma = 1;

  // The number of data events
  double norm = normIn;

  // The factor of more MC that we have over data
  double MCfactor = MCfactorIn;

  std::vector<double> MC;

  // Make some new MC following the Gaussian
  for (int i = 0; i < int(norm*MCfactor); i++) MC.push_back(rnd->Gaus(mean, sigma));
  // Sort the MC
  sort(MC.begin(), MC.end()); 

  // Vector of the binedges
  std::vector<double> binedges;
  int minx = -4;
  int maxx = 4;
  int currentindex=0;
  // At least 20 MC events per bin
  int minnum = 0;
  // Minimum width of each bin
  double mindist = 0.2;

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
  TH1D* hMC = new TH1D("hMC", ";Some Parameter", binedges.size()-1, &binedges[0]);

  // Fill the histogram with our MC
  for(int i=0; i<MC.size(); i++) hMC->Fill(MC[i]);
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

  TF1* gausnominal = new TF1("gausnom","gaus",minx,maxx);
  gausnominal->SetParameters(1, mean, sigma);
  TF1* gauspred = new TF1("gauspred","gaus",minx,maxx);
  gauspred->SetParameters(1, mean, sigma);

  // The proposed parameters
  double meanp = rnd->Gaus(mean, 0.5);
  double sigmap = rnd->Gaus(sigma, 0.5);
  //double sigmap = sigma;
  double normp = norm;

  gauspred->SetParameters(1, meanp, sigmap);
  // Do the first reweight
  if (fittype < 4) reweight(hMCnorm, MC, gausnominal, gauspred);
  else if (fittype == 4) reweightwithnorm(hMCnorm, MC, gausnominal, gauspred, binnormp);
  hMCnorm->Scale(norm/double(MC.size()));
  double LLH;

  // Calculate the likelihood
  if      (fittype == 0) LLH = calcLLH(hData, hMCnorm);
  else if (fittype == 1) LLH = calcLLHwithPen(hData, hMCnorm, true);
  else if (fittype == 2) LLH = calcLLHwithPen(hData, hMCnorm, false);
  else if (fittype == 3) LLH = calcLLHwithvar(hData, hMCnorm, hMC, rnd);
  else if (fittype == 4) LLH = calcLLHwithnorm(hData, hMCnorm, binnormp, binuncert);

  // Proposed likelihood
  double LLHp = LLH;

  // Number of MCMC steps
  const int nsteps = 1E4;

  TString filename = "MCStats_fittype_";
  filename += fittype;
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
  if (fittype == 4) {
    for(int i=0; i<binnorm.size(); i++) {
      TString name="bin";
      name+=i;
      TString des = name;
      des+="/D";
      tree->Branch(name,&binnorm[i],des);
    }
  }

  //std::cout << meanp << " " << sigmap << " " << normp << std::endl;

  // Run the MCMC
  int naccept = 0;
  for(int i=0; i<nsteps; i++) {

    //std::cout << meanp << " " << sigmap << " " << normp << std::endl;

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
    if (fittype == 4) {
      for (int j = 0; j < binnormp.size(); j++) {
        binnormp[j] = rnd->Gaus(binnorm[j], binuncert[j]*0.1);
      }
    }

    gauspred->SetParameters(1, meanp, sigmap);
    // reweight the monte carlo prediction
    if (fittype < 4) reweight(hMCnorm, MC, gausnominal, gauspred);
    else if (fittype == 4) reweightwithnorm(hMCnorm, MC, gausnominal, gauspred, binnormp);
    hMCnorm->Scale(norm/double(MC.size()));
    // Scale the MC back to data "POT"

    // Calculate the likelihood for the proposed step
    if      (fittype == 0) LLHp = calcLLH(hData, hMCnorm);
    else if (fittype == 1) LLHp = calcLLHwithPen(hData, hMCnorm, true);
    else if (fittype == 2) LLHp = calcLLHwithPen(hData, hMCnorm, false);
    else if (fittype == 3) LLHp = calcLLHwithvar(hData, hMCnorm, hMC, rnd);
    else if (fittype == 4) LLHp = calcLLHwithnorm(hData, hMCnorm, binnormp, binuncert);
    // Full Barlow-Beeston, which solves the normalisation in each bin

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
      if (fittype == 4) {
        for (int j = 0; j < binnorm.size(); j++) binnorm[j] = binnormp[j];
      }
    }
    tree->Fill();
  }
  std::cout << "fittype " << fittype << " had " << double(naccept)/double(nsteps)*100.0 << "% acceptance" << std::endl;
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
