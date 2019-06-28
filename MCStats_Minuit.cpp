#include "MCStats_Minuit.h"

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
    hmcpred->Fill(MC[i], fill);
  }
}

// Reweight by changing the bin normalisation and doing a reweight
void reweightwithnorm(TH1D* hmcpred, std::vector<double> MC, TF1* nom, TF1* pred, std::vector<double> binnorm) {
  hmcpred->Reset();
  for (int i = 0; i < MC.size(); i++) {
    int bin = hmcpred->FindBin(MC[i]);
    double fill = pred->Eval(MC[i])/nom->Eval(MC[i]);
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

storer::storer() {
  rnd = new TRandom3(1234);
}

// Do the evaluation
double storer::DoEval(const double *x) {
  // Update internal parameters
  meanp = x[0];
  sigmap = x[1];

  // Do the reweight
  ReWeight();

  return GetLikelihood();
}

void storer::Setup() {

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

  // The parameter values
  hMC = new TH1D("hMC", ";p_{#mu} Fake", binedges.size()-1, &binedges[0]);

  // Fill the histogram with our MC
  for (int i=0; i<MC.size(); i++) {
    // Throw away rubbish
    if (MC[i] < binedges.front() || MC[i] > binedges.back()) continue;
    hMC->Fill(MC[i]);
  }

  // Now set the errors
  for (int i = 0; i < hMC->GetNbinsX()+1; ++i) {
    hMC->SetBinError(i+1, sqrt(hMC->GetBinContent(i+1)));
  }

  // The bin normalisation, and uncerainty
  for(int i=1; i<=hMC->GetNbinsX(); i++) {
    binnorm.push_back(1);
    if (hMC->GetBinContent(i) > 0) binuncert.push_back(hMC->GetBinError(i)/hMC->GetBinContent(i));
    else binuncert.push_back(0);
  }

  // The scaled MC histogram
  hMCnorm = (TH1D*)hMC->Clone("hMCnorm");
  hMCnorm->Scale(1./MCfactor);

  // The data histogram
  hData = (TH1D*)hMCnorm->Clone("hData");

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

  // Set up the weighting
  SetWeighting();

  // Update the data
  for (int i = 0; i < hData->GetXaxis()->GetNbins()+1; ++i) {
    double center = hData->GetBinCenter(i+1);
    //double weight = gauspred->Eval(center)/gausnom->Eval(center);
    double weight = gauspred->Eval(center);
    hData->SetBinContent(i+1, hData->GetBinContent(i+1)*(1+weight));
  }

}

void storer::SetWeighting() {
  gausnom = new TF1("gausnom","gaus", MC.front(), MC.back());
  gausnom->SetParameters(1, mean, sigma);

  gauspred = new TF1("gauspred","gaus", MC.front(), MC.back());
  gauspred->SetParameters(1, mean, sigma);
}

void storer::ReWeight() {

  // Set the parameters of the Gaussian update
  gauspred->SetParameters(1, meanp, sigmap);

  // Reweight
  if (testtype != LLHStat_StatWithNorm) reweight(hMCnorm, MC, gausnom, gauspred);
  else reweightwithnorm(hMCnorm, MC, gausnom, gauspred, binnorm);

  // Scale back
  hMCnorm->Scale(norm/double(MC.size()));
}


double storer::GetLikelihood() {
  // Calculate the likelihood
  double LLH = 0;
  if      (testtype == LLHStat_StatOnly) LLH = calcLLH(hData, hMCnorm, hMC, false);
  else if (testtype == LLHStat_StatOnlyPen) LLH = calcLLH(hData, hMCnorm, hMC, true);
  else if (testtype == LLHStat_PearsonNoPen) LLH = calcLLHwithPen(hData, hMCnorm, hMC, false);
  else if (testtype == LLHStat_PearsonPen) LLH = calcLLHwithPen(hData, hMCnorm, hMC, true);
  else if (testtype == LLHStat_StatWithFluc) LLH = calcLLHwithvar(hData, hMCnorm, hMC, rnd);
  else if (testtype == LLHStat_StatWithNorm) LLH = calcLLHwithnorm(hData, hMCnorm, binnorm, binuncert);
  else if (testtype == LLH_BarlowBeestonGauss) LLH = LLHBarlowBeeston(hData, hMCnorm, binuncert);
  return LLH;
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

  storer Storage;
  // Set up the defaults (the production values)
  Storage.SetMeanProd(mean_prod);
  Storage.SetSigmaProd(sigma_prod);

  // The fitting values we're trying to find
  Storage.SetMean(2);
  Storage.SetSigma(0.9);
  Storage.SetNorm(norm);
  Storage.SetMCfact(MCfactor);

  // Now set up the MC, the scaled MC and the data
  Storage.Setup();

  TString filename = "MCStats_minuit2_fittype_";
  filename += testtype;
  filename += Form("_nData%2.2f_MCscale%1.1f_muSigmaOnly_dataclone", normIn, MCfactorIn);
  filename += ".root";
  TFile* file = new TFile(filename,"RECREATE");

  std::string fit = "Minuit2";
  std::string type = "Migrad";
  ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer(fit, type);

  Storage.GetData()->Write("data_pre");
  Storage.GetMC()->Write("mc_pre");
  Storage.GetUnscaledMC()->Write("mcunscaled_pre");

  const int nparams = 2;
  ROOT::Math::Functor f(&Storage, &storer::DoEval, 2);

  min->SetFunction(f);
  double var[nparams] = {0, 1};
  double step[nparams] = {0.01, 0.01};
  min->SetVariable(0, "x", var[0], step[0]);
  min->SetVariable(1, "y", var[1], step[1]);
  min->Minimize();

  const double *xs = min->X();
  std::cout << "Minimum: f(" << xs[0] << ", " << xs[1] << "): " << min->MinValue() << std::endl;

  file->cd();
  Storage.GetData()->Write();
  Storage.GetMC()->Write();
  Storage.GetUnscaledMC()->Write();
  file->Close();

  std::cout << "Wrote to " << file->GetName() << std::endl;
}

int main(int argc, char **argv) {
  if (argc != 4) {
    std::cerr << "Need three arguments: fittype, n_data, MC_scaling" << std::endl;
    return -1;
  }
  MCStats(std::atoi(argv[1]), std::atof(argv[2]), std::atof(argv[3]));
  return 0;
}
