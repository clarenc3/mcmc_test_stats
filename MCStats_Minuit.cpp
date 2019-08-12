#include "MCStats_Minuit.h"

// Standard LLH calculator
// Add in MC statistical penalty
double calcLLH(TH1D* hdata, TH1D* hmc) {
  double LLH=0;
  for(int i=1; i<=hdata->GetNbinsX(); i++) {
    double dat=hdata->GetBinContent(i);
    double mc=hmc->GetBinContent(i);
    if(dat>0) LLH+=2*(mc-dat+dat*TMath::Log(dat/mc));
    else LLH+=2*mc;
  }
  return LLH;
}

// Simple LLH calculator with a random variation on the bin content
// herrors is the nominal MC histogram
double calcLLHwithvar(TH1D* hdata, TH1D* hmc, TH1D* herrors, TRandom3* rnd) {
  double LLH=0;
  for(int i = 0; i < hdata->GetNbinsX()+1; i++) {
    double dat = hdata->GetBinContent(i+1);
    // Throw the MC within the fractional MC statistical error
    double origmc = hmc->GetBinContent(i+1);
    // Get the fractional uncertainty: standard deviation 
    double uncert = sqrt(hmc->GetSumw2()->GetArray()[i+1]);
    double errors = herrors->GetBinContent(i+1);
    //std::cout << uncert << " " << errors << std::endl;
    double mc;

    if (uncert > 0) mc = rnd->Gaus(origmc, uncert);
    else mc = hmc->GetBinContent(i+1);

    if (dat > 0 && mc > 0) LLH += 2*(mc-dat+dat*TMath::Log(dat/mc));
    else LLH += 2*mc;
  }
  return LLH;
}


double calcLLHwithnorm(TH1D* hdata, TH1D* hmc, std::vector<double> norm, std::vector<double> uncert) {
  double LLH=0;
  for(int i = 0; i < hdata->GetNbinsX()+1; i++) {
    double dat = hdata->GetBinContent(i+1);
    double mc = hmc->GetBinContent(i+1);
    if(dat>0) LLH += 2*(mc-dat+dat*TMath::Log(dat/mc));
    else LLH += 2*mc;
    if (uncert[i] > 0) LLH += pow((norm[i]-1)/uncert[i],2);
  }
  return LLH;
}

// Likelihood a la Pearson, which assumes a Gaussian distributed statistic with an error of sqrt(N) --- i.e. its error is Poissonian
// withPen = true penalises on both data and MC statistics
// withPen = false penalises on data only
double calcLLHwithPen(TH1D* hdata, TH1D* hmc, TH1D *herrors, bool withPen = true) {
  double LLH=0;
  for(int i = 0; i < hdata->GetNbinsX()+1; i++) {
    double dat = hdata->GetBinContent(i+1);
    double mc = hmc->GetBinContent(i+1);
    //double mcraw = herrors->GetBinContent(i+1);
    double sumw2 = hmc->GetSumw2()->GetArray()[i+1];
    if (withPen && dat+mc>0) LLH += (dat-mc)*(dat-mc)/(dat+sumw2);
    else if (!withPen) {
      if (dat != 0 && dat > 0) LLH += (dat-mc)*(dat-mc)/(dat);
      else if (mc > 0) LLH += (dat-mc)*(dat-mc)/(sumw2);
    }
  }
  return LLH;
}

// Try calculating LLH with Tianlu's method
double calcLLHTianlu(TH1D *hdata, TH1D *hmc) {
  double llh = 0.0; 
  for (int i = 0; i < hdata->GetXaxis()->GetNbins()+1; ++i) {
    double data = hdata->GetBinContent(i+1);
    double mc = hmc->GetBinContent(i+1);
    double w2 = hmc->GetSumw2()->GetArray()[i+1];
    if (w2 == 0) continue;

    double b = mc/w2;
    double a = mc*b+1;
    double k = data;

    llh += a * std::log(b) + std::lgamma(k+a)-std::lgamma(k+1)-((k+a)*std::log1p(b))-std::lgamma(a);
  }

  return -2*llh;
}

void reweight(TH1D* hmcpred, std::vector<double> MC, TF1* sig, TF1 *bkg, double scale) {
  // Reset the Monte Carlo prediction
  hmcpred->Reset();
  double nominal = scale/MC.size();
  // Reweight each MC event
  for (int i = 0; i < MC.size(); i++) {
    double pred = sig->Eval(MC[i]);
    double nom = bkg->Eval(MC[i]);
    //double weight = (1+pred/nom)*nominal;
    double weight = (1+pred/nom)*nominal;
    hmcpred->Fill(MC[i], weight);
  }
}

// Reweight by changing the bin normalisation and doing a reweight
void reweightwithnorm(TH1D* hmcpred, std::vector<double> MC, TF1* nom, TF1* pred, std::vector<double> binnorm) {
  hmcpred->Reset();
  for (int i = 0; i < MC.size(); i++) {
    int bin = hmcpred->FindBin(MC[i]);
    double fill = 1+pred->Eval(MC[i])/nom->Eval(MC[i]);
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
  if (dat == 0) stat = 2*newmc;
  else if (newmc > 0) stat = 2*(newmc-dat+dat*TMath::Log(dat/newmc));
  double pen = (beta-1)*(beta-1)/(fracerror*fracerror);
  return stat+pen;
}

double LLHBarlowBeeston(TH1D *hData, TH1D *hMC, std::vector<double> binuncert) {
  double llht = 0.0;

  for (int i = 0; i < hData->GetXaxis()->GetNbins()+1; ++i) {
    // Calculate the pi for each i bin
    //double llh = calcp(hData->GetBinContent(i+1), hMC->GetBinContent(i+1), hMCraw->GetBinContent(i+1));
    double data = hData->GetBinContent(i+1);
    double mc = hMC->GetBinContent(i+1);
    double mcuncert = sqrt(hMC->GetSumw2()->GetArray()[i+1])/mc;
    //std::cout << mc << " " << mcuncert << std::endl;
    //std::cout << mc << " " << mcuncert << std::endl;
    //std::cout << binuncert[i] << std::endl;
    //std::cout << sqrt(mcuncert)/mc << std::endl;
    double llh = calcp_lite(data, mc, mcuncert);
    //std::cout << i << " = " << llh << std::endl;
    // Add to the total likelihood
    llht += llh;
  }
  //throw;

  return llht;
}

storer::storer(int seed=0) {
  rnd = new TRandom3(seed);
  mean = -999;
  sigma = -999;
  norm = -999;
  MCfactor = -999;
  meanp = -999;
  sigmap = -999;
  RefLLH = -999;
}

// Destructor
storer::~storer() {
  delete hData;
  delete hMC;
  delete hMCnorm;
  delete sig;
  delete bkg;
  delete rnd;
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

  double minx = -4;
  double maxx = 4;

  // Make a reference TF1 to draw from
  bkg = new TF1("bkg", "0.5*exp(-x/2)", minx, maxx);

  // Make some new MC following the Gaussian
  for (int i = 0; i < int(norm*MCfactor); i++) {
    double num = bkg->GetRandom();
    if (num < minx || num > maxx) {
      i--;
      continue;
    }
    MC.push_back(num);
  }
  // Sort the MC
  sort(MC.begin(), MC.end()); 


  // Vector of the binedges
  std::vector<double> binedges;

  // At least 10 MC events per bin
  int minimum = 10;
  // Minimum width of each bin
  double mindist = 0.1;
  int currentindex = 0;
  for (int i = 0; i < MC.size(); ++i) {
    // Check the number of events is greater than threshold
    if (i - currentindex < minimum) continue;
    if (MC[i] - MC[currentindex] < mindist) continue;
    binedges.push_back(MC[i]);
    currentindex = i;
  }

  // The parameter values
  hMC = new TH1D("hMC", ";p_{#mu} Fake", binedges.size()-1, &binedges[0]);
  hMC->Sumw2();

  // Fill the histogram with our MC
  for (int i=0; i<MC.size(); i++) {
    // Throw away rubbish
    if (MC[i] < binedges.front() || MC[i] > binedges.back()) continue;
    hMC->Fill(MC[i]);
  }

  // Now set the errors
  //for (int i = 0; i < hMC->GetNbinsX()+1; ++i) {
  //hMC->SetBinError(i+1, sqrt(hMC->GetBinContent(i+1)));
  //}

  // The bin normalisation, and uncerainty
  for (int i = 0; i < hMC->GetNbinsX()+1; i++) {
    binnorm.push_back(1);
    if (hMC->GetBinContent(i+1) > 0) {
      binuncert.push_back(hMC->GetBinError(i+1)/hMC->GetBinContent(i+1));
    } else {
      binuncert.push_back(0);
    }
  }

  // The scaled MC histogram
  hMCnorm = (TH1D*)hMC->Clone("hMCnorm");
  hMCnorm->Sumw2();
  // Scale down by the MC factor
  hMCnorm->Scale(1./MCfactor);

  // The data histogram
  hData = (TH1D*)hMCnorm->Clone("hData");
  hData->Sumw2();

  // Number of data points as a Gaussian around the N with error sqrt(N)
  //int ndatapts = rnd->Gaus(norm, sqrt(norm));
  int ndatapts = norm;
  hData->Reset();
  // Make the TF1 that the data is sampling from
  sig = new TF1("sig", "-1/sqrt(2*3.1415*[1])*exp(-pow(x-[0], 2)/(2*[1]))", minx, maxx);
  sig->SetParameter(0, mean);
  sig->SetParameter(1, sigma);

  for (int i = 0; i < ndatapts; i++) {
    // First get the random number from the TF1
    double ran = bkg->GetRandom();    
    hData->Fill(ran);
  }

  //TFile *temp = new TFile("temp.root", "recreate");
  //temp->cd();
  //bkg->Write("bkg");
  //sig->Write("sig");

  //hData->Write("pre_weight");

  // Update the data
  for (int i = 0; i < hData->GetXaxis()->GetNbins()+1; ++i) {
    double center = hData->GetBinCenter(i+1);
    double numer = sig->Eval(center);
    double denom = bkg->Eval(center);
    double weight = numer/denom;
    //double weight = numer;
    //std::cout << center << " = " << " numer: " << numer << " denom: " << denom << " ratio: " << numer/denom << std::endl;
    //hData->SetBinContent(i+1, hData->GetBinContent(i+1)*(1+weight));
    hData->SetBinContent(i+1, hData->GetBinContent(i+1)*(1+weight));
  }
  //temp->cd();
  //hData->Write("post_weight");
  //throw;

  // Calculate the reference distribution
  meanp = mean;
  sigmap = sigma;
  //hMCnorm->Write("hMCnorm_pre");
  ReWeight();
  RefLLH = GetLikelihood(); 
  // Reset the variables
  meanp = mean;
  sigmap = sigma;

  //hMC->Write("hMC");
  //hMCnorm->Write("hMCnorm_post");
  //temp->Close();
  //throw;

  // Check the values
  CheckMe();
}

bool storer::CheckMe() {

  bool good = true;
  if (mean == -999) {
    std::cerr << "mean not set" << std::endl;
    good = false;
  }
  if (sigma == -999) {
    std::cerr << "sigma not set" << std::endl;
    good = false;
  }
  if (norm == -999) {
    std::cerr << "norm not set" << std::endl;
    good = false;
  }
  if (MCfactor == -999) {
    std::cerr << "MCfactor not set" << std::endl;
    good = false;
  }
  if (meanp == -999) {
    std::cerr << "meanp not set" << std::endl;
    good = false;
  }
  if (sigmap == -999) {
    std::cerr << "sigmap not set" << std::endl;
    good = false;
  }
  if (RefLLH == -999) {
    std::cerr << "RefLLH not set" << std::endl;
    good = false;
  }

  if (!good) {
    std::cerr << "Found error, will quit" << std::endl;
    throw;
  }

  return good;
}

void storer::ReWeight() {
  // Set the parameters of the Gaussian update
  sig->SetParameters(meanp, sigmap);

  // Reweight
  reweight(hMCnorm, MC, sig, bkg, norm);
}


double storer::GetLikelihood() {
  // Calculate the likelihood
  double LLH = 0;
  if      (testtype == LLHStat_StatOnly) LLH = calcLLH(hData, hMCnorm);
  else if (testtype == LLHStat_PearsonNoPen) LLH = calcLLHwithPen(hData, hMCnorm, hMC, false);
  else if (testtype == LLHStat_PearsonPen) LLH = calcLLHwithPen(hData, hMCnorm, hMC, true);
  //else if (testtype == LLHStat_StatWithFluc) LLH = calcLLHwithvar(hData, hMCnorm, hMC, rnd);
  //else if (testtype == LLHStat_StatWithNorm) LLH = calcLLHwithnorm(hData, hMCnorm, binnorm, binuncert);
  else if (testtype == LLH_BarlowBeestonGauss) LLH = LLHBarlowBeeston(hData, hMCnorm, binuncert);
  else if (testtype == LLH_Tianlu) LLH = calcLLHTianlu(hData, hMCnorm);

  //std::cout << "Poisson: " << calcLLH(hData, hMCnorm, hMC, false) << std::endl;
  //std::cout << "Poisson penalty: " << calcLLH(hData, hMCnorm, hMC, true) << std::endl;
  //std::cout << "Pearson no penalty: " << calcLLHwithPen(hData, hMCnorm, hMC, false) << std::endl;
  //std::cout << "Pearson penalty: " << calcLLHwithPen(hData, hMCnorm, hMC, true) << std::endl;
  //std::cout << "Poisson with var: " << calcLLHwithvar(hData, hMCnorm, hMC, rnd) << std::endl;
  //std::cout << "Poisson with norm: " << calcLLHwithnorm(hData, hMCnorm, binnorm, binuncert) << std::endl;
  //std::cout << "BB-lite: " << LLHBarlowBeeston(hData, hMCnorm, binuncert) << std::endl;
  //std::cout << "Tianlu: " << calcLLHTianlu(hData, hMCnorm) << std::endl;
  //std::cout << "*" << std::endl;
  //
  hMCnorm->SetTitle(Form("#chi^{2}/nbins=%2.2f/%i", LLH, hMCnorm->GetXaxis()->GetNbins()));

  return LLH;
}



void MCStats(int fittype, double normIn=200., double MCfactorIn = 10) {
  if (fittype > nLLH || fittype < 0) {
    std::cerr << "I take ";
    for (int i = 0; i < nLLH; ++i) std::cout << i << " ";
    std::cout << "as arguments" << std::endl;
    std::cerr << "  0 = LLHStat_StatOnly" << std::endl;
    std::cerr << "  1 = LLHStat_PearsonNoPen" << std::endl;
    std::cerr << "  2 = LLHStat_PearsonPen" << std::endl;
    //std::cerr << "  3 = LLHStat_StatWithFluc" << std::endl;
    //std::cerr << "  4 = LLHStat_StatWithNorm" << std::endl;
    std::cerr << "  3 = LLH_BarlowBeestonGaus" << std::endl;
    std::cerr << "  4 = LLH_Tianlu" << std::endl;
    exit(-1);
  }

  // The number of data events
  double norm = normIn;

  // The factor of more MC that we have over data
  double MCfactor = MCfactorIn;

  TString filename = "MCStats_minuit2_fittype_";
  filename += fittype;
  filename += Form("_nData%2.2f_MCscale%1.1f_muSigmaOnly_newshape_lots", normIn, MCfactorIn);
  filename += ".root";
  TFile* file = new TFile(filename, "RECREATE");

  TTree *tree = new TTree("output", "");
  double mean;
  double sigma;
  double mean_error;
  double sigma_error;
  double edm2;
  double min2;
  double min2_true;
  int status;
  tree->Branch("mean", &mean, "mean/D");
  tree->Branch("sigma", &sigma, "sigma/D");
  tree->Branch("error_mean", &mean_error, "error_mean/D");
  tree->Branch("error_sigma", &sigma_error, "error_sigma/D");
  tree->Branch("minimum_llh", &min2, "minimum_llh/D");
  tree->Branch("minimum_llh_true", &min2_true, "minimum_llh_true/D");
  tree->Branch("edm", &edm2, "edm/D");
  tree->Branch("status", &status, "status/I");

  const int ntest = 500;
  int nBadFits = 0;
  // The parameters that we're fitting
  double mean_value = -4.2;
  double sigma_value = 0.2;
  for (int i = 0; i < ntest; ++i) {

    storer Storage(i);
    // Set up the defaults (the production values)
    Storage.SetTestType(static_cast<test_type>(fittype));

    // The fitting values we're trying to find
    Storage.SetMean(mean_value);
    Storage.SetSigma(sigma_value);
    // How many data events
    Storage.SetNorm(norm);
    // Factor how many more MC than data we have
    Storage.SetMCfact(MCfactor);

    // Now set up the MC, the scaled MC and the data
    Storage.Setup();

    // Get the true reference LLH
    min2_true = Storage.GetRefLLH();

    std::string fit = "Minuit2";
    std::string type = "Migrad";
    ROOT::Math::Minimizer *min = ROOT::Math::Factory::CreateMinimizer(fit, type);

    file->cd();
    Storage.GetData()->Write(Form("data_pre_%i", i));
    Storage.GetMC()->Write(Form("mc_pre_%i", i));
    Storage.GetUnscaledMC()->Write(Form("mcunscaled_pre_%i", i));

    const int nparams = 2;
    ROOT::Math::Functor f(&Storage, &storer::DoEval, 2);

    min->SetFunction(f);
    double var[nparams] = {0, 1};
    double step[nparams] = {0.01, 0.01};
    min->SetVariable(0, "x", var[0], step[0]);
    min->SetVariable(1, "y", var[1], step[1]);
    min->Minimize();

    //min->PrintResults();
    const double *postfit = min->X();
    const double *errors = min->Errors();
    const double minimum = min->MinValue();
    const double edm = min->Edm();
    std::cout << "Minimum: f(" << postfit[0] << "+/-" << errors[0] << ", " << postfit[1] << "+/-" << errors[1] <<  "): " << minimum << ", edm = " << edm << std::endl;

    mean = postfit[0];
    sigma = postfit[1];
    mean_error = errors[0];
    sigma_error = errors[1];
    edm2 = edm;
    min2 = minimum;
    status = min->Status();
    if (status != 0) nBadFits++;

    // Run Minos errors
    double errup, errdown;
    int runopt;
    min->GetMinosError(0, errup, errdown, runopt);
    mean_error = std::max(errup, errdown);
    min->GetMinosError(1, errup, errdown, runopt);
    sigma_error = std::max(errup, errdown);

    /*
       bool ok = min->GetMinosError(0, errup, errdown, runopt);
       std::cout << ok << " for getting Minos Error for param 1" << std::endl;
       std::cout << "mean_error=" << mean_error << std::endl;
       std::cout << "minos + = " << errup << std::endl;
       std::cout << "minos - = " << errdown << std::endl;
       std::cout << "runopt = " << runopt << std::endl;
       */

    tree->Fill();

    Storage.GetData()->Write(Form("data_post_%i", i));
    Storage.GetMC()->Write(Form("mc_post_%i", i));
    Storage.GetUnscaledMC()->Write(Form("mcunscaled_post_%i", i));

  }
  tree->Write();
  file->Close();

  std::cout << "Found " << nBadFits << "/" << ntest << " failed fits (" << double(nBadFits)/ntest*100.0 << "%)" << std::endl;
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
