const int nFiles = 5;
TFile *file[nFiles];
TTree *tree[nFiles];
TH1D *plots[nFiles];
TCanvas *canv;
std::string output;

void SetStyle(TH1D *plot, int i) {
  plot->SetLineWidth(2);
  if (i == 0) {
    plot->SetLineColor(kBlack);
    plot->SetTitle("Poisson");
  } else if (i == 1) {
    plot->SetTitle("Pearson");
    plot->SetLineColor(kBlue-7);
  } else if (i == 2) {
    plot->SetTitle("Pearson+MC stat");
    plot->SetLineColor(kGreen-2);
  } else if (i == 3) {
    plot->SetTitle("BB-lite");
    plot->SetLineColor(kYellow-2);
  } else if (i == 4) {
    plot->SetTitle("IceCube");
    plot->SetLineColor(kRed-3);
  }
  plot->GetXaxis()->SetLabelSize(plot->GetXaxis()->GetLabelSize()*2.5);
  plot->GetYaxis()->SetLabelSize(plot->GetYaxis()->GetLabelSize()*2.5);
  plot->GetYaxis()->SetNdivisions(205);
}

void CalcCoverage() {

  TPaveText *text[nFiles];
  for (int i = 0; i < nFiles; ++i) {

    text [i] = new TPaveText(0., 0., 1., 1., "NDC");

    int nentries = tree[i]->GetEntries();

    double mean;
    double error_mean;
    double sigma;
    double error_sigma;
    double status;
    double edm;
    tree[i]->SetBranchAddress("mean", &mean);
    tree[i]->SetBranchAddress("error_mean", &error_mean);
    tree[i]->SetBranchAddress("sigma", &sigma);
    tree[i]->SetBranchAddress("error_sigma", &error_sigma);
    tree[i]->SetBranchAddress("status", &status);
    tree[i]->SetBranchAddress("edm", &edm);

    double true_mean = 1.0;
    double true_sigma = 0.1;
    int ngood_mean = 0;
    int ngood_sigma = 0;
    int nbad_fits = 0;
    for (int j = 0; j < nentries; ++j) {
      tree[i]->GetEntry(j);
      if (status != 0) {
        nbad_fits++;
        //continue;
      }
      if (mean + error_mean > true_mean && mean - error_mean < true_mean) ngood_mean++;
      if (sigma + error_sigma > true_sigma && sigma - error_sigma < true_sigma) ngood_sigma++;

    }

    std::string title;
    int color;
    if (i == 0) {
      color = kBlack;
      title = "Poisson";
    } else if (i == 1) {
      title = "Pearson";
      color = kBlue-7;
    } else if (i == 2) {
      title = "Pearson+MC stat";
      color = kGreen-2;
    } else if (i == 3) {
      title = "BB-lite";
      color  = kYellow-2;
    } else if (i == 4) {
      title = "IceCube";
      color = kRed-3;
    }

    text[i]->AddText(title.c_str());
    ((TText*)text[i]->GetListOfLines()->Last())->SetTextColor(color);
    text[i]->AddText(Form("#mu=%i/%i (%2.1f%)", ngood_mean, nentries, double(ngood_mean)/nentries*100.));
    ((TText*)text[i]->GetListOfLines()->Last())->SetTextColor(color);
    text[i]->AddText(Form("#sigma=%i/%i (%2.1f%)", ngood_sigma, nentries, double(ngood_sigma)/nentries*100.));
    ((TText*)text[i]->GetListOfLines()->Last())->SetTextColor(color);
    text[i]->AddText(Form("n_{bad}=%i/%i (%2.1f%)", nbad_fits, nentries, double(nbad_fits)/nentries*100.));
    ((TText*)text[i]->GetListOfLines()->Last())->SetTextColor(color);
  }

  canv->cd();
  canv->Clear();
  canv->Divide(2,3);
  for (int i = 0; i < nFiles; ++i) {
    canv->cd(i+1);
    text[i]->Draw();
  }
  canv->Print(Form("%s.pdf", output.c_str()));
  canv->Clear();
}

TH1D* FillMe(TH1D *plot) {
  TH1D *plotcp = (TH1D*)plot->Clone(Form("%s_clone", plot->GetName()));
  plotcp->Reset();
  plotcp->SetFillColor(plot->GetLineColor());
  plotcp->SetLineWidth(0);
  plotcp->SetFillStyle(1001);
  // Find maxmimum bin
  int maxbin = plot->GetMaximumBin();
  int binup = 0;
  int bindown = 0;
  double integral = 0;
  if (plot->Integral() == 0) return plotcp;

  while (integral < 0.68 && maxbin-bindown > 0 && maxbin+binup < plot->GetXaxis()->GetNbins()) {
    if (plot->GetBinContent(maxbin+binup) > plot->GetBinContent(maxbin-bindown)) {
      integral += (plot->GetBinContent(maxbin+binup))/plot->Integral();
      plotcp->SetBinContent(maxbin+binup, plot->GetBinContent(maxbin+binup));
      binup++;
    } else {
      integral += (plot->GetBinContent(maxbin-bindown))/plot->Integral();
      plotcp->SetBinContent(maxbin-bindown, plot->GetBinContent(maxbin-bindown));
      bindown++;
    }
  }

  return plotcp;
}

void PlotOne(TH1D **plots, TCanvas *canv) {
  canv->Clear();
  for (int i = 0; i < nFiles; ++i) {
    if (i == 0) plots[i]->Draw();
    else plots[i]->Draw("same");
  }
  canv->Print("hi.pdf");
}

void Draw(std::string drawcmd, int rebin=1) {

  double maxx = 0;
  double minx = 100000;
  // Find the minima and maxima for drawing
  /*
  // Find maximum and minimum for plots
  for (int i = 0; i < nFiles; ++i) {

    if (maxx < tree[i]->GetMaximum(drawcmd.c_str())) {
      maxx = tree[i]->GetMaximum(drawcmd.c_str());
    }

    if (minx > tree[i]->GetMinimum(drawcmd.c_str())) {
      minx = tree[i]->GetMinimum(drawcmd.c_str());
    }

  }

  if (maxx == minx) continue;
  */

  if (drawcmd == "mean") {
    minx = 1-0.5;
    maxx = 1+0.5;
  } else if (drawcmd == "error_mean") {
    minx = 0.;
    maxx = 0.5;
  } else if (drawcmd == "sigma") {
    minx = 0;
    maxx = 1.0;
  } else if (drawcmd == "error_sigma") {
    minx = 0.;
    maxx = 0.5;
  } else if (drawcmd == "minimum_llh") {
    minx = 0;
    maxx = 350;
  }

  canv->Clear();
  // Draw the means
  double max = 0;
  std::string BurnIn = "status==0 && edm < 0.1";
  //std::string BurnIn;
  for (int i = 0; i < nFiles; ++i) {
    file[i]->cd();
    tree[i]->Draw(Form("%s>>%s_%i(200, %f, %f)", drawcmd.c_str(), drawcmd.c_str(), i, minx, maxx), BurnIn.c_str(), "");
    plots[i] = (TH1D*)gDirectory->Get(Form("%s_%i", drawcmd.c_str(), i));
    plots[i]->SetDirectory(0);
    plots[i]->SetTitle(plots[i]->GetName());
    plots[i]->Rebin(rebin);
    if (tree[i]->GetMinimum(drawcmd.c_str()) != tree[i]->GetMaximum(drawcmd.c_str())) {
      if (plots[i]->GetMaximum() > max) max = plots[i]->GetMaximum();
    }
  }
  plots[0]->SetTitle(drawcmd.c_str());

  TLine *line = new TLine();
  line->SetY1(0);
  line->SetY2(1.2*max);
  if (drawcmd == "mean") {
    line->SetX1(0);
    line->SetX2(0);
  } else if (drawcmd == "sigma") {
    line->SetX1(1);
    line->SetX2(1);
  }
  line->SetLineWidth(3);
  line->SetLineColor(kRed);
  line->SetLineStyle(kDashed);

  TLine *line2[nFiles];

  canv->Clear();
  canv->Divide(1, nFiles);
  TPaveText *text[nFiles];
  for (int i = 0; i < nFiles; ++i) {
    canv->cd(i+1);
    (canv->cd(i+1))->SetLeftMargin((canv->cd(i+1))->GetLeftMargin()*0.6);
    (canv->cd(i+1))->SetRightMargin((canv->cd(i+1))->GetRightMargin()*0.4);
    (canv->cd(i+1))->SetBottomMargin((canv->cd(i+1))->GetBottomMargin()*0.6);
    (canv->cd(i+1))->SetTopMargin((canv->cd(i+1))->GetTopMargin()*0.2);
    SetStyle(plots[i], i);
    // Make a little stats box
    text[i] = new TPaveText(0.7, 0.2, 0.95, 0.9, "NDC");
    text[i]->AddText(0, 0, plots[i]->GetTitle());
    text[i]->AddText(0, 0, Form("#mu=%.3f", plots[i]->GetMean()));
    text[i]->AddText(0, 0, Form("#sigma=%.3f", plots[i]->GetRMS()));
    plots[i]->Draw();
    plots[i]->SetTitle("");
    // Make the filler plot too
    FillMe(plots[i])->Draw("same");
    (canv->cd(i+1))->SetGridy(false);
    line->Draw("same");
    line2[i] = new TLine(plots[i]->GetMean(), 0, plots[i]->GetMean(), 1.2*max);
    line2[i]->SetX2(plots[i]->GetMean());
    line2[i]->SetLineColor(kWhite);
    line2[i]->SetLineWidth(3);
    line2[i]->Draw("same");
    text[i]->SetBorderSize(0);
    text[i]->Draw("same");
  }
  canv->Print(Form("%s.pdf", output.c_str()));
  canv->Clear();
  delete line;
  for (int i = 0; i < nFiles; ++i) {
    delete text[i];
    delete line2[i];
  }
}

void plotem(double ndata, double scale) {

  std::cout << "Setting up the files" << std::endl;
  for (int i = 0; i < nFiles; ++i) {
    file[i] = new TFile(Form("MCStats_minuit2_fittype_%i_nData%2.2f_MCscale%1.1f_muSigmaOnly_lots.root", i, ndata, scale));
  }

  // Get the TTree
  for (int i = 0; i < nFiles; ++i) {
    file[i]->cd();
    tree[i] = (TTree*)file[i]->Get("output")->Clone(Form("tree_%i", i));
  }

  canv = new TCanvas("canv", "canv", 1024, 1024);
  output = std::string(file[0]->GetName());
  output = output.substr(0, output.find(".root"));
  canv->Print(Form("%s.pdf[", output.c_str()));

  CalcCoverage(); 

  Draw(std::string("mean"), 1);
  Draw(std::string("error_mean"), 1);
  Draw(std::string("sigma"), 1);
  Draw(std::string("error_sigma"), 1);
  Draw(std::string("minimum_llh"), 1);

  canv->Print(Form("%s.pdf]", output.c_str()));
}
