void DoRatioPlot(TH1* , TH1* , string , string, string);

template <size_t N>
void DoRatioPlot2(TH1D* (&hists)[N], string (&names)[N], string path);

void fit_summary(string production = "train_801453/unambiguous") {
    // trains: 797807, 798097, 802001, 801453, 806703, 806704

    string histName = "hDimuonMass";
    if (production.find("unambiguous") != std::string::npos) {
        histName = "hDimuonMassUnamb";
    }

    vector<double> ptMins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0};
    vector<double> ptMaxs = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    vector<double> ptBins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    string fInList[] = {
        Form("streaming_request/%s/pt_0.0_0.5/output__%s_pt_0.0_0.5_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_0.5_1.0/output__%s_pt_0.5_1.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_1.0_1.5/output__%s_pt_1.0_1.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_1.5_2.0/output__%s_pt_1.5_2.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_2.0_2.5/output__%s_pt_2.0_2.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_2.5_3.0/output__%s_pt_2.5_3.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_3.0_3.5/output__%s_pt_3.0_3.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_3.5_4.0/output__%s_pt_3.5_4.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_4.0_4.5/output__%s_pt_4.0_4.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_4.5_5.0/output__%s_pt_4.5_5.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_5.0_6.0/output__%s_pt_5.0_6.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_6.0_7.0/output__%s_pt_6.0_7.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_7.0_8.0/output__%s_pt_7.0_8.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_8.0_10.0/output__%s_pt_8.0_10.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_10.0_12.0/output__%s_pt_10.0_12.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_12.0_15.0/output__%s_pt_12.0_15.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/pt_15.0_20.0/output__%s_pt_15.0_20.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str(), histName.c_str()),
    };

    TH1D *hJpsiMean = new TH1D("hJpsiMean", "J/#psi mean;#it{p}_{T} (GeV/#it{c});#mu_{J/#psi}", ptBins.size()-1, ptBins.data());
    TH1D *hJpsiWidth = new TH1D("hJpsiWidth", "J/#psi width;#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi}", ptBins.size()-1, ptBins.data());
    TH1D *hJpsiYield = new TH1D("hJpsiYield", "J/#psi yield;#it{p}_{T} (GeV/#it{c});Y_{J/#psi}", ptBins.size()-1, ptBins.data());

    vector<double> jpsiMeans, jpsiWidths, jpsiYields, jpsiYieldErrs;
    for (int i = 0; i < sizeof(fInList)/sizeof(fInList[0]); i++) {
        TFile *fIn = TFile::Open(fInList[i].c_str(), "READ");
        fIn->ls();
        TH1D *hFitResults = (TH1D*) fIn->Get(Form("fit_results_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails_%s_pt_%2.1f_%2.1f", histName.c_str(), ptMins[i], ptMaxs[i]));
        if (!hFitResults) {
            hFitResults = (TH1D*) fIn->Get(Form("fit_results_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails_%s_pt_%2.1f_%2.1f", histName.c_str(), ptMins[i], ptMaxs[i]));
        }

        hJpsiMean->SetBinContent(i+1, hFitResults->GetBinContent(hFitResults->GetXaxis()->FindBin("mean_Jpsi")));
        hJpsiMean->SetBinError(i+1, hFitResults->GetBinError(hFitResults->GetXaxis()->FindBin("mean_Jpsi")));
        hJpsiWidth->SetBinContent(i+1, hFitResults->GetBinContent(hFitResults->GetXaxis()->FindBin("width_Jpsi")));
        hJpsiWidth->SetBinError(i+1, hFitResults->GetBinError(hFitResults->GetXaxis()->FindBin("width_Jpsi")));
        hJpsiYield->SetBinContent(i+1, hFitResults->GetBinContent(hFitResults->GetXaxis()->FindBin("sig_Jpsi")));
        hJpsiYield->SetBinError(i+1, hFitResults->GetBinError(hFitResults->GetXaxis()->FindBin("sig_Jpsi")));
    }

    hJpsiMean->SetMarkerStyle(20);
    hJpsiMean->SetMarkerColor(kBlack);
    hJpsiMean->SetLineColor(kBlack);
    hJpsiWidth->SetMarkerStyle(20);
    hJpsiWidth->SetMarkerColor(kBlack);
    hJpsiWidth->SetLineColor(kBlack);

    hJpsiYield->Scale(1./hJpsiYield->Integral(), "WIDTH");
    hJpsiYield->SetMarkerStyle(20);
    hJpsiYield->SetMarkerColor(kBlack);
    hJpsiYield->SetLineColor(kBlack);

    TCanvas *canvas = new TCanvas("canvas", "Fit Summary", 1800, 600);
    canvas->Divide(3, 1);
    canvas->cd(1);
    hJpsiMean->GetYaxis()->SetRangeUser(3.0, 3.2);
    hJpsiMean->Draw("EP");
    canvas->cd(2);
    hJpsiWidth->GetYaxis()->SetRangeUser(0.05, 0.12);
    hJpsiWidth->Draw("EP");
    canvas->cd(3);
    gPad->SetLogy(true);
    hJpsiYield->Draw("EP");
    canvas->SaveAs(Form("streaming_request/%s/fit_summary.pdf", production.c_str()));

    TFile *fOut = TFile::Open(Form("streaming_request/%s/fit_summary.root", production.c_str()), "RECREATE");
    hJpsiMean->Write();
    hJpsiWidth->Write();
    hJpsiYield->Write();



    vector<double> ptMinsHighMass = {0.0, 2.0, 5.0};
    vector<double> ptMaxsHighMass = {2.0, 5.0, 10.0};
    vector<double> ptBinsHighMass = {0.0, 2.0, 5.0, 10.0};

    string fInListHighMass[] = {
        Form("streaming_request/%s/high_mass/pt_0.0_2.0/output__%s_pt_0.0_2.0_CB2_CB2_CB2_VWG__8_12_0.8949_1.0794__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/high_mass/pt_2.0_5.0/output__%s_pt_2.0_5.0_CB2_CB2_CB2_VWG__8_12_0.8949_1.0794__data_tails.root", production.c_str(), histName.c_str()),
        Form("streaming_request/%s/high_mass/pt_5.0_10.0/output__%s_pt_5.0_10.0_CB2_CB2_CB2_VWG__8_12_0.8949_1.0794__data_tails.root", production.c_str(), histName.c_str()),
    };

    TH1D *hUpsilon1sMean = new TH1D("hUpsilon1sMean", "J/#psi mean;#it{p}_{T} (GeV/#it{c});#mu_{J/#psi}", ptBinsHighMass.size()-1, ptBinsHighMass.data());
    TH1D *hUpsilon1sWidth = new TH1D("hUpsilon1sWidth", "J/#psi width;#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi}", ptBinsHighMass.size()-1, ptBinsHighMass.data());
    TH1D *hUpsilon1sYield = new TH1D("hUpsilon1sYield", "J/#psi yield;#it{p}_{T} (GeV/#it{c});Y_{J/#psi}", ptBinsHighMass.size()-1, ptBinsHighMass.data());

    vector<double> Upsilon1sMeans, Upsilon1sWidths, Upsilon1sYields, Upsilon1sYieldErrs;
    for (int i = 0; i < sizeof(fInListHighMass)/sizeof(fInListHighMass[0]); i++) {
        TFile *fIn = TFile::Open(fInListHighMass[i].c_str(), "READ");
        fIn->ls();
        TH1D *hFitResultsHighMass = (TH1D*) fIn->Get(Form("fit_results_CB2_CB2_CB2_VWG__8_12_0.8949_1.0794__data_tails_%s_pt_%2.1f_%2.1f", histName.c_str(), ptMinsHighMass[i], ptMaxsHighMass[i]));
        hUpsilon1sMean->SetBinContent(i+1, hFitResultsHighMass->GetBinContent(hFitResultsHighMass->GetXaxis()->FindBin("mean_Upsilon1s")));
        hUpsilon1sMean->SetBinError(i+1, hFitResultsHighMass->GetBinError(hFitResultsHighMass->GetXaxis()->FindBin("mean_Upsilon1s")));
        hUpsilon1sWidth->SetBinContent(i+1, hFitResultsHighMass->GetBinContent(hFitResultsHighMass->GetXaxis()->FindBin("width_Upsilon1s")));
        hUpsilon1sWidth->SetBinError(i+1, hFitResultsHighMass->GetBinError(hFitResultsHighMass->GetXaxis()->FindBin("width_Upsilon1s")));
        hUpsilon1sYield->SetBinContent(i+1, hFitResultsHighMass->GetBinContent(hFitResultsHighMass->GetXaxis()->FindBin("sig_Upsilon1s")));
        hUpsilon1sYield->SetBinError(i+1, hFitResultsHighMass->GetBinError(hFitResultsHighMass->GetXaxis()->FindBin("sig_Upsilon1s")));
    }

    hUpsilon1sMean->SetMarkerStyle(20);
    hUpsilon1sMean->SetMarkerColor(kBlack);
    hUpsilon1sMean->SetLineColor(kBlack);
    hUpsilon1sWidth->SetMarkerStyle(20);
    hUpsilon1sWidth->SetMarkerColor(kBlack);
    hUpsilon1sWidth->SetLineColor(kBlack);

    hUpsilon1sYield->Scale(1./hUpsilon1sYield->Integral(), "WIDTH");
    hUpsilon1sYield->SetMarkerStyle(20);
    hUpsilon1sYield->SetMarkerColor(kBlack);
    hUpsilon1sYield->SetLineColor(kBlack);

    TCanvas *canvasHighMass = new TCanvas("canvasHighMass", "Fit Summary", 1800, 600);
    canvasHighMass->Divide(3, 1);
    canvasHighMass->cd(1);
    hUpsilon1sMean->GetYaxis()->SetRangeUser(9.0, 10);
    hUpsilon1sMean->Draw("EP");
    canvasHighMass->cd(2);
    hUpsilon1sWidth->GetYaxis()->SetRangeUser(0.05, 0.2);
    hUpsilon1sWidth->Draw("EP");
    canvasHighMass->cd(3);
    gPad->SetLogy(true);
    hUpsilon1sYield->Draw("EP");
    canvasHighMass->SaveAs(Form("streaming_request/%s/fit_summary_high_mass.pdf", production.c_str()));

    TFile *fOutHighMass = TFile::Open(Form("streaming_request/%s/fit_summary_high_mass.root", production.c_str()), "RECREATE");
    hUpsilon1sMean->Write();
    hUpsilon1sWidth->Write();
    hUpsilon1sYield->Write();
}

void compare(string prod1 = "train_797807", string prod2 = "train_802001") {
    string paths[] = {"train_806704/all", "train_806704/unambiguous", "train_806703/all", "train_806703/unambiguous"};
    string names[] = {"LHC26ac_pass1_fwdalign (all)", "LHC26ac_pass1_fwdalign (unamb)", "LHC26ac_pass1 + GMM (all)", "LHC26ac_pass1 + GMM (unamb)"};
    const int nTrains = sizeof(paths)/sizeof(paths[0]);
    TFile *fIns[nTrains];
    TH1D *hJpsiMeans[nTrains], *hJpsiWidths[nTrains], *hJpsiYields[nTrains];
    for (int i = 0;i < nTrains;i++) {
        std::cout << Form("streaming_request/%s/fit_summary.root", paths[i].c_str()) << std::endl;
        fIns[i] = TFile::Open(Form("streaming_request/%s/fit_summary.root", paths[i].c_str()));
        hJpsiMeans[i] = (TH1D*) fIns[i]->Get("hJpsiMean");
        hJpsiWidths[i] = (TH1D*) fIns[i]->Get("hJpsiWidth");
        hJpsiYields[i] = (TH1D*) fIns[i]->Get("hJpsiYield");
    }

    DoRatioPlot2(hJpsiMeans, names, "streaming_request/mean_jpsi_comparison.pdf");
    DoRatioPlot2(hJpsiWidths, names, "streaming_request/width_jpsi_comparison.pdf");
    DoRatioPlot2(hJpsiYields, names, "streaming_request/yield_jpsi_comparison.pdf");


    TH1D *hUpsilon1sMeans[nTrains], *hUpsilon1sWidths[nTrains], *hUpsilon1sYields[nTrains];
    for (int i = 0;i < nTrains;i++) {
        std::cout << Form("streaming_request/%s/fit_summary_high_mass.root", paths[i].c_str()) << std::endl;
        fIns[i] = TFile::Open(Form("streaming_request/%s/fit_summary_high_mass.root", paths[i].c_str()));
        hUpsilon1sMeans[i] = (TH1D*) fIns[i]->Get("hUpsilon1sMean");
        hUpsilon1sMeans[i]->GetYaxis()->SetTitle("#mu_{#varUpsilon(1S)}");
        hUpsilon1sMeans[i]->SetTitle("#varUpsilon(1S) mean");
        hUpsilon1sWidths[i] = (TH1D*) fIns[i]->Get("hUpsilon1sWidth");
        hUpsilon1sWidths[i]->GetYaxis()->SetRangeUser(0.05, 0.30);
        hUpsilon1sWidths[i]->GetYaxis()->SetTitle("#sigma_{#varUpsilon(1S)}");
        hUpsilon1sWidths[i]->SetTitle("#varUpsilon(1S) width");
        hUpsilon1sYields[i] = (TH1D*) fIns[i]->Get("hUpsilon1sYield");
        hUpsilon1sYields[i]->GetYaxis()->SetRangeUser(0.0, 0.30);
        hUpsilon1sYields[i]->GetYaxis()->SetTitle("Y_{#varUpsilon(1S)}");
        hUpsilon1sYields[i]->SetTitle("#varUpsilon(1S) yield");
    }

    DoRatioPlot2(hUpsilon1sMeans, names, "streaming_request/mean_upsilon1s_comparison.pdf");
    DoRatioPlot2(hUpsilon1sWidths, names, "streaming_request/width_upsilon1s_comparison.pdf");
    DoRatioPlot2(hUpsilon1sYields, names, "streaming_request/yield_upsilon1s_comparison.pdf");
}

template <size_t N>
void DoRatioPlot2(TH1D* (&hists)[N], string (&names)[N], string path) {

    Color_t colors[] = {kBlack, kRed+1, kAzure+4, kGreen+1};
    //const int nHists = sizeof(&hists)/sizeof(&(hists[0]));
    const int nHists = N;
    for (int i = 0; i < nHists; i++) {
        if (!hists[i]) {
            std::cerr << "Errore: istogrammi non validi!" << std::endl;
            return;
        }
    }

    // Canvas
    TCanvas* c = new TCanvas("c", "Histogram comparison", 800, 800);
    gStyle->SetOptStat(false);

    // Pads
    TPad* pad1 = new TPad("pad1", "Histograms", 0, 0.30, 1, 1);
    TPad* pad2 = new TPad("pad2", "Ratio", 0, 0, 1, 0.30);

    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);

    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.30);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.05);

    pad1->Draw();
    pad2->Draw();

    // -------------------------
    // Upper pad: histograms
    // -------------------------
    pad1->cd();
    gStyle->SetOptStat(false);
    //gPad->SetLogy(true);
    gStyle->SetOptStat(false);

    for (int i = 0;i < nHists;i++) {
        hists[i]->SetLineColor(colors[i]);
        hists[i]->SetMarkerColor(colors[i]);
        hists[i]->SetLineWidth(2);

    }

    hists[0]->GetXaxis()->SetLabelSize(0);
    hists[0]->GetXaxis()->SetTitleSize(0);

    TLegend* leg = new TLegend(0.35, 0.65, 0.95, 0.88);
    for (int i = 0;i < nHists;i++) {
        hists[i]->Draw(i==0 ? "E" : "E SAME");
        leg->AddEntry(hists[i], names[i].c_str(), "lep");
    }
    leg->SetTextSize(0.045);
    leg->Draw();

    // -------------------------
    // Lower pad: ratio h1 / h2
    // -------------------------
    pad2->cd();
    gStyle->SetOptStat(false);

    TH1D* hRatios[nHists];
    for (int i = 1;i < nHists;i++) {
        hRatios[i-1] = (TH1D*) hists[i]->Clone(Form("hRatio_%d", i));
        hRatios[i-1]->SetDirectory(nullptr);
        hRatios[i-1]->Divide(hists[0]);

        hRatios[i-1]->SetTitle("");
        hRatios[i-1]->SetLineColor(colors[i]);
        hRatios[i-1]->SetMarkerColor(colors[i]);
        hRatios[i-1]->SetMarkerStyle(20);

        hRatios[i-1]->GetYaxis()->SetTitle("Ratio");
        hRatios[i-1]->GetYaxis()->SetNdivisions(505);
        hRatios[i-1]->GetYaxis()->SetTitleSize(0.10);
        hRatios[i-1]->GetYaxis()->SetTitleOffset(0.45);
        hRatios[i-1]->GetYaxis()->SetLabelSize(0.09);

        hRatios[i-1]->GetXaxis()->SetTitle(hists[0]->GetXaxis()->GetTitle());
        hRatios[i-1]->GetXaxis()->SetTitleSize(0.12);
        hRatios[i-1]->GetXaxis()->SetLabelSize(0.10);
        hRatios[i-1]->GetXaxis()->SetTitleOffset(1.0);

        hRatios[i-1]->SetMinimum(0.8);
        hRatios[i-1]->SetMaximum(1.2);

        hRatios[i-1]->Draw(i==0 ? "E" : "E SAME");
    }

    // Reference line at ratio = 1
    TLine* line = new TLine(
        hRatios[0]->GetXaxis()->GetXmin(), 1.0,
        hRatios[0]->GetXaxis()->GetXmax(), 1.0
    );

    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw("SAME");

    c->cd();
    c->Update();
    c->SaveAs(Form("%s.pdf", path.c_str()));
}






void DoRatioPlot(TH1* h1, TH1* h2, string nameH1, string nameH2, string path) {
    if (!h1 || !h2) {
        std::cerr << "Errore: istogrammi non validi!" << std::endl;
        return;
    }

    // Canvas
    TCanvas* c = new TCanvas("c", "Histogram comparison", 800, 800);

    // Pads
    TPad* pad1 = new TPad("pad1", "Histograms", 0, 0.30, 1, 1);
    TPad* pad2 = new TPad("pad2", "Ratio", 0, 0, 1, 0.30);

    pad1->SetBottomMargin(0.02);
    pad1->SetLeftMargin(0.12);
    pad1->SetRightMargin(0.05);

    pad2->SetTopMargin(0.03);
    pad2->SetBottomMargin(0.30);
    pad2->SetLeftMargin(0.12);
    pad2->SetRightMargin(0.05);

    pad1->Draw();
    pad2->Draw();

    // -------------------------
    // Upper pad: histograms
    // -------------------------
    pad1->cd();
    gStyle->SetOptStat(false);

    h1->SetLineColor(kRed);
    h1->SetMarkerColor(kRed);
    h1->SetLineWidth(2);

    h2->SetLineColor(kBlue);
    h2->SetMarkerColor(kBlue);
    h2->SetLineWidth(2);

    h1->GetXaxis()->SetLabelSize(0);
    h1->GetXaxis()->SetTitleSize(0);

    h1->Draw("E");
    h2->Draw("E SAME");

    TLegend* leg = new TLegend(0.45, 0.75, 0.95, 0.88);
    leg->AddEntry(h1, nameH1.c_str(), "lep");
    leg->AddEntry(h2, nameH2.c_str(), "lep");
    leg->SetTextSize(0.045);
    leg->Draw();

    // -------------------------
    // Lower pad: ratio h1 / h2
    // -------------------------
    pad2->cd();
    gStyle->SetOptStat(false);

    TH1* hRatio = (TH1*)h1->Clone("hRatio");
    hRatio->SetDirectory(nullptr);
    hRatio->Divide(h2);

    hRatio->SetTitle("");

    hRatio->SetLineColor(kBlack);
    hRatio->SetMarkerColor(kBlack);
    hRatio->SetMarkerStyle(20);

    hRatio->GetYaxis()->SetTitle("Ratio");
    hRatio->GetYaxis()->SetNdivisions(505);
    hRatio->GetYaxis()->SetTitleSize(0.10);
    hRatio->GetYaxis()->SetTitleOffset(0.45);
    hRatio->GetYaxis()->SetLabelSize(0.09);

    hRatio->GetXaxis()->SetTitle(h1->GetXaxis()->GetTitle());
    hRatio->GetXaxis()->SetTitleSize(0.12);
    hRatio->GetXaxis()->SetLabelSize(0.10);
    hRatio->GetXaxis()->SetTitleOffset(1.0);

    hRatio->SetMinimum(0.8);
    hRatio->SetMaximum(1.2);

    hRatio->Draw("E");

    // Reference line at ratio = 1
    TLine* line = new TLine(
        hRatio->GetXaxis()->GetXmin(), 1.0,
        hRatio->GetXaxis()->GetXmax(), 1.0
    );

    line->SetLineStyle(2);
    line->SetLineColor(kRed);
    line->Draw("SAME");

    c->cd();
    c->Update();
    c->SaveAs(Form("%s.pdf", path.c_str()));
}