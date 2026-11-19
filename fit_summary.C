void DoRatioPlot(TH1* , TH1* , string , string, string);

template <size_t N>
void DoRatioPlot2(TH1D* (&hists)[N], string (&names)[N], string path);

void fit_summary(string production = "train_801453/unambiguous") {
    // trains: 797807, 798097, 802001, 801453

    vector<double> ptMins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0};
    vector<double> ptMaxs = {0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};
    vector<double> ptBins = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0};

    string fInList[] = {
        Form("streaming_request/%s/pt_0.0_0.5/output__hDimuonMassUnamb_pt_0.0_0.5_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_0.5_1.0/output__hDimuonMassUnamb_pt_0.5_1.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_1.0_1.5/output__hDimuonMassUnamb_pt_1.0_1.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_1.5_2.0/output__hDimuonMassUnamb_pt_1.5_2.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_2.0_2.5/output__hDimuonMassUnamb_pt_2.0_2.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_2.5_3.0/output__hDimuonMassUnamb_pt_2.5_3.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_3.0_3.5/output__hDimuonMassUnamb_pt_3.0_3.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_3.5_4.0/output__hDimuonMassUnamb_pt_3.5_4.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_4.0_4.5/output__hDimuonMassUnamb_pt_4.0_4.5_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_4.5_5.0/output__hDimuonMassUnamb_pt_4.5_5.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_5.0_6.0/output__hDimuonMassUnamb_pt_5.0_6.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_6.0_7.0/output__hDimuonMassUnamb_pt_6.0_7.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_7.0_8.0/output__hDimuonMassUnamb_pt_7.0_8.0_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_8.0_10.0/output__hDimuonMassUnamb_pt_8.0_10.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_10.0_12.0/output__hDimuonMassUnamb_pt_10.0_12.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_12.0_15.0/output__hDimuonMassUnamb_pt_12.0_15.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str()),
        Form("streaming_request/%s/pt_15.0_20.0/output__hDimuonMassUnamb_pt_15.0_20.0_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails.root", production.c_str()),
    };

    TH1D *hJpsiMean = new TH1D("hJpsiMean", "J/#psi mean;#it{p}_{T} (GeV/#it{c});#mu_{J/#psi}", ptBins.size()-1, ptBins.data());
    TH1D *hJpsiWidth = new TH1D("hJpsiWidth", "J/#psi width;#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi}", ptBins.size()-1, ptBins.data());
    TH1D *hJpsiYield = new TH1D("hJpsiYield", "J/#psi yield;#it{p}_{T} (GeV/#it{c});Y_{J/#psi}", ptBins.size()-1, ptBins.data());

    vector<double> jpsiMeans, jpsiWidths, jpsiYields, jpsiYieldErrs;
    for (int i = 0; i < sizeof(fInList)/sizeof(fInList[0]); i++) {
        TFile *fIn = TFile::Open(fInList[i].c_str(), "READ");
        fIn->ls();
        TH1D *hFitResults = (TH1D*) fIn->Get(Form("fit_results_CB2_CB2_VWG__2.1_4.9_0.584_1.05__data_tails_hDimuonMassUnamb_pt_%2.1f_%2.1f", ptMins[i], ptMaxs[i]));
        if (!hFitResults) {
            hFitResults = (TH1D*) fIn->Get(Form("fit_results_CB2_CB2_RatioVWG__2.1_4.9_0.584_1.08__data_tails_hDimuonMassUnamb_pt_%2.1f_%2.1f", ptMins[i], ptMaxs[i]));
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
}

void compare(string prod1 = "train_797807", string prod2 = "train_802001") {
    string paths[] = {"train_802001", "train_801453/all", "train_801453/unambiguous"};
    string names[] = {"LHC26ac_pass1_fwdalign", "LHC26ac_pass1 + GMM (all)", "LHC26ac_pass1 + GMM (unamb)"};
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

    DoRatioPlot2(hJpsiMeans, names, "streaming_request/mean_comparison.pdf");
    DoRatioPlot2(hJpsiWidths, names, "streaming_request/width_comparison.pdf");
    DoRatioPlot2(hJpsiYields, names, "streaming_request/yield_comparison.pdf");
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

    for (int i = 0;i < nHists;i++) {
        hists[i]->SetLineColor(colors[i]);
        hists[i]->SetMarkerColor(colors[i]);
        hists[i]->SetLineWidth(2);

    }

    hists[0]->GetXaxis()->SetLabelSize(0);
    hists[0]->GetXaxis()->SetTitleSize(0);

    TLegend* leg = new TLegend(0.35, 0.75, 0.95, 0.88);
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