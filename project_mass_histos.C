void SetHist(TH1D *hist, Color_t color, int style, double alpha = 1) {
    hist->SetMarkerStyle(style);
    hist->SetMarkerColorAlpha(color, alpha);
    hist->SetMarkerSize(1);
    hist->SetLineColorAlpha(color, alpha);
    hist->Scale(1./hist->Integral());
}

void SetLegend(TLegend *legend){
    legend->SetBorderSize(0);
    legend->SetFillColor(10);
    legend->SetFillStyle(1);
    legend->SetLineStyle(0);
    legend->SetLineColor(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.045);
}

void project_mass_histos(double minPt = 0, double maxPt = 20) {
    string label = "LHC26ac_pass1_skimmed";
    string fInName = "streaming_request/train_801453/AnalysisResults.root";
    //string label = "LHC26ak_apass1";
    //string fInName = "tag_comparison/AnalysisResults_LHC26ak_apass1.root";
    string fOutName = Form("streaming_request/train_801453/dimuMass_pt_%2.1f_%2.1f_%s.root", minPt, maxPt, label.c_str());
    //string fOutName = Form("tag_comparison/dimuMass_pt_%1.0f_%1.0f_%s.root", minPt, maxPt, label.c_str());
    string muonCutTr = "muonQualityCutsMUONStandalone"; // matchedMchMid, muonQualityCutsMUONStandalone muonLowPt210SigmaPDCA
    string histBlock = "Mass_Pt_Cent_Rap"; // Mass_Pt_Cent_Rap Mass_Pt_Rapidity
    
    TFile *fIn = new TFile(fInName.c_str(), "READ");
    TList *hlistTrDimuons = (TList*) fIn->Get("analysis-same-event-pairing/output");
    TList *listTrDimuons = (TList*) hlistTrDimuons->FindObject(Form("PairsMuonSEPM_%s", muonCutTr.c_str()));
    THnSparseD *hSparseDimuon = (THnSparseD*) listTrDimuons->FindObject(histBlock.c_str());
    hSparseDimuon->GetAxis(1)->SetRangeUser(minPt, maxPt);
    hSparseDimuon->GetAxis(3)->SetRangeUser(2.5, 4); // switch to 2 in case of Mass_Pt_Rapidity
    TH1D *hDimuonMass = (TH1D*) hSparseDimuon->Projection(0, "hDimuonMass");
    hDimuonMass->Draw();

    TFile *fOut = new TFile(fOutName.c_str(), "RECREATE");
    hDimuonMass->Write("hDimuonMass");
    fOut->Close();
}

void run_projections() {
    //double minPt[] = {0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 10.0};
    //double maxPt[] = {20.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 6.0, 7.0, 8.0, 10.0, 12.0, 15.0, 20.0, 20.0};

    double minPt[] = {0.0, 0.0, 2.0, 5.0};
    double maxPt[] = {20.0, 2.0, 5.0, 10.0};

    string label = "LHC26ac_pass1_skimmed";
    string fInName = "streaming_request/train_806703/AnalysisResults.root";
    string fOutName = Form("streaming_request/train_806703/dimuHighMass_%s.root", label.c_str());
    //string fOutName = Form("tag_comparison/dimuMass_pt_%1.0f_%1.0f_%s.root", minPt, maxPt, label.c_str());
    string muonCutTr = "muonQualityCutsMUONStandalone"; // matchedMchMid, muonQualityCutsMUONStandalone muonLowPt210SigmaPDCA
    string histBlock = "Mass_Pt_Rapidity"; // Mass_Pt_Cent_Rap Mass_Pt_Rapidity
        
    TFile *fIn = new TFile(fInName.c_str(), "READ");

    TFile *fOut = new TFile(fOutName.c_str(), "RECREATE");
    for (int i = 0; i < sizeof(minPt)/sizeof(minPt[0]); i++) {
        std::cout << minPt[i] << " - " << maxPt[i] << std::endl;
        //project_mass_histos(minPt[i], maxPt[i]);

        ////////////
        TList *hlistTrDimuons = (TList*) fIn->Get("analysis-same-event-pairing/output");
        TList *listTrDimuons = (TList*) hlistTrDimuons->FindObject(Form("PairsMuonSEPM_%s", muonCutTr.c_str()));
        TList *listTrDimuonsUnamb = (TList*) hlistTrDimuons->FindObject(Form("PairsMuonSEPM_unambiguous_%s", muonCutTr.c_str()));
        THnSparseD *hSparseDimuon = (THnSparseD*) listTrDimuons->FindObject(histBlock.c_str());
        THnSparseD *hSparseDimuonUnamb = (THnSparseD*) listTrDimuonsUnamb->FindObject(histBlock.c_str());

        hSparseDimuon->GetAxis(1)->SetRangeUser(minPt[i], maxPt[i]);
        if (histBlock == "Mass_Pt_Rapidity") {
            hSparseDimuon->GetAxis(2)->SetRangeUser(2.5, 4); // switch to 2 in case of Mass_Pt_Rapidity
        } else {
            hSparseDimuon->GetAxis(3)->SetRangeUser(2.5, 4);
        }
        TH1D *hDimuonMass = (TH1D*) hSparseDimuon->Projection(0, Form("hDimuonMass_pt_%2.1f_%2.1f", minPt[i], maxPt[i]));

        hSparseDimuonUnamb->GetAxis(1)->SetRangeUser(minPt[i], maxPt[i]);
        if (histBlock == "Mass_Pt_Rapidity") {
            hSparseDimuonUnamb->GetAxis(2)->SetRangeUser(2.5, 4); // switch to 2 in case of Mass_Pt_Rapidity
        } else {
            hSparseDimuonUnamb->GetAxis(3)->SetRangeUser(2.5, 4);
        }
        TH1D *hDimuonMassUnamb = (TH1D*) hSparseDimuonUnamb->Projection(0, Form("hDimuonMassUnamb_pt_%2.1f_%2.1f", minPt[i], maxPt[i]));

        fOut->cd();
        hDimuonMass->Write(Form("hDimuonMass_pt_%2.1f_%2.1f", minPt[i], maxPt[i]));
        hDimuonMassUnamb->Write(Form("hDimuonMassUnamb_pt_%2.1f_%2.1f", minPt[i], maxPt[i]));
        ////////////
    }
    fOut->Close();
}

void compare_histos(double minPt = 0, double maxPt = 2) {
    TFile *fIn1 = new TFile("~/GITHUB/dq_streaming_qa/streaming_request/train_731427/QA.root", "READ");
    TH1D *hDimuonMass1 = (TH1D*) fIn1->Get("hDimuonMass");
    TH1D *hPt1 = (TH1D*) fIn1->Get("hPt");
    TH1D *hEta1 = (TH1D*) fIn1->Get("hEta");
    TH1D *hPhi1 = (TH1D*) fIn1->Get("hPhi");
    TH1D *hPdca1 = (TH1D*) fIn1->Get("hPdca");
    TH1D *hRabs1 = (TH1D*) fIn1->Get("hRabs");
    TH1D *hChi2MCHMID1 = (TH1D*) fIn1->Get("hChi2MCHMID");

    SetHist(hDimuonMass1, kRed, 20);
    SetHist(hPt1, kRed, 20);
    SetHist(hEta1, kRed, 20);
    SetHist(hPhi1, kRed, 20);
    SetHist(hPdca1, kRed, 20);
    SetHist(hRabs1, kRed, 20);
    SetHist(hChi2MCHMID1, kRed, 20);

    TFile *fIn2 = new TFile("~/GITHUB/dq_streaming_qa/streaming_request/train_731428/QA.root", "READ");
    TH1D *hDimuonMass2 = (TH1D*) fIn2->Get("hDimuonMass");
    TH1D *hPt2 = (TH1D*) fIn2->Get("hPt");
    TH1D *hEta2 = (TH1D*) fIn2->Get("hEta");
    TH1D *hPhi2 = (TH1D*) fIn2->Get("hPhi");
    TH1D *hPdca2 = (TH1D*) fIn2->Get("hPdca");
    TH1D *hRabs2 = (TH1D*) fIn2->Get("hRabs");
    TH1D *hChi2MCHMID2 = (TH1D*) fIn2->Get("hChi2MCHMID");

    SetHist(hDimuonMass2, kBlue, 24);
    SetHist(hPt2, kBlue, 24);
    SetHist(hEta2, kBlue, 24);
    SetHist(hPhi2, kBlue, 24);
    SetHist(hPdca2, kBlue, 24);
    SetHist(hRabs2, kBlue, 24);
    SetHist(hChi2MCHMID2, kBlue, 24);

    auto canvasDimuonMass = new TCanvas("canvasDimuonMass", "A ratio example");
    gStyle->SetOptStat(false);
    hDimuonMass1->SetTitle("");
    hDimuonMass1->GetXaxis()->SetTitle("#it{M}_{#mu#mu} (GeV/#it{c}^{2})");
    auto rp = new TRatioPlot(hDimuonMass1, hDimuonMass2);
    canvasDimuonMass->SetTicks(0, 1);
    rp->Draw();
    rp->GetUpperPad()->SetLogy();
    rp->GetLowerRefGraph()->SetMinimum(0.5);
    rp->GetLowerRefGraph()->SetMaximum(1.5);
    rp->GetLowYaxis()->SetNdivisions(505);

    rp->GetUpperPad()->cd();
    TLegend *legendDimuonMass = new TLegend(0.45, 0.7, 0.85, 0.85);
    SetLegend(legendDimuonMass);
    legendDimuonMass->AddEntry(hDimuonMass1, "Standard association", "P");
    legendDimuonMass->AddEntry(hDimuonMass2, "Time association", "P");
    legendDimuonMass->Draw();

    canvasDimuonMass->Modified();
    canvasDimuonMass->Update();

    TCanvas *CanvasMuonKine = new TCanvas("CanvasMuonKine", "", 1800, 1200);
    CanvasMuonKine->Divide(3,2);
    CanvasMuonKine->cd(1);
    gPad->SetLogy(true);
    hPt1->Draw("EP");
    hPt2->Draw("EP SAME");

    TLegend *legendMuonKine = new TLegend(0.45, 0.7, 0.85, 0.85);
    legendMuonKine->AddEntry(hPt1, "Standard assoc.", "P");
    legendMuonKine->AddEntry(hPt2, "Time assoc.", "P");
    legendMuonKine->Draw();

    CanvasMuonKine->cd(2);
    gPad->SetLogy(true);
    hEta1->Draw("EP");
    hEta2->Draw("EP SAME");
    CanvasMuonKine->cd(3);
    
    gPad->SetLogy(true);
    hPhi1->Draw("EP");
    hPhi2->Draw("EP SAME");

    CanvasMuonKine->cd(4);
    gPad->SetLogy(true);
    hPdca1->Draw("EP");
    hPdca2->Draw("EP SAME");

    CanvasMuonKine->cd(5);
    gPad->SetLogy(true);
    hRabs1->Draw("EP");
    hRabs2->Draw("EP SAME");

    CanvasMuonKine->cd(6);
    gPad->SetLogy(true);
    hChi2MCHMID1->Draw("EP");
    hChi2MCHMID2->Draw("EP SAME");

    canvasDimuonMass->SaveAs("streaming_request/DimuonMass_LHC24.pdf");
    CanvasMuonKine->SaveAs("streaming_request/MuonKine_LHC24.pdf");

    //canvas->SaveAs(Form("~/GITHUB/dq_streaming_qa/tag_comparison/mass_spectrum_comparison_pt_%1.0f_%1.0f.pdf", minPt, maxPt));
}

void plot_fit_quality_variables() {
    double jpsiMean[] = {3.0877, 3.0880, 3.0921, 3.0877, 3.0834};
    double errJpsiMean[] = {0.0126, 0.0048, 0.0001, 0.0040, 0.0059};
    double jpsiWidth[] = {0.0671, 0.0691, 0.0706, 0.0693, 0.0795};
    double errJpsiWidth[] = {0.0390, 0.0059, 0.0026, 0.0040, 0.0078};
    double jpsiRawYield[] = {25009, 11599, 6017, 4209, 1855};
    double errJpsiRawYield[] = {5361, 898, 309, 247, 185};

    double ptBins[] = {0, 2, 3, 4, 6, 10};
    TH1D *hJpsiMean = new TH1D("hJpsiMean", ";#it{p}_{T} (GeV/#it{c});#mu_{J/#psi}", 5, ptBins);
    TH1D *hJpsiWidth = new TH1D("hJpsiWidth", ";#it{p}_{T} (GeV/#it{c});#sigma_{J/#psi}", 5, ptBins);
    TH1D *hJpsiRawYield = new TH1D("hJpsiRawYield", ";#it{p}_{T} (GeV/#it{c});d#it{N}_{J/#psi}/d#it{p}_{T}", 5, ptBins);

    SetHist(hJpsiMean, kRed+1, 20);
    SetHist(hJpsiWidth, kRed+1, 20);
    SetHist(hJpsiRawYield, kRed+1, 20);

    for (int iPt = 0;iPt < 5;iPt++) {
        hJpsiMean->SetBinContent(iPt+1, jpsiMean[iPt]);
        hJpsiMean->SetBinError(iPt+1, errJpsiMean[iPt]);

        hJpsiWidth->SetBinContent(iPt+1, jpsiWidth[iPt]);
        hJpsiWidth->SetBinError(iPt+1, errJpsiWidth[iPt]);

        hJpsiRawYield->SetBinContent(iPt+1, jpsiRawYield[iPt]);
        hJpsiRawYield->SetBinError(iPt+1, errJpsiRawYield[iPt]);
    }

    TCanvas *canvasJpsiMean = new TCanvas("canvasJpsiMean", "", 800, 600);
    gStyle->SetOptStat(false);
    hJpsiMean->GetYaxis()->SetRangeUser(2.6, 3.6);
    hJpsiMean->Draw();

    TCanvas *canvasJpsiWidth = new TCanvas("canvasJpsiWidth", "", 800, 600);
    gStyle->SetOptStat(false);
    hJpsiWidth->GetYaxis()->SetRangeUser(0, 0.15);
    hJpsiWidth->Draw();

    TCanvas *canvasJpsiRawYield = new TCanvas("canvasJpsiRawYield", "", 800, 600);
    gStyle->SetOptStat(false);
    gPad->SetLogy(true);
    hJpsiRawYield->Scale(1, "WIDTH");
    hJpsiRawYield->Draw();

    canvasJpsiMean->SaveAs("tag_comparison/jpsiMean.pdf");
    canvasJpsiWidth->SaveAs("tag_comparison/jpsiWidth.pdf");
    canvasJpsiRawYield->SaveAs("tag_comparison/jpsiRawYield.pdf");
}