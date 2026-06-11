#include <iostream>
#include <fstream>
#include <vector>
#include <filesystem>
#include <algorithm>
#include <string>

#include <TFile.h>
#include <TDirectoryFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TLorentzVector.h>
#include <TLegend.h>
#include <THnSparse.h>
#include <TCanvas.h>
#include <TLine.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TGrid.h>
#include <TH1.h>
#include <TH2.h>
#include <TF1.h>
#include <TRandom3.h>

using namespace std;

void SetLegend(TLegend *legend){
    legend->SetBorderSize(0);
    legend->SetFillColor(10);
    legend->SetFillStyle(1);
    legend->SetLineStyle(0);
    legend->SetLineColor(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.045);
}
void SetHist(TH1D *hist, Color_t color, int size, double alpha = 1) {
    hist->SetMarkerStyle(size);
    hist->SetMarkerColorAlpha(color, alpha);
    hist->SetMarkerSize(0.5);
    hist->SetLineColorAlpha(color, alpha);
    hist->Scale(1./hist->Integral());
}
TCanvas* DoPlot(TH1D *hist, string canvasName, string varName) {
    int firstBin = hist->FindFirstBinAbove(0.0);
    int lastBin  = hist->FindLastBinAbove(0.0);
    double xRangeMin = hist->GetBinLowEdge(firstBin);
    double xRangeMax = hist->GetBinCenter(lastBin) + hist->GetBinWidth(lastBin);

    TCanvas *canvas = new TCanvas(canvasName.c_str(), "", 800, 600);
    gPad->SetLogy(true);
    SetHist(hist, kBlack, 20);
    hist->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
    hist->SetTitle("");
    hist->SetName(canvasName.c_str());
    hist->Draw("EP");

    return canvas;
}
TCanvas* DoRatioPlot(TH1D *histSum, std::vector<TH1D*> hist, string canvasName, string varName) {
    const int nRuns = int(hist.size());
    gStyle->SetPalette(kRainBow);

    int firstBin = histSum->FindFirstBinAbove(0.0);
    int lastBin  = histSum->FindLastBinAbove(0.0);
    double xRangeMin = histSum->GetBinLowEdge(firstBin);
    double xRangeMax = histSum->GetBinCenter(lastBin) + histSum->GetBinWidth(lastBin);

    TCanvas *canvas = new TCanvas(canvasName.c_str(), "", 600, 800);

    TPad *pad1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    TPad *pad2 = new TPad("pad2", "pad2", 0, 0.0, 1, 0.3);

    pad1->SetBottomMargin(0);
    pad2->SetTopMargin(0);
    pad2->SetBottomMargin(0.3);

    pad1->Draw();
    pad2->Draw();

    //gPad->BuildLegend();
    pad1->cd();
    gPad->SetLogy(true);
    SetHist(histSum, kBlack, 20);
    histSum->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
    histSum->SetTitle("");
    histSum->SetName("Sum");
    histSum->Draw("EP");
    for (int iRun = 0;iRun < nRuns;iRun++) {
        hist[iRun]->Scale(1./hist[iRun]->Integral());
        hist[iRun]->Draw("SAME PLC PMC");
    }
    gPad->BuildLegend();

    pad2->cd();
    std::vector<TH1D*> histRatio(nRuns);
    for (int iRun = 0;iRun < nRuns;iRun++) {
        histRatio[iRun] = (TH1D*) hist[iRun]->Clone(Form("histRatio_%d", iRun));

        if (iRun == 0) {
            histRatio[iRun]->SetTitle("");
            histRatio[iRun]->GetYaxis()->SetTitle("Data/MC");
            histRatio[iRun]->GetXaxis()->SetTitle(varName.c_str());
            histRatio[iRun]->GetXaxis()->SetLabelSize(0.08);
            histRatio[iRun]->GetXaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetNdivisions(505);
            histRatio[iRun]->GetYaxis()->SetTitleSize(0.1);
            histRatio[iRun]->GetYaxis()->SetLabelSize(0.08);
        }

        histRatio[iRun]->SetTitle("");
        histRatio[iRun]->GetXaxis()->SetRangeUser(xRangeMin, xRangeMax);
        histRatio[iRun]->Divide(histSum);
        histRatio[iRun]->GetYaxis()->SetRangeUser(0.5, 1.5);
        histRatio[iRun]->Draw("SAME PLC PMC");
    }
    TLine *line = new TLine(xRangeMin, 1, xRangeMax, 1);
    line->SetLineStyle(2);
    line->Draw("SAME");

    canvas->cd();

    return canvas;
}

// process: table-maker / table-reader
// production: minBias / skimmed / heavyIon

void streaming_qa(string train = "686270", string process = "table-maker", string production = "minBias") {
    gStyle->SetPalette(kRainBow);
    string pathIn = Form("/Users/lucamicheletti/cernbox/DQ_streaming/train_%s", train.c_str());
    string muonCutTm = "muonMinimalCuts10SigmaPDCA";
    string muonCutTr = "muonQualityCutsMUONStandalone";
    string histBlock = "Mass_Pt_centrFT0C_V2";
    string triggerMask = "fSingleMuLow";

    TH1D *hSumVtxZ, *hSumCentFT0C, *hSumCentFT0M;
    TH2D *hSumPsi2A_CentFT0C, *hSumPsi2B_CentFT0C, *hSumPsi2C_CentFT0C;
    TH1D *hSumPt, *hSumEta, *hSumPhi, *hSumpdca, *hSumRAtAbsorberEnd, *hSumChi2MCHMID, *hSumChi2MCHMFT;

    if (process == "table-maker") {
        std::vector<int> runList;
        int entry;
        std::ifstream fInRunList(Form("%s/run_list.txt", pathIn.c_str()));
        while (fInRunList >> entry) {runList.push_back(entry);}
        const int nRuns = int(runList.size());
        std::cout << nRuns << std::endl;

        int iRun = 0;
        int nClasses = production == "skimmed" ? 11 : 6; // 9 = skimmed, 6 = minBias
        
        TH2D *hCountersVsRuns = new TH2D("hCountersVsRuns", ";Run;Counter", nRuns, 0, nRuns, nClasses, 0, nClasses);
        hCountersVsRuns->GetYaxis()->SetBinLabel(1, "BCs TVX");
        hCountersVsRuns->GetYaxis()->SetBinLabel(2, "BCs TVX after BC cuts");
        hCountersVsRuns->GetYaxis()->SetBinLabel(3, "Collisions All");
        hCountersVsRuns->GetYaxis()->SetBinLabel(4, "Collisions TVX");
        hCountersVsRuns->GetYaxis()->SetBinLabel(5, "Collisions Acc");
        hCountersVsRuns->GetYaxis()->SetBinLabel(6, "Collisions after cuts");
        if (production == "skimmed") {
            hCountersVsRuns->GetYaxis()->SetBinLabel(7, "inspected TVX");
            hCountersVsRuns->GetYaxis()->SetBinLabel(8, Form("%s scalers", triggerMask.c_str()));
            hCountersVsRuns->GetYaxis()->SetBinLabel(9, Form("%s selections", triggerMask.c_str()));
            hCountersVsRuns->GetYaxis()->SetBinLabel(10, Form("%s AnalysedTriggers", triggerMask.c_str()));
            hCountersVsRuns->GetYaxis()->SetBinLabel(11, Form("%s", triggerMask.c_str()));
        }

        std::vector<TH1D*> hVtxZ(nRuns);
        std::vector<TH1D*> hCentFT0C(nRuns);
        std::vector<TH1D*> hCentFT0M(nRuns);
        std::vector<TH2D*> hPsi2A_CentFT0C(nRuns);
        std::vector<TH2D*> hPsi2B_CentFT0C(nRuns);
        std::vector<TH2D*> hPsi2C_CentFT0C(nRuns);
        std::vector<TH1D*> hPt(nRuns);
        std::vector<TH1D*> hEta(nRuns);
        std::vector<TH1D*> hPhi(nRuns);
        std::vector<TH1D*> hpdca(nRuns);
        std::vector<TH1D*> hRAtAbsorberEnd(nRuns);
        std::vector<TH1D*> hChi2MCHMID(nRuns);
        std::vector<TH1D*> hChi2MCHMFT(nRuns);

        TH2D *hZorroInfo, *hZorroSel; // Histograms for skimmed datasets
        for (auto run : runList) {
            std::cout << "-------- Processing run " << run << " --------" << std::endl;
            TFile *fIn = new TFile(Form("%s/%d/AnalysisResults.root", pathIn.c_str(), run));

            TH1D *hCounterTVX = (TH1D*) fIn->Get("eventselection-run3/luminosity/hCounterTVX");
            TH1D *hCounterTVXafterBCcuts = (TH1D*) fIn->Get("eventselection-run3/luminosity/hCounterTVXafterBCcuts");
            TH1D *hColCounterAll = (TH1D*) fIn->Get("eventselection-run3/eventselection/hColCounterAll");
            TH1D *hColCounterTVX = (TH1D*) fIn->Get("eventselection-run3/eventselection/hColCounterTVX");
            TH1D *hColCounterAcc = (TH1D*) fIn->Get("eventselection-run3/eventselection/hColCounterAcc");

            if (production == "skimmed") {
                TList *hlistZorro = (TList*) fIn->Get("table-maker/Statistics");
                hZorroInfo = (TH2D*) hlistZorro->FindObject("ZorroInfo");
                hZorroSel = (TH2D*) hlistZorro->FindObject("ZorroSel");
            }

            TList *hlistTm = (TList*) fIn->Get("table-maker/output");
            TList *listTmEvents = (TList*) hlistTm->FindObject("Event_AfterCuts");
            hVtxZ[iRun] = (TH1D*) listTmEvents->FindObject("VtxZ");
            hCentFT0C[iRun] = (TH1D*) listTmEvents->FindObject("CentFT0C");
            hCentFT0M[iRun] = (TH1D*) listTmEvents->FindObject("CentFT0M");
            hPsi2A_CentFT0C[iRun] = (TH2D*) listTmEvents->FindObject("Psi2A_CentFT0C");
            hPsi2B_CentFT0C[iRun] = (TH2D*) listTmEvents->FindObject("Psi2B_CentFT0C");
            hPsi2C_CentFT0C[iRun] = (TH2D*) listTmEvents->FindObject("Psi2C_CentFT0C");

            hCountersVsRuns->SetBinContent(iRun+1, 1, hCounterTVX->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 2, hCounterTVXafterBCcuts->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 3, hColCounterAll->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 4, hColCounterTVX->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 5, hColCounterAcc->GetEntries());
            hCountersVsRuns->SetBinContent(iRun+1, 6, hVtxZ[iRun]->GetEntries());
            if (production == "skimmed") {
                hCountersVsRuns->SetBinContent(iRun+1, 7, hZorroInfo -> GetBinContent(1, hZorroInfo -> GetYaxis() -> FindBin("inspected TVX")));
                hCountersVsRuns->SetBinContent(iRun+1, 8, hZorroInfo -> GetBinContent(1, hZorroInfo -> GetYaxis() -> FindBin(Form("%s scalers", triggerMask.c_str()))));
                hCountersVsRuns->SetBinContent(iRun+1, 9, hZorroInfo -> GetBinContent(1, hZorroInfo -> GetYaxis() -> FindBin(Form("%s selections", triggerMask.c_str()))));
                hCountersVsRuns->SetBinContent(iRun+1, 10, hZorroSel -> GetBinContent(1, hZorroSel -> GetYaxis() -> FindBin(Form("%s AnalysedTriggers", triggerMask.c_str()))));
                hCountersVsRuns->SetBinContent(iRun+1, 11, hZorroSel -> GetBinContent(1, hZorroSel -> GetYaxis() -> FindBin(triggerMask.c_str())));
            }

            hCountersVsRuns->GetXaxis()->SetBinLabel(iRun+1, Form("%d", run));

            TList *listTmMuons = (TList*) hlistTm->FindObject(Form("Muons_%s", muonCutTm.c_str()));
            hPt[iRun] = (TH1D*) listTmMuons->FindObject("Pt"); hPt[iRun]->Rebin(10);
            hEta[iRun] = (TH1D*) listTmMuons->FindObject("Eta");
            hPhi[iRun] = (TH1D*) listTmMuons->FindObject("Phi");
            hpdca[iRun] = (TH1D*) listTmMuons->FindObject("pdca");
            hRAtAbsorberEnd[iRun] = (TH1D*) listTmMuons->FindObject("RAtAbsorberEnd");
            hChi2MCHMID[iRun] = (TH1D*) listTmMuons->FindObject("Chi2MCHMID");
            hChi2MCHMFT[iRun] = (TH1D*) listTmMuons->FindObject("Chi2MCHMFT");

            if (iRun == 0) {
                hSumVtxZ = (TH1D*) hVtxZ[iRun]->Clone("hSumVtxZ");
                hSumCentFT0C = (TH1D*) hCentFT0C[iRun]->Clone("hSumCentFT0C");
                hSumCentFT0M = (TH1D*) hCentFT0M[iRun]->Clone("hSumCentFT0M");
                if (production == "heavyIon") {
                    hSumPsi2A_CentFT0C = (TH2D*) hPsi2A_CentFT0C[iRun]->Clone("hSumPsi2A_CentFT0C");
                    hSumPsi2B_CentFT0C = (TH2D*) hPsi2B_CentFT0C[iRun]->Clone("hSumPsi2B_CentFT0C");
                    hSumPsi2C_CentFT0C = (TH2D*) hPsi2C_CentFT0C[iRun]->Clone("hSumPsi2C_CentFT0C");
                }
                hSumPt = (TH1D*) hPt[iRun]->Clone("hSumPt");
                hSumEta = (TH1D*) hEta[iRun]->Clone("hSumEta");
                hSumPhi = (TH1D*) hPhi[iRun]->Clone("hSumPhi");
                hSumpdca = (TH1D*) hpdca[iRun]->Clone("hSumpdca");
                hSumRAtAbsorberEnd = (TH1D*) hRAtAbsorberEnd[iRun]->Clone("hSumRAtAbsorberEnd");
                hSumChi2MCHMID = (TH1D*) hChi2MCHMID[iRun]->Clone("hSumChi2MCHMID");
                hSumChi2MCHMFT = (TH1D*) hChi2MCHMFT[iRun]->Clone("hSumChi2MCHMFT");
            } else {
                hSumVtxZ->Add(hVtxZ[iRun]);
                hSumCentFT0C->Add(hCentFT0C[iRun]);
                hSumCentFT0M->Add(hCentFT0M[iRun]);
                if (production == "heavyIon") {
                    hSumPsi2A_CentFT0C->Add(hPsi2A_CentFT0C[iRun]);
                    hSumPsi2B_CentFT0C->Add(hPsi2B_CentFT0C[iRun]);
                    hSumPsi2C_CentFT0C->Add(hPsi2C_CentFT0C[iRun]);
                }
                hSumPt->Add(hPt[iRun]);
                hSumEta->Add(hEta[iRun]);
                hSumPhi->Add(hPhi[iRun]);
                hSumpdca->Add(hpdca[iRun]);
                hSumRAtAbsorberEnd->Add(hRAtAbsorberEnd[iRun]);
                hSumChi2MCHMID->Add(hChi2MCHMID[iRun]);
                hSumChi2MCHMFT->Add(hChi2MCHMFT[iRun]);
            }

            hVtxZ[iRun]->SetTitle(Form("%d", run));
            hCentFT0C[iRun]->SetTitle(Form("%d", run));
            hCentFT0M[iRun]->SetTitle(Form("%d", run));
            if (production == "heavyIon") {
                hPsi2A_CentFT0C[iRun]->SetTitle(Form("%d", run));
                hPsi2B_CentFT0C[iRun]->SetTitle(Form("%d", run));
                hPsi2C_CentFT0C[iRun]->SetTitle(Form("%d", run));
            }
            hPt[iRun]->SetTitle(Form("%d", run));
            hEta[iRun]->SetTitle(Form("%d", run));
            hPhi[iRun]->SetTitle(Form("%d", run));
            hpdca[iRun]->SetTitle(Form("%d", run));
            hRAtAbsorberEnd[iRun]->SetTitle(Form("%d", run));
            hChi2MCHMID[iRun]->SetTitle(Form("%d", run));
            hChi2MCHMFT[iRun]->SetTitle(Form("%d", run));

            iRun++;
        }

        TCanvas *canvasNormalization = new TCanvas("canvasNormalization", "", 800, 600);
        gStyle->SetOptStat(false);
        hCountersVsRuns->Draw("COLZ TEXT");

        TCanvas *canvasVtxZ = DoRatioPlot(hSumVtxZ, hVtxZ, "canvasVtxZ", "vtx_{Z} (cm)");
        TCanvas *canvasCentFT0C = DoRatioPlot(hSumCentFT0C, hCentFT0C, "canvasCentFT0C", "Centr FT0C (%)");
        TCanvas *canvasCentFT0M = DoRatioPlot(hSumCentFT0M, hCentFT0M, "canvasCentFT0M", "Centr FT0M (%)");

        // Rebin EP distribution
        TCanvas *canvasPsi2A_CentFT0C = nullptr;
        TCanvas *canvasPsi2B_CentFT0C = nullptr;
        TCanvas *canvasPsi2C_CentFT0C = nullptr;
        if (production == "heavyIon") {
            hSumPsi2A_CentFT0C->RebinX(5);
            hSumPsi2B_CentFT0C->RebinX(5);
            hSumPsi2C_CentFT0C->RebinX(5);
            for (int iRun = 0;iRun < nRuns;iRun++) {
                hPsi2A_CentFT0C[iRun]->RebinX(5);
                hPsi2B_CentFT0C[iRun]->RebinX(5);
                hPsi2C_CentFT0C[iRun]->RebinX(5);
            }

            canvasPsi2A_CentFT0C = new TCanvas("canvasPsi2A_CentFT0C", "", 1800, 1800);
            canvasPsi2A_CentFT0C->Divide(3, 3);

            for (int iCentr = 0;iCentr < 9;iCentr++) {
                canvasPsi2A_CentFT0C->cd(iCentr+1);
                TH1D *hTmpSumPsi2A_CentFT0C = (TH1D*) hSumPsi2A_CentFT0C->ProjectionY(Form("hTmpSumPsi2A_CentFT0C_%d_%d", iRun, iCentr), iCentr+1, iCentr+1);
                SetHist(hTmpSumPsi2A_CentFT0C, kBlack, 20);
                hTmpSumPsi2A_CentFT0C->Draw("EP");
                for (int iRun = 0;iRun < nRuns;iRun++) {
                    TH1D *hTmpPsi2A_CentFT0C = (TH1D*) hPsi2A_CentFT0C[iRun]->ProjectionY(Form("hTmpPsi2A_CentFT0C_%d_%d", iRun, iCentr), iCentr+1, iCentr+1);
                    hTmpPsi2A_CentFT0C->Scale(1./hTmpPsi2A_CentFT0C->Integral());
                    hTmpPsi2A_CentFT0C->Draw("SAME PLC PMC");
                }
            }

            canvasPsi2B_CentFT0C = new TCanvas("canvasPsi2B_CentFT0C", "", 1800, 1800);
            canvasPsi2B_CentFT0C->Divide(3, 3);

            for (int iCentr = 0;iCentr < 9;iCentr++) {
                canvasPsi2B_CentFT0C->cd(iCentr+1);
                TH1D *hTmpSumPsi2B_CentFT0C = (TH1D*) hSumPsi2B_CentFT0C->ProjectionY(Form("hTmpSumPsi2B_CentFT0C_%d_%d", iRun, iCentr), iCentr+1, iCentr+1);
                SetHist(hTmpSumPsi2B_CentFT0C, kBlack, 20);
                hTmpSumPsi2B_CentFT0C->Draw("EP");
                for (int iRun = 0;iRun < nRuns;iRun++) {
                    TH1D *hTmpPsi2B_CentFT0C = (TH1D*) hPsi2B_CentFT0C[iRun]->ProjectionY(Form("hTmpPsi2B_CentFT0C_%d_%d", iRun, iCentr), iCentr+1, iCentr+1);
                    hTmpPsi2B_CentFT0C->Scale(1./hTmpPsi2B_CentFT0C->Integral());
                    hTmpPsi2B_CentFT0C->Draw("SAME PLC PMC");
                }
            }

            canvasPsi2C_CentFT0C = new TCanvas("canvasPsi2C_CentFT0C", "", 1800, 1800);
            canvasPsi2C_CentFT0C->Divide(3, 3);

            for (int iCentr = 0;iCentr < 9;iCentr++) {
                canvasPsi2C_CentFT0C->cd(iCentr+1);
                TH1D *hTmpSumPsi2C_CentFT0C = (TH1D*) hSumPsi2C_CentFT0C->ProjectionY(Form("hTmpSumPsi2C_CentFT0C_%d_%d", iRun, iCentr), iCentr+1, iCentr+1);
                SetHist(hTmpSumPsi2C_CentFT0C, kBlack, 20);
                hTmpSumPsi2C_CentFT0C->Draw("EP");
                for (int iRun = 0;iRun < nRuns;iRun++) {
                    TH1D *hTmpPsi2C_CentFT0C = (TH1D*) hPsi2C_CentFT0C[iRun]->ProjectionY(Form("hTmpPsi2C_CentFT0C_%d_%d", iRun, iCentr), iCentr+1, iCentr+1);
                    hTmpPsi2C_CentFT0C->Scale(1./hTmpPsi2C_CentFT0C->Integral());
                    hTmpPsi2C_CentFT0C->Draw("SAME PLC PMC");
                }
            }
        }

        TCanvas *canvasPt = DoRatioPlot(hSumPt, hPt, "canvasPt", "#it{p}_{T} (GeV/#it{c})");
        TCanvas *canvasEta = DoRatioPlot(hSumEta, hEta, "canvasEta", "#eta");
        TCanvas *canvasPhi = DoRatioPlot(hSumPhi, hPhi, "canvasPhi", "#varphi");
        TCanvas *canvaspdca = DoRatioPlot(hSumpdca, hpdca, "canvaspdca", "p x DCA");
        TCanvas *canvasRAtAbsorberEnd = DoRatioPlot(hSumRAtAbsorberEnd, hRAtAbsorberEnd, "canvasRAtAbsorberEnd", "R_{abs}");
        TCanvas *canvasChi2MCHMID = DoRatioPlot(hSumChi2MCHMID, hChi2MCHMID, "canvasChi2MCHMID", "#chi^{2}_{MCH-MID}");
        TCanvas *canvasChi2MCHMFT = DoRatioPlot(hSumChi2MCHMFT, hChi2MCHMFT, "canvasChi2MCHMFT", "#chi^{2}_{MFT-MCH}");

        canvasVtxZ->SaveAs(Form("%s/figures/VtxZ.pdf", pathIn.c_str()));
        canvasCentFT0C->SaveAs(Form("%s/figures/CentFT0C.pdf", pathIn.c_str()));
        canvasCentFT0M->SaveAs(Form("%s/figures/CentFT0M.pdf", pathIn.c_str()));
        if (production == "heavyIon") {canvasPsi2A_CentFT0C->SaveAs(Form("%s/figures/Psi2A_CentFT0C.pdf", pathIn.c_str()));}
        canvasPt->SaveAs(Form("%s/figures/Pt.pdf", pathIn.c_str()));
        canvasEta->SaveAs(Form("%s/figures/Eta.pdf", pathIn.c_str()));
        canvasPhi->SaveAs(Form("%s/figures/Phi.pdf", pathIn.c_str()));
        canvaspdca->SaveAs(Form("%s/figures/pdca.pdf", pathIn.c_str()));
        canvasRAtAbsorberEnd->SaveAs(Form("%s/figures/RAtAbsorberEnd.pdf", pathIn.c_str()));
        canvasChi2MCHMID->SaveAs(Form("%s/figures/Chi2MCHMID.pdf", pathIn.c_str()));
        canvasChi2MCHMFT->SaveAs(Form("%s/figures/Chi2MCHMFT.pdf", pathIn.c_str()));


        TFile *fOut = new TFile(Form("%s/streamingResults.root", pathIn.c_str()), "RECREATE");
        hCountersVsRuns->Write();
        canvasVtxZ->Write();
        canvasCentFT0C->Write();
        canvasCentFT0M->Write();
        if (production == "heavyIon") {
            canvasPsi2A_CentFT0C->Write();
            canvasPsi2B_CentFT0C->Write();
            canvasPsi2C_CentFT0C->Write();
        }
        canvasPt->Write();
        canvasEta->Write();
        canvasPhi->Write();
        canvaspdca->Write();
        canvasRAtAbsorberEnd->Write();
        canvasChi2MCHMID->Write();
        canvasChi2MCHMFT->Write();
        fOut->Close();
    } else if (process == "table-reader") {
        TFile *fIn = new TFile(Form("%s/AnalysisResults.root", pathIn.c_str()));

        TList *hlistTrEvents = (TList*) fIn->Get("analysis-event-selection/output");
        TList *listTrEventsBeforeCuts = (TList*) hlistTrEvents->FindObject("Event_BeforeCuts");
        TH1D *hSumBeforeCutsVtxZ = (TH1D*) listTrEventsBeforeCuts->FindObject("VtxZ");

        TList *listTrEvents = (TList*) hlistTrEvents->FindObject("Event_AfterCuts");
        hSumVtxZ = (TH1D*) listTrEvents->FindObject("VtxZ");
        hSumCentFT0C = (TH1D*) listTrEvents->FindObject("CentFT0C");
        hSumCentFT0M = (TH1D*) listTrEvents->FindObject("CentFT0M");
        hSumPsi2A_CentFT0C = (TH2D*) listTrEvents->FindObject("Psi2A_CentFT0C");
        hSumPsi2B_CentFT0C = (TH2D*) listTrEvents->FindObject("Psi2B_CentFT0C");
        hSumPsi2C_CentFT0C = (TH2D*) listTrEvents->FindObject("Psi2C_CentFT0C");

        TList *hlistTrMuons = (TList*) fIn->Get("analysis-muon-selection/output");
        TList *listTrMuons = (TList*) hlistTrMuons->FindObject(Form("TrackMuon_%s", muonCutTr.c_str()));
        hSumPt = (TH1D*) listTrMuons->FindObject("Pt"); hSumPt->Rebin(10);
        hSumEta = (TH1D*) listTrMuons->FindObject("Eta");
        hSumPhi = (TH1D*) listTrMuons->FindObject("Phi");
        hSumpdca = (TH1D*) listTrMuons->FindObject("pdca");
        hSumRAtAbsorberEnd = (TH1D*) listTrMuons->FindObject("RAtAbsorberEnd");
        hSumChi2MCHMID = (TH1D*) listTrMuons->FindObject("Chi2MCHMID");
        hSumChi2MCHMFT = (TH1D*) listTrMuons->FindObject("Chi2MCHMFT");

        TList *hlistTrDimuons = (TList*) fIn->Get("analysis-same-event-pairing/output");
        TList *listTrDimuons = (TList*) hlistTrDimuons->FindObject(Form("PairsMuonSEPM_%s", muonCutTr.c_str()));
        THnSparseD *hSparseDimuon = (THnSparseD*) listTrDimuons -> FindObject(histBlock.c_str());
        hSparseDimuon -> GetAxis(1) -> SetRangeUser(0., 20.);
        hSparseDimuon -> GetAxis(3) -> SetRangeUser(2.5, 4);
        TH1D *hDimuonMass = (TH1D*) hSparseDimuon -> Projection(0, "hDimuonMass");
        TH1D *hCounters = new TH1D("hCounters", ";Run;Counter", 2, 0, 2);
        hCounters->GetXaxis()->SetBinLabel(1, "Collisions before cuts");
        hCounters->GetXaxis()->SetBinLabel(2, "Collisions after cuts");
        hCounters->SetBinContent(1, hSumBeforeCutsVtxZ->GetEntries());
        hCounters->SetBinContent(2, hSumVtxZ->GetEntries());

        TCanvas *canvasVtxZ = DoPlot(hSumVtxZ, "canvasVtxZ", "vtx_{Z} (cm)");
        TCanvas *canvasCentFT0C = DoPlot(hSumCentFT0C, "canvasCentFT0C", "Centr FT0C (%)");
        TCanvas *canvasCentFT0M = DoPlot(hSumCentFT0M, "canvasCentFT0M", "Centr FT0M (%)");


        // Rebin EP distribution
        int minCentrs[] = {0, 10, 20, 30, 40, 50, 60, 70, 80};
        int maxCentrs[] = {10, 20, 30, 40, 50, 60, 70, 80, 90};

        hSumPsi2A_CentFT0C->RebinX(5);
        hSumPsi2B_CentFT0C->RebinX(5);
        hSumPsi2C_CentFT0C->RebinX(5);

        std::vector<TH1D*> hPsi2As_CentFT0C(9);
        std::vector<TH1D*> hPsi2Bs_CentFT0C(9);
        std::vector<TH1D*> hPsi2Cs_CentFT0C(9);

        TCanvas *canvasPsi2A_CentFT0C = new TCanvas("canvasPsi2A_CentFT0C", "", 1800, 1800);
        canvasPsi2A_CentFT0C->Divide(3, 3);

        for (int iCentr = 0;iCentr < 9;iCentr++) {
            canvasPsi2A_CentFT0C->cd(iCentr+1);
            TH1D *hTmpSumPsi2A_CentFT0C = (TH1D*) hSumPsi2A_CentFT0C->ProjectionY(Form("hTmpSumPsi2A_CentFT0C_%d", iCentr), iCentr+1, iCentr+1);
            SetHist(hTmpSumPsi2A_CentFT0C, kBlack, 20);
            hPsi2As_CentFT0C[iCentr] = (TH1D*) hTmpSumPsi2A_CentFT0C->Clone(Form("hPsi2As_CentFT0C_%d_%d", minCentrs[iCentr], maxCentrs[iCentr]));
            SetHist(hPsi2As_CentFT0C[iCentr], kRed+1, 20);
            hTmpSumPsi2A_CentFT0C->Draw("EP");
        }

        TCanvas *canvasPsi2B_CentFT0C = new TCanvas("canvasPsi2B_CentFT0C", "", 1800, 1800);
        canvasPsi2B_CentFT0C->Divide(3, 3);

        for (int iCentr = 0;iCentr < 9;iCentr++) {
            canvasPsi2B_CentFT0C->cd(iCentr+1);
            TH1D *hTmpSumPsi2B_CentFT0C = (TH1D*) hSumPsi2B_CentFT0C->ProjectionY(Form("hTmpSumPsi2B_CentFT0C_%d", iCentr), iCentr+1, iCentr+1);
            SetHist(hTmpSumPsi2B_CentFT0C, kBlack, 20);
            hPsi2Bs_CentFT0C[iCentr] = (TH1D*) hTmpSumPsi2B_CentFT0C->Clone(Form("hPsi2Bs_CentFT0C_%d_%d", minCentrs[iCentr], maxCentrs[iCentr]));
            SetHist(hPsi2Bs_CentFT0C[iCentr], kBlue+1, 24, 0.5);
            hTmpSumPsi2B_CentFT0C->Draw("EP");
        }

        TCanvas *canvasPsi2C_CentFT0C = new TCanvas("canvasPsi2C_CentFT0C", "", 1800, 1800);
        canvasPsi2C_CentFT0C->Divide(3, 3);

        for (int iCentr = 0;iCentr < 9;iCentr++) {
            canvasPsi2C_CentFT0C->cd(iCentr+1);
            TH1D *hTmpSumPsi2C_CentFT0C = (TH1D*) hSumPsi2C_CentFT0C->ProjectionY(Form("hTmpSumPsi2C_CentFT0C_%d", iCentr), iCentr+1, iCentr+1);
            SetHist(hTmpSumPsi2C_CentFT0C, kBlack, 20);
            hPsi2Cs_CentFT0C[iCentr] = (TH1D*) hTmpSumPsi2C_CentFT0C->Clone(Form("hPsi2Cs_CentFT0C_%d_%d", minCentrs[iCentr], maxCentrs[iCentr]));
            SetHist(hPsi2Cs_CentFT0C[iCentr], kGreen+1, 24, 0.5);
            hTmpSumPsi2C_CentFT0C->Draw("EP");
        }

        // Produced comparison plot for different detectors
        TLatex *latexTitle = new TLatex();
        latexTitle->SetTextSize(0.045);
        latexTitle->SetNDC();
        latexTitle->SetTextFont(42);

        TLegend *legend = new TLegend(0.45,0.30,0.65,0.50);
        SetLegend(legend);
        legend->AddEntry(hPsi2As_CentFT0C[0],"TPC","PL");
        legend->AddEntry(hPsi2Bs_CentFT0C[0],"FT0C","PL");
        legend->AddEntry(hPsi2Cs_CentFT0C[0],"FT0A","PL");
        legend->Draw();

        TCanvas *canvasPsi2All_CentFT0C = new TCanvas("canvasPsi2All_CentFT0C", "", 1800, 1800);
        canvasPsi2All_CentFT0C->Divide(3, 3);
        for (int iCentr = 0;iCentr < 9;iCentr++) {
            gStyle->SetOptStat(false);
            canvasPsi2All_CentFT0C->cd(iCentr+1);
            hPsi2As_CentFT0C[iCentr]->GetXaxis()->SetTitle("#Psi_{2}");
            hPsi2As_CentFT0C[iCentr]->GetYaxis()->SetRangeUser(0, 0.020);
            hPsi2As_CentFT0C[iCentr]->Draw("EP");
            hPsi2Bs_CentFT0C[iCentr]->Draw("EP SAME");
            hPsi2Cs_CentFT0C[iCentr]->Draw("EP SAME");
            legend->Draw("SAME");
            latexTitle->DrawLatex(0.65, 0.85, Form("FTOC centrality %d#minus%d%%", minCentrs[iCentr], maxCentrs[iCentr]));
        }
        canvasPsi2All_CentFT0C->SaveAs(Form("%s/figures/Psi2ABC_CentFT0C.pdf", pathIn.c_str()));

        TCanvas *canvasPt = DoPlot(hSumPt, "canvasPt", "#it{p}_{T} (GeV/#it{c})");
        TCanvas *canvasEta = DoPlot(hSumEta, "canvasEta", "#eta");
        TCanvas *canvasPhi = DoPlot(hSumPhi, "canvasPhi", "#varphi");
        TCanvas *canvaspdca = DoPlot(hSumpdca, "canvaspdca", "p x DCA");
        TCanvas *canvasRAtAbsorberEnd = DoPlot(hSumRAtAbsorberEnd, "canvasRAtAbsorberEnd", "R_{abs}");
        TCanvas *canvasChi2MCHMID = DoPlot(hSumChi2MCHMID, "canvasChi2MCHMID", "#chi^{2}_{MCH-MID}");

        TFile *fOut = new TFile(Form("%s/QA.root", pathIn.c_str()), "RECREATE");
        hCounters->Write();
        hDimuonMass->Write("hDimuonMass");
        canvasVtxZ->Write();
        canvasCentFT0C->Write();
        canvasCentFT0M->Write();
        canvasPsi2A_CentFT0C->Write();
        canvasPsi2B_CentFT0C->Write();
        canvasPsi2C_CentFT0C->Write();
        canvasPt->Write();
        canvasEta->Write();
        canvasPhi->Write();
        canvaspdca->Write();
        canvasRAtAbsorberEnd->Write();
        canvasChi2MCHMID->Write();
        fOut->Close();

    } else {
        std::cout << "-------- No process selected --------" << std::endl;
        return;
    }
}