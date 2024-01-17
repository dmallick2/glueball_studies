#include <iostream>
#include <cmath>
#include "src/style.h"
#include "src/fitfunc.h"

using namespace std;

void compare2()
{
    TLatex *t2 = new TLatex();
    t2->SetNDC(); // to self adjust the text so that it remains in the box
    t2->SetTextSize(0.04);
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);
    gStyle->SetOptFit(0);

    TFile *file[3];
    TH1F *hmass[3];
    TH1F *hwidth[3];
    TH1F *hyield[3];

    file[0] = new TFile("/home/sawan/check_k892/output/pp/kstar/LHC23_pass1_lowB_lowIR/LHC23_pass1_lowB_lowIR.root", "READ");
    file[1] = new TFile("/home/sawan/check_k892/output/pp/kstar/LHC23_pass1_lowB_highIR_sampling/LHC23_pass1_lowB_highIR_sampling.root", "READ");
    file[2] = new TFile("/home/sawan/check_k892/output/pbpb/kstar/LHC23zzh_pass1/fixed_width/pbpb.root", "READ");

    if (file[0]->IsZombie() || file[1]->IsZombie() || file[2]->IsZombie())
    {
        cerr << "File not found" << endl;
        return;
    }

    TCanvas *cmass = new TCanvas("", "", 1200, 1000);
    SetCanvasStyle2(cmass, 0.15, 0.05, 0.05, 0.15);
    int color[3] = {1, 2, 4};
    int marker[3] = {53,55,54};

    // mass vs pT

    for (int i = 0; i < 3; i++)
    {
        hmass[i] = (TH1F *)file[i]->Get("mass");
        hwidth[i] = (TH1F *)file[i]->Get("width");
        hyield[i] = (TH1F *)file[i]->Get("yield");
        hmass[i]->SetMarkerStyle(marker[i]);
        hmass[i]->SetMarkerColor(color[i]);
        hmass[i]->SetLineColor(color[i]);
        hmass[i]->SetLineWidth(2);
        hmass[i]->SetMarkerSize(2);
        hwidth[i]->SetMarkerStyle(marker[i]);
        hwidth[i]->SetMarkerColor(color[i]);
        hwidth[i]->SetLineColor(color[i]);
        hwidth[i]->SetLineWidth(2);
        hwidth[i]->SetMarkerSize(2);
        hyield[i]->SetMarkerStyle(marker[i]);
        hyield[i]->SetMarkerColor(color[i]);
        hyield[i]->SetLineColor(color[i]);
        hyield[i]->SetMarkerSize(2);
        hyield[i]->SetLineWidth(2);
    }

    TH1F *dummy = new TH1F("", "", 16 * 5, -1, 16);
    SetHistoStyle(dummy, 1, 20, 1, 0.05, 0.05, 0.045, 0.045, 1.13, 1.8);
    dummy->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    dummy->GetYaxis()->SetTitle("Mass (GeV/c^{2})");
    dummy->GetYaxis()->SetTitleOffset(1.4);

    dummy->SetMaximum(0.945);
    dummy->SetMinimum(0.855);
    dummy->GetXaxis()->SetRangeUser(-0.01, 16.1);
    dummy->Draw();
    hmass[1]->Draw("same");
    hmass[0]->Draw("same");
    hmass[2]->Draw("same");
    TLegend *massleg = new TLegend(0.45, 0.68, 0.9, 0.91);
    massleg->SetFillStyle(0);
    massleg->SetBorderSize(0);
    massleg->SetTextFont(42);
    massleg->SetTextSize(0.04);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    // TLine *line = new TLine(hmass[1]->GetXaxis()->GetXmin(), 0.895, hmass[1]->GetXaxis()->GetXmax(), 0.895);
    TLine *line = new TLine(-0.01, 0.895, 16.1, 0.895);
    SetLineStyle(line, kGreen+2);
    line->SetLineWidth(3);
    massleg->AddEntry(line, "PDG Mass (0.895 GeV/c^{2})", "l");
    massleg->AddEntry(hmass[0], "pp 13.6 TeV low IR", "lep");
    massleg->AddEntry(hmass[1], "pp 13.6 TeV high IR", "lep");
    massleg->AddEntry(hmass[2], "Pb-Pb 5.36 TeV", "lep");
    massleg->Draw("l");
    cmass->SaveAs("/home/sawan/check_k892/output/compare/mass.png");

    TCanvas *cwidth = new TCanvas("", "", 1200, 1000);
    SetCanvasStyle2(cwidth, 0.15, 0.05, 0.05, 0.15);

    // // width vs pT
    TH1F *dummy2 = new TH1F("", "", 17 * 5, -1, 16);
    SetHistoStyle(dummy2, 1, 20, 1, 0.05, 0.05, 0.045, 0.045, 1.13, 1.8);
    dummy2->GetXaxis()->SetTitle("p_{T} (GeV/c)");
    dummy2->GetYaxis()->SetTitle("Width (GeV)");
    dummy2->GetYaxis()->SetTitleOffset(1.4);

    dummy2->SetMaximum(0.21);
    dummy2->SetMinimum(0);
    dummy2->GetXaxis()->SetRangeUser(-0.01, 16.1);
    dummy2->Draw();
    hwidth[1]->Draw("same");
    hwidth[0]->Draw("same");
    // hwidth[2]->Draw("same");
    TLegend *widthleg = new TLegend(0.48, 0.74, 0.9, 0.91);
    widthleg->SetFillStyle(0);
    widthleg->SetBorderSize(0);
    widthleg->SetTextFont(42);
    widthleg->SetTextSize(0.04);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
    // TLine *line = new TLine(hmass[1]->GetXaxis()->GetXmin(), 0.895, hmass[1]->GetXaxis()->GetXmax(), 0.895);
    TLine *line2 = new TLine(0, 0.047, 16, 0.047);
    SetLineStyle(line2, 4);
    widthleg->AddEntry(line2, "PDG Width (0.047 GeV)", "l");
    widthleg->AddEntry(hmass[0], "pp 13.6 TeV low IR", "lep");
    widthleg->AddEntry(hmass[1], "pp 13.6 TeV high IR", "lep");
    // widthleg->AddEntry(hmass[2], "Pb-Pb 5.36 TeV", "lep");
    widthleg->Draw("l");
    cwidth->SaveAs("/home/sawan/check_k892/output/compare/width.png");


    //yield plots
      TCanvas *cyield = new TCanvas("", "", 1200, 1000);
    SetCanvasStyle2(cyield, 0.15, 0.05, 0.05, 0.15);
    gPad->SetLogy();
     TH1F *dummy3 = new TH1F("", "", 17 * 5, -1, 16);
    SetHistoStyle(dummy3, 1, 20, 1, 0.05, 0.05, 0.045, 0.045, 1.13, 1.8);
     dummy3->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    dummy3->GetYaxis()->SetTitle("1/#it{N}_{Ev}d^{2}#it{N}/(d#it{y}d#it{p}_{T}) [(GeV/#it{c})^{-1}]");
    dummy3->GetYaxis()->SetTitleOffset(1.4);
    dummy3->SetMaximum(5);
    dummy3->SetMinimum(15e-8);
    dummy3->GetXaxis()->SetRangeUser(-0.01, 16.1);
    dummy3->Draw();
    hyield[1]->Draw("same");
    hyield[0]->Draw("same");
    hyield[2]->Draw("same");
    TLegend *yieldleg = new TLegend(0.48, 0.74, 0.9, 0.91);
    yieldleg->SetFillStyle(0);
    yieldleg->SetBorderSize(0);
    yieldleg->SetTextFont(42);
    yieldleg->SetTextSize(0.04);
    t2->DrawLatex(0.28, 0.96, "#bf{K(892)^{0} #rightarrow #pi + K}");
   
    yieldleg->AddEntry(hmass[0], "pp 13.6 TeV low IR", "lep");
    yieldleg->AddEntry(hmass[1], "pp 13.6 TeV high IR", "lep");
    yieldleg->AddEntry(hmass[2], "Pb-Pb 5.36 TeV", "lep");
    yieldleg->Draw("l");
    cyield->SaveAs("/home/sawan/check_k892/output/compare/yield.png");


    // ratio of yields
     TCanvas *cyield_ratio = new TCanvas("", "", 1200, 1000);
    SetCanvasStyle2(cyield_ratio, 0.15, 0.05, 0.05, 0.15);
     TH1F *dummy4 = new TH1F("", "", 17 * 5, -1, 16);
     TH1F *hratio = new TH1F ("", "", hyield[0]->GetNbinsX(), hyield[0]->GetXaxis()->GetXmin(), hyield[0]->GetXaxis()->GetXmax());
    SetHistoStyle(dummy4, 1, 20, 1, 0.05, 0.05, 0.045, 0.045, 1.13, 1.8);
     dummy4->GetXaxis()->SetTitle("#it{p}_{T} (GeV/#it{c})");
    dummy4->GetYaxis()->SetTitle("High IR/ Low IR");
    dummy4->GetYaxis()->SetTitleOffset(1.4);
    dummy4->SetMaximum(1.8);
    dummy4->SetMinimum(0.5);
    dummy4->GetXaxis()->SetRangeUser(-0.01, 9);
    for (int i = 0; i < hyield[0]->GetNbinsX(); i++)
    {
        float ratio = hyield[1]->GetBinContent(i+1)/hyield[0]->GetBinContent(i+1);
        float ratio_err = ratio*sqrt(pow(hyield[0]->GetBinError(i+1)/hyield[0]->GetBinContent(i+1),2)+pow(hyield[1]->GetBinError(i+1)/hyield[1]->GetBinContent(i+1),2));
        cout<<"ratio "<<ratio<<endl;
        hratio->SetBinContent(i+1,ratio);
        hratio->SetBinError(i+1,ratio_err);
    }
    dummy4->Draw();
    hratio->SetMarkerStyle(53);
    hratio->SetMarkerColor(1);
    hratio->SetLineColor(1);
    hratio->SetLineWidth(2);
    hratio->SetMarkerSize(2);
    hratio->Draw("same");
    TLatex lat;
    lat.SetNDC();
    lat.SetTextSize(0.06);
    lat.SetTextFont(42);
    lat.DrawLatex(0.5,0.88,"pp 13.6 TeV");
    cyield_ratio->SaveAs("/home/sawan/check_k892/output/compare/yield_ratio.png");
    
   
}