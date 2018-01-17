
//root
#include <TH1D.h>
//#include <TStyle.h>

// #include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2D.h>		
#include <TEfficiency.h>


//C, C++
#include <vector>


void createWeight(){

  gStyle->SetErrorX(0);
  
  
TFile* file_WOM = new TFile("./runs/7/out.root");
TFile* file_noWOM = new TFile("./runs/10/out.root");
TTree* tree_WOM; 
TTree* tree_noWOM;

file_WOM->GetObject("T",tree_WOM);
file_noWOM->GetObject("T",tree_noWOM);
  
TCanvas* c = new TCanvas("c");
c->Divide(2,1);

Int_t nbins = 100;

TH1D* h_WOM = new TH1D("h_WOM","WOM spectrum",nbins,0,1000);
  h_WOM->SetLineColor(2);
  h_WOM->SetLineWidth(2);
TH1D* h_noWOM = new TH1D("h_noWOM","no WOM spectrum",nbins,0,1000);
  h_noWOM->SetLineWidth(2);
TH1D* h_factor = new TH1D("h_factor","weight factor WOM/noWOM",nbins,0,1000);
  h_factor->SetMarkerStyle(20);

      c->cd(1);
      tree_noWOM->Draw("Integral_0_300[6]>>h_noWOM","isTrig&&runNr==10");
      tree_WOM->Draw("Integral_0_300[6]>>h_WOM","isTrig&&runNr==7","same");
      
      TLegend* leg = new TLegend(0.7,0.4,0.85,0.88);
      leg->AddEntry(h_WOM,"WOM","l");
      leg->AddEntry(h_noWOM,"no WOM","l");
      
      h_noWOM->Scale(1./h_noWOM->GetEntries());
      h_WOM->Scale(1./h_WOM->GetEntries());
      h_noWOM->Draw("hist");
      h_WOM->Draw("histsame");
      leg->Draw("same");
      c->cd(2);
      
      h_factor->Divide(h_WOM,h_noWOM);
      h_factor->Draw("");
      
  c->Print("./createWeight.pdf","pdf");
  
  
  TFile f("weight.root","RECREATE");
  h_factor->Write();
  f.Close();
}


