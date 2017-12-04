
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


void effAllruns(){

  gStyle->SetErrorX(0);
  
TFile* file = new TFile("./LabData.root");
TTree* tree;
file->GetObject("T",tree);
  
TCanvas* c = new TCanvas("c");
TCanvas* c0 = new TCanvas("c0");
vector<int> runs={2,3,4,5,6,7,8};
vector<int> mCol={4,4,4,4,4,4,4};
vector<int> mStyle={20,21,22,23,24,25,26};


TH1D* hNum = new TH1D("hNum","eff. vs  #p.e.; #p.e.; eff.",200,0,600);
TH1D* hDeNum = new TH1D("hDeNum","eff. vs  #p.e.; #p.e.; eff.",200,0,600);

hDeNum->SetStats(0);
//hDeNum->GetYaxis()->SetRangeUser(0.5,1.01);

hDeNum->SetTitle("eff. vs  #p.e.; #p.e.; eff.");
c->cd();
hDeNum->DrawCopy();

TEfficiency* eff;

TLegend* leg = new TLegend(0.7,0.4,0.85,0.88);

 for(int i=0;i<(int)runs.size();i++){
  TString cut(""); cut.Form("isTrig==1&&runNr==%d",runs.at(i));
  TH1D hIntegral("hIntegral","; #p.e.;",200,0,600);
  c0->cd();
  tree->Draw("Integral[6]>>hIntegral",cut);
  cout << "runNr: " << runs.at(i) << endl;
  for(int bin=1;bin<=200;bin++){
    hNum->SetBinContent(bin,hIntegral.Integral(bin,200));
    hDeNum->SetBinContent(bin,hIntegral.Integral(1,200));
    //cout << hIntegral.Integral(bin,200)/hIntegral.Integral(1,200) << endl;
    
  }
      eff=new TEfficiency(*hNum,*hDeNum);
      eff->SetMarkerColor(mCol.at(i));
      eff->SetMarkerStyle(mStyle.at(i));
      c->cd();
      eff->DrawClone("same,e1");
      
      TString legText(""); legText.Form("run %d",runs.at(i));
      leg->AddEntry(eff,legText,"p");
 }
  
  leg->Draw("same");
  c->Print("./effAllruns.pdf","pdf");
  
}
 
// // //   h1->SetLineColor(color[0]);
// // //   leg1->AddEntry(h1,"PMT-1 1000 V","l");



