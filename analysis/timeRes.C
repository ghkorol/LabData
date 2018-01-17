
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
#include <TTreeReader.h>


//C, C++
#include <vector>


void timeRes(){

  gStyle->SetErrorX(0);
  
TFile* file = new TFile("./LabData.root");
TTree* tree= nullptr;
TTreeReader theReader("T",file);
TTreeReaderValue<Float_t> tPMT(theReader, "tPMT");
TTreeReaderValue<Float_t> w(theReader, "w");
TTreeReaderValue<Int_t> isTrig(theReader, "isTrig");
TTreeReaderValue<Int_t> runNr(theReader, "runNr");
TTreeReaderArray<Float_t> Integral_0_300(theReader, "Integral_0_300");

TCanvas* c = new TCanvas("c");

//isTrig&&runNr==10&&Integral_0_300[6]>20&&Integral_0_300[6]<60
TH1D* hTimeRes = new TH1D("hTimeRes",";tPMT, ns;",100,-40,0);


   while(theReader.Next()){
      if (!(*isTrig && *runNr==10 && Integral_0_300[6]>0 && Integral_0_300[6]<100)) {
         continue; 
      }
      hTimeRes->Fill(*tPMT,*w) ;
     
  }

  hTimeRes->Draw();
  gStyle->SetOptFit(1);
  hTimeRes->Fit("gaus");
  c->Print("./timeRes.pdf","pdf");
  
}


