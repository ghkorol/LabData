//root
#include <TLine.h>
//#include <TStyle.h>
#include <TVirtualFFT.h>
#include <TString.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TTree.h>
#include <TH1D.h>
#include <TMath.h>
#include <TF1.h>

//C, C++
#include <stdio.h>
#include <assert.h>
//#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <vector>
//#include <string>
//#include <iomanip>

using namespace std;

const int n_Ch = 7;

float SP = 0.3125;
float pe = 47.46;//mV*ns
int wavesPrintRate = 1000000;
int ftPrintRate = 1000000;
int trigPrintRate = 1000000;//100
int signalPrintRate = 100000;//100
double coef = 2.5 / (4096 * 10);

int runNr = -999;

void cleanEventMemory(std::vector<TObject*>& trash);
float CDF(TH1F* hWave,TF1* fTrigFit,float thr);
float CDFinvert(TH1F* hWave,float thr);
float iCFD(TH1F* hWave,float t,float thrNpe,float BL);
float integral(TH1F* hWave,float t1,float t2,float BL);
 
float* getBL(TH1F* hWave, float* BL, float t1, float t2); 
TString getRunName(TString inDataFolder);
void read(TString inFileList, TString inDataFolder, TString outFile);

int main(int argc, char *argv[]){
  TString inFileList;
  TString inDataFolder;
  TString outFile;

  if(argc == 5){
    inFileList = argv[1];
    inDataFolder = argv[2];
    outFile = argv[3];
    runNr=atoi(argv[4]);
    cout<<"In data file list : "<<inFileList<<endl
        <<"In data path      : "<<inDataFolder<<endl
        <<"Out root file     : "<<outFile<<endl;
    read(inFileList, inDataFolder, outFile);
  }
  else{
    cout<<" ERROR --->  in input arguments "<<endl
        <<"        [1] - in data file list"<<endl
        <<"        [2] - in data path"<<endl
        <<"        [3] - out root file"<<endl;
  }
  return 0;
}

void read(TString _inFileList, TString _inDataFolder, TString _outFile){

  TH1D* h_weight = nullptr;
  TFile f("weight.root");
   f.GetObject("h_factor",h_weight);
   h_weight->SetDirectory(0);
  f.Close();

  
  TF1* fTrigFit = new TF1("fTrigFit","gaus");
  fTrigFit->SetParameter(0,800);
  fTrigFit->SetParameter(2,1);
  fTrigFit->SetLineWidth(1);
  
  ///////////////////Root file with data/////////////////
  TFile *rootFile = new TFile( _outFile, "RECREATE");
  if (rootFile->IsZombie()) {
    cout << "PROBLEM with the initialization of the output ROOT ntuple "
	 << _outFile << ": check that the path is correct!!!"
	 << endl;
    exit(-1);
  }
  TTree *tree = new TTree("T", "USBWC Data Tree");
  //rootFile->SetCompressionLevel(2);
  //tree->SetAutoSave(1000000);
  // Create new event
  TTree::SetBranchStyle(0);

  Int_t EventNumber=-999;
  Float_t SamplingPeriod = -999;
  Double_t EpochTime = -999;
  Int_t Year = -999;
  Int_t Month = -999;
  Int_t Day = -999;
  Int_t Hour = -999;
  Int_t Minute = -999;
  Int_t Second = -999;
  Int_t Millisecond = -999;
  Float_t trigT = -999;//t_trig = (t0+t1+t2+t3)/4
  Float_t tPMT = -999;
  Float_t w = -999;
  Float_t tPMTi = -999;
  Float_t trigTp = -999;//t_trig' = [(t0+t1)-(t2+t3)]/4
  Int_t isTrig = -999;
  Float_t trigGate = -999;
  Int_t nCh = -1;
  int nActiveCh = -1;
  Int_t ChannelNr[n_Ch];
  std::vector<float> amp(n_Ch,-999);
  std::vector<float> max(n_Ch,-999);
  std::vector<float> min(n_Ch,-999);
  Float_t t[n_Ch];
  Float_t BL[n_Ch];//store baseline for n_Ch channels
  Float_t BL_RMS[n_Ch];//store rms of baseline for n_Ch channels
  float BL_output[2];//array used for output getBL-function
  float Integral_0_300[n_Ch];//array used to store Integral of signal from 0 to 300ns
  float Integral[n_Ch];
  int NumberOfBins;
  Int_t EventIDsamIndex[n_Ch];
  Int_t FirstCellToPlotsamIndex[n_Ch];
  
  std::vector<TH1F*> hChSum;
  for(int i=0;i<n_Ch;i++){
    TString name("");
    name.Form("hChSum_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChSum.push_back(h);
  }
  
  std::vector<TH1F*> hChShift;
  for(int i=0;i<n_Ch;i++){
    TString name("");
    name.Form("hChShift_%d",i);
    TH1F* h = new TH1F("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h->SetName(name);
    hChShift.push_back(h);
  }
  
  
  std::vector<TH1F> hChtemp;
  for(int i=0;i<n_Ch;i++){
    TString name("");
    name.Form("hChtemp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChtemp.push_back(h);
  }
  
  std::vector<TH1F> hChShift_temp;
  for(int i=0;i<n_Ch;i++){
    TString name("");
    name.Form("hChShift_temp_%d",i);
    TH1F h("h",";ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
    h.SetName(name);
    hChShift_temp.push_back(h);
  }
  
  
  Short_t amplValues[n_Ch][1024];
  TH1F hCh("hCh","dummy;ns;Amplitude, mV",1024,-0.5*SP,1023.5*SP);
  TString plotSaveFolder  = _outFile;
  plotSaveFolder.ReplaceAll("out.root","");
  TCanvas cWaves("cWaves","cWaves",1000,700);
  cWaves.Divide(3,3);
  TCanvas cCh0("cCh0","cCh0",1500,900);
  cCh0.Divide(2,2);
  TCanvas cTrig("cTrig","cTrig",1500,900);
  cTrig.Divide(2,2);
  TCanvas cSignal("cSignal","cSignal",1500,900);
  cSignal.Divide(2,2);

   //Event USBWC
  tree->Branch("EventNumber",&EventNumber, "EventNumber/I");
  tree->Branch("SamplingPeriod", &SamplingPeriod,  "SamplingPeriod/F");
  tree->Branch("EpochTime",&EpochTime, "EpochTime/D");
  tree->Branch("Year",&Year, "Year/I");
  tree->Branch("Month",&Month, "Month/I");
  tree->Branch("Day",&Day, "Day/I");
  tree->Branch("Hour",&Hour, "Hour/I");
  tree->Branch("Minute",&Minute, "Minute/I");
  tree->Branch("Second",&Second, "Second_/I");
  tree->Branch("Millisecond",&Millisecond, "Millisecond/I");
  tree->Branch("trigT",&trigT, "trigT/F");
  tree->Branch("tPMT",&tPMT, "tPMT/F");
  tree->Branch("w",&w, "w/F");
  tree->Branch("tPMTi",&tPMTi, "tPMTi/F");

  tree->Branch("trigGate",&trigGate,"trigGate/F");
  
  tree->Branch("trigTp",&trigTp, "trigTp/F");
  tree->Branch("isTrig",&isTrig,"isTrig/I");
  
  tree->Branch("runNr",&runNr, "runNr/I");
  tree->Branch("nCh",&nCh, "nCh/I");
  tree->Branch("ch",ChannelNr, "ch[nCh]/I");
  tree->Branch("amp",amp.data(), "amp[nCh]/F");
  tree->Branch("max",max.data(), "max[nCh]/F");
  tree->Branch("min",min.data(), "min[nCh]/F");
  tree->Branch("t",t, "t[nCh]/F");
  tree->Branch("BL", BL, "BL[nCh]/F");
  tree->Branch("BL_RMS", BL_RMS, "BL_RMS[nCh]/F");
  tree->Branch("Integral_0_300", Integral_0_300, "Integral_0_300[nCh]/F");
  tree->Branch("Integral", Integral, "Integral[nCh]/F");
  tree->Branch("EventIDsamIndex",EventIDsamIndex, "EventIDsamIndex[nCh]/I");
  tree->Branch("FirstCellToPlotsamIndex",FirstCellToPlotsamIndex, "FirstCellToPlotsamIndex[nCh]/I");

    int nitem = 1;
    ifstream inList;
    TString fileName;
    inList.open(_inFileList);
    assert(inList.is_open());

    int wavePrintStatus=-1;
    int ch0PrintStatus=-1;
    int trigPrintStatus=-1;
    int signalPrintStatus=-1;
    while(inList >> fileName){
      fileName = _inDataFolder + fileName;
      cout << endl;
      cout << fileName << endl;
      FILE* pFILE = fopen(fileName.Data(),"rb");
      if (pFILE==NULL) {fputs ("File error",stderr); assert(0);}
      //cout<<" ---> File to convert : " << fileName << endl;
      fseek (pFILE , 0 , SEEK_END);
      int totFileSizeByte = ftell (pFILE);
      rewind (pFILE);
      cout<<"totFileSizeByte = "<<totFileSizeByte<<endl;
      char header[328];
      nitem=fread(header,1,328,pFILE);
      cout << "Header:\n" << header << endl;
      

      
      
      char* word;
      word = strtok(header," \n");
      while(word != NULL){
	  if(strcmp("ACQUIRED:",word) == 0){
	    word = strtok(NULL, " \n");
	    nActiveCh = atoi(word);
	    break;
	  }
	 //printf ("%s\n",word);
	 word = strtok(NULL, " \n");
      }

      if(nActiveCh>9){
	cout << endl;
	char dummy;
	nitem=fread(&dummy,1,1,pFILE);
      }


      
      int whileCounter = 0;
      while(nitem>0){ //event loop
      std::vector<TObject*> eventTrash;
      
      whileCounter++;
	nitem = fread (&EventNumber	,sizeof(int), 1,pFILE);
	nitem = fread (&EpochTime	,sizeof(double)      , 1,pFILE);
	nitem = fread (&Year		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Month		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Day		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Hour		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Minute		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Second		,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&Millisecond	,sizeof(unsigned int), 1,pFILE);
	nitem = fread (&nCh	,sizeof(unsigned int),1,pFILE); // since V2.8.14 the number of stored channels is written for each event


	if(EventNumber%10==0)printf("POS, ev, y-m-d-h-min-s-ms, nActive-nCh: %ld, %d, %d-%d-%d-%d-%d-%d-%d, %d-%d \n", ftell(pFILE), EventNumber,Year,Month,Day,Hour,Minute,Second,Millisecond,nActiveCh,nCh);

	float	MeasuredBaseline[n_Ch];
	float	AmplitudeValue[n_Ch];
	float	ComputedCharge[n_Ch];
	float	RiseTimeInstant[n_Ch];
	float	FallTimeInstant[n_Ch];
	float	RawTriggerRate[n_Ch];
	float floatR=-1;
        for(int i = 0;i<nCh;i++){
	  //printf("i, currentPositionByte %d %ld\n",i,ftell(pFILE));
	  nitem = fread (&ChannelNr[i]	       ,sizeof(int),1,pFILE);
          nitem = fread (&EventIDsamIndex[i]        ,sizeof(int),1,pFILE);
	  nitem = fread (&FirstCellToPlotsamIndex[i],sizeof(int),1,pFILE);
	  nitem = fread (&floatR,1,4,pFILE); MeasuredBaseline[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); AmplitudeValue[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); ComputedCharge[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); RiseTimeInstant[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); FallTimeInstant[i] = floatR;
	  nitem = fread (&floatR,1,4,pFILE); RawTriggerRate[i] = floatR;
	  ChannelNr[i]=i;

	  TString title("");
	  title.Form("ch %d, ev %d",i,EventNumber);
	  hCh.Reset();
	  hCh.SetTitle(title);
	  
	    for(int j = 0;j<1024;j++){
	      nitem = fread (&amplValues[i][j],sizeof(short),1,pFILE);
	      hCh.SetBinContent(j+1,-(amplValues[i][j]*coef*1000));
	    }//for 1024
	  
	  //for(int t=0;t<nActiveCh-nCh;t++){
	  //  int dummy;
	  //  nitem = fread(&dummy,sizeof(int),1,pFILE);
	  //  printf("trigger channel number: %d\n",dummy);
	  //}

	  max[i]=hCh.GetMaximum();
	  min[i]=hCh.GetMinimum();
	  amp[i]=hCh.GetMaximum();
	  
	  getBL(&hCh, BL_output,0,75);
	  BL[i] = BL_output[0];
	  BL_RMS[i] = BL_output[1];
	  
	  for(int j=1;j<=hCh.GetXaxis()->GetNbins();j++){
	    hCh.SetBinError(j,BL_RMS[i]);
	  }
	  hChtemp.at(i) = hCh;
	  if(i<=5)t[i]=CDF(&hCh,fTrigFit,0.4);
	  else t[i]=CDF(&hCh,fTrigFit,0.1);
	  Integral_0_300[i] = (hCh.Integral(1, 1024, "width")-BL[i]*1024*SP)/pe;//Calculating Integral of histogram from 0 to 300ns; starting from bin 1 (0 is the overflow bin) to bin corresponding to 300ns. Option "width" multiplies by bin-width such that the integral is independant of the binning
	  
	  
	  if(EventNumber%wavesPrintRate==0){
	    if(i==0)cWaves.cd(3);
	    if(i==1)cWaves.cd(1);
	    if(i==2)cWaves.cd(9);
	    if(i==3)cWaves.cd(7);
	    if(i==4)cWaves.cd(2);
	    if(i==5)cWaves.cd(8);
	    if(i==6)cWaves.cd(5);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000);
	    ln->SetLineColor(2);
	    ln->Draw("same");
	    //fTrigFit->DrawCopy("same");
	    eventTrash.push_back(ln);
	  }
	  

	 if(EventNumber%trigPrintRate==0&&(i<4)){
	    cTrig.cd(i+1);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000);
	    ln->SetLineColor(2);
	    ln->Draw("same");
	    //fTrigFit->DrawCopy("same");
	    eventTrash.push_back(ln);
	  }

	 if(EventNumber%signalPrintRate==0&&(i>=4&&i<=7)){
	    cSignal.cd(i+1-4);
	    hCh.DrawCopy();
	    TLine* ln = new TLine(t[i],-2000,t[i],2000); 
	    ln->SetLineColor(2);
	    ln->Draw("same");
	    eventTrash.push_back(ln);
	  }



       }//for nCh

      trigT = (t[0]+t[1]+t[2]+t[3]+t[4]+t[5])/6;
      trigTp = (t[0]+t[1]+t[4]-t[2]-t[3]-t[5])/6;
      tPMT = t[6]-trigT;
      if(tPMT<-52){
	t[6]=CDFinvert(&hChtemp.at(6),0.1);
	tPMT = t[6]-trigT;
      }
      //tPMTi = iCFD(&hChtemp.at(5),trigT-55,2,BL[5])-trigT;
      w = h_weight->GetBinContent(h_weight->FindBin(Integral_0_300[6]));
      
      Integral[6] = integral(&hChtemp.at(6),t[6]-5,t[6]+65,BL[6])/pe;
      

      trigGate = abs(*(std::max_element(t,t+6))-*(std::min_element(t,t+6)));  
      
      if(max[0]<1240&&max[1]<1240&&max[2]<1240&&max[3]<1240&&max[4]<1240&&max[5]<1240){
	isTrig=1;
	if(isTrig&&Integral_0_300[0]>2&&Integral_0_300[1]>1.8&&Integral_0_300[2]>1.3&&Integral_0_300[3]>1.4&&Integral_0_300[4]>0.8&&Integral_0_300[5]>1){
	  isTrig=1;
// // 	  if(trigT<140&&trigT>90&&trigGate<15){
// // 	    isTrig=1;
// // 	  }
// // 	  else isTrig=0;
	}
	else if(min[0]>-5&&min[1]>-5&&min[2]>-5&&min[3]>-5&&min[4]>-5&&min[5]>-5) isTrig=1;
	else isTrig=0;
      }
      else isTrig=0;
      if(isTrig==1){
	int shift = (int)((140-trigT)/SP);
	for(int j=0;j<(int)hChtemp.size();j++){
	  hChSum.at(j)->Add(&hChtemp.at(j),1);
	  hChShift_temp.at(j).Reset();
	  for(int bin=1;bin<=hCh.GetXaxis()->GetNbins()-shift;bin++){
	   hChShift_temp.at(j).SetBinContent(shift+bin,hChtemp.at(j).GetBinContent(bin));
	  }
	  hChShift.at(j)->Add(&hChShift_temp.at(j),1);
	}
      }
      
      
      //if(EventNumber%wavesPrintRate==0&&(tPMT>-10&&isTrig)){
      //if(EventNumber%wavesPrintRate==0&&(tPMT<-20&&isTrig)){
      //if(EventNumber%wavesPrintRate==0&&isTrig&&tPMT>-2&&tPMT<-1){
      if(EventNumber%wavesPrintRate==0){
	    //TString plotSaveName("");
	    //plotSaveName.Form("%s/wave-%d.png",plotSaveFolder.Data(),EventNumber);
	    if(wavePrintStatus<0){cWaves.Print((TString)(plotSaveFolder+"/waves.pdf("),"pdf");wavePrintStatus=0;}
	    else cWaves.Print((TString)(plotSaveFolder+"/waves.pdf"),"pdf");
      }
      if(EventNumber%trigPrintRate==0){
	    if(trigPrintStatus<0){cTrig.Print((TString)(plotSaveFolder+"/trig.pdf("),"pdf");trigPrintStatus=0;}
	    else cTrig.Print((TString)(plotSaveFolder+"/trig.pdf"),"pdf");
      }
      if(EventNumber%signalPrintRate==0){
	    if(signalPrintStatus<0){cSignal.Print((TString)(plotSaveFolder+"/signal.pdf("),"pdf");signalPrintStatus=0;}
	    else cSignal.Print((TString)(plotSaveFolder+"/signal.pdf"),"pdf");
      }

      tree->Fill();
      //cleanEventMemory(eventTrash);
      }//while events

    }//while
    inList.close();
    cWaves.Clear();
    cWaves.Print((TString)(plotSaveFolder+"/waves.pdf)"),"pdf");
    cCh0.Print((TString)(plotSaveFolder+"/ch0.pdf)"),"pdf");
    cTrig.Print((TString)(plotSaveFolder+"/trig.pdf)"),"pdf");
    cSignal.Print((TString)(plotSaveFolder+"/signal.pdf)"),"pdf");

  rootFile = tree->GetCurrentFile();
  rootFile->Write();
  rootFile->Close();
}

float CDF(TH1F* hWave, TF1* fTrigFit,float thr){
  float peak=hWave->GetMaximum();
  int timePos=1;
  float val = 0;
  while(val<thr*peak){
    timePos+=1;
    val = hWave->GetBinContent(timePos);
  }
  
  double x1 = SP*(timePos-1);
  double x2 = SP*(timePos);
  double y1 = hWave->GetBinContent(timePos-1);
  double y2 = hWave->GetBinContent(timePos);
  double k = (x2-x1)/(y2-y1);
  return  x1+k*(thr*peak-y1);
  
  //fit procedure
  //fTrigFit->SetParameter(1,SP*(timePos)+1);
  //fTrigFit->SetRange(SP*(timePos-6),SP*(timePos+2));
  //hWave->Fit(fTrigFit,"RNQ");
  //double p0=fTrigFit->GetParameter(0);
  //double p1=fTrigFit->GetParameter(1);
  //double p2=fTrigFit->GetParameter(2);
  //return p1-sqrt(2*p2*p2*log(p0/(0.1*abs(peak))));
}

float CDFinvert(TH1F* hWave,float thr){
  float peak=hWave->GetMaximum();
  int timePos=hWave->GetMaximumBin();
  float val = peak;
  while(val>thr*peak){
    val = hWave->GetBinContent(timePos);
    timePos-=1;
  }
  
  double x1 = SP*(timePos);
  double x2 = SP*(timePos+1);
  double y1 = hWave->GetBinContent(timePos);
  double y2 = hWave->GetBinContent(timePos+1);
  double k = (x2-x1)/(y2-y1);
  return  x1+k*(thr*peak-y1);
  
}


float iCFD(TH1F* hWave,float t,float thrNpe,float BL){
  int bin1 = hWave->GetXaxis()->FindBin(t);
  float bin1_UpEdge = hWave->GetXaxis()->GetBinUpEdge(bin1);
  float sum = (bin1_UpEdge-t)*(hWave->GetBinContent(bin1)-BL);
  int bin=bin1+1;
  while(sum<thrNpe){
      sum+=(hWave->GetBinContent(bin)-BL)*SP;
      bin+=1;
      if(bin>1024)return -999;
  }
  float bin_LowEdge = hWave->GetXaxis()->GetBinUpEdge(bin);
  return bin_LowEdge + (sum-thrNpe)/(hWave->GetBinContent(bin)-BL); 
}


float integral(TH1F* hWave,float t1,float t2,float BL){
  float BW = hWave->GetXaxis()->GetBinWidth(1);
  int bin1 = hWave->FindBin(t1);
  int bin2 = hWave->FindBin(t2);
  float c1 = hWave->GetBinContent(bin1);
  float c2 = hWave->GetBinContent(bin2);
  return hWave->Integral(bin1,bin2,"width")-BL*(t2-t1)-c1*(t1-hWave->GetXaxis()->GetBinLowEdge(bin1))-c2*(hWave->GetXaxis()->GetBinUpEdge(bin2)-t2);
}



float* getBL(TH1F* hWave, float* BL, float t1, float t2){
  /*
  Function to calculate the baseline of the given TH1F-Object.
  Input: TH1F-Object to calculate baseline for; bool isNegative; float-array BL for the output
  Output: baseline and rms of baseline written to 1st and 2nd component of BL-array
  The baseline is calculated as the mean of the values in range of (t1,t2) of the TH1F-Object. The rms
  value is also calculated from the same values.
  
  The float-array BL that is used for the output must be declared before the function call using 'float BL[2];'.
  This is to insure that the output is stored on the heap and not deleted when the memory on the stack is freed up.
   
  Dependencies: function uses C++ vector-class and needs the TMath-header
  */
  
  vector<float> amp;
  for (int i = int(t1/SP); i < int(t2/SP); i++){
    
    amp.push_back(hWave->GetBinContent(i+1));
  }
  BL[0] = TMath::Mean(amp.begin(), amp.end());
  BL[1] = TMath::RMS(amp.begin(), amp.end());
  return BL;
}


TString getRunName(TString inDataFolder){
      char* word;
      char* lastWord;
      word = strtok((char*)inDataFolder.Data(),"/");
      while(word != NULL){
	 //printf ("%s\n",word);
	 lastWord = word;
	 word = strtok(NULL, "/");
      }

      return (TString)lastWord;
}

void cleanEventMemory(std::vector<TObject*>& trash){
  for(int i=0;i<(int)trash.size();i++){
    trash.at(i)->Delete();
  }
  trash.clear();
}
