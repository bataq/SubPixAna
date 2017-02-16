// C++ file
// Auther KAWABATA Tomoki
// Date 2017-02-16

#include <TROOT.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TTree.h>
#include <TFile.h>
#include <TH2.h>
#include <TString.h>
#include <TRandom.h>

#include <vector>
#include <math.h>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <boost/program_options.hpp>
using namespace boost::program_options;

//#define __CINT__


void GetRandomRange(Double_t *r,Double_t *Rs,Double_t *Rw,int SINGLE,int DOUBLE);

#ifdef __CINT__
int SubPixAna(const TString rootFile){
#else
int SubPixAna(const TString rootFile);

int main(int argc, char* argv[]){

  // Option Analysis //
  TString dataFile;
  options_description opt("Options");
  opt.add_options()
    ("help,h","HELP")
    ("file,f",value<std::string>(),"Event List ROOT File");
  variables_map vm;
  store(parse_command_line(argc, argv, opt), vm);
  notify(vm);
  if( vm.count("help") || !vm.count("file") ){
    std::cout << "Usage  :./SubPixAna -f <Event List Root File>" << std::endl;
    std::cout << opt << std::endl;
    return 0;
  }else{
    std::string file = vm["file"].as<std::string>();
    dataFile=file;
    std::cout << "dataFile is "<<dataFile<<std::endl;
      return 0;
  }
 SubPixAna(dataFile);

  return 0;
}

 int SubPixAna(const TString rootFile){

#endif

  const int SINGLE=1;
  const int DOUBLE=1;
  const int PH1=0;

  if(!rootFile){
    std::cerr<<"Usage: SubPixAna(root file name)"<<std::endl;
    return 0;
  }

  std::cout<<"rootFile name is "<<rootFile<<std::endl;
  TFile *f = new TFile(rootFile,"read");
  if(!f){
    std::cerr<<"Error: no such file." << std::endl;
    return 0;
  }

  gROOT->ProcessLine("#include <vector>");
  
  TTree *etr=(TTree*)f->Get("tree_evlist");
  UShort_t        ra;
  UShort_t        ca;
  UShort_t        type;
  Float_t         ph_merge;
  std::vector<float>   *vortex_ph=0;
  TBranch        *b_ra;
  TBranch        *b_ca;
  TBranch        *b_type;
  TBranch        *b_ph_merge;
  TBranch        *b_vortex_ph;
  etr->SetBranchAddress("ra", &ra, &b_ra);
  etr->SetBranchAddress("ca", &ca, &b_ca);
  etr->SetBranchAddress("type", &type, &b_type);
  etr->SetBranchAddress("ph_merge", &ph_merge, &b_ph_merge);
  etr->SetBranchAddress("vortex_ph", &vortex_ph, &b_vortex_ph);
  etr->SetBranchStatus("*",0);
  etr->SetBranchStatus("ra",1);
  etr->SetBranchStatus("ca",1);
  etr->SetBranchStatus("type",1);
  etr->SetBranchStatus("ph_merge",1);
  etr->SetBranchStatus("vortex_ph",1);

  Int_t allEntries=etr->GetEntries();
    
  TTree *htr=(TTree*)f->Get("tree_header");
  std::vector<std::string>  *value=0;
  TBranch        *b_value;
  htr->SetBranchAddress("value", &value, &b_value);
  htr->SetBranchStatus("*",0);
  htr->SetBranchStatus("value",1);

  TFile *sf = new TFile("subpix.root","recreate");
  std::cout<<"## make file subpix.root ##"<<std::endl;  
  TTree *tr = new TTree("tree_subpix","Sub Pixel Analyzed Tree");
  
  // variables for SubPixTree 
  Double_t sca,sra;
  Int_t stype;
  Double_t r,Rs,Rw;
  Double_t randSX,randSY; // Rondom single x & y direction
  Double_t randD;
  htr->GetEntry(0);
  //  Double_t spTh=std::stod(value->at(5)); gcc v4.4.7 では使用不可
  std::stringstream ss(value->at(5));
  Double_t spTh;
  ss>>spTh;
  Double_t cor=TMath::Abs((PH1+spTh)/(PH1-spTh));// correction factor
  //Double_t cor=1;
  tr->Branch("sca",&sca);
  tr->Branch("sra",&sra);
  tr->Branch("stype",&stype);

  gRandom->SetSeed(time(NULL));
  GetRandomRange(&r,&Rs,&Rw,SINGLE,DOUBLE);//初期化
  randSX = gRandom->Uniform(-Rs,Rs);
  randSY = gRandom->Uniform(-Rs,Rs);
  randD = gRandom->Uniform(-Rw,Rw);
  for(int i=0; i<allEntries;i++){
    etr->GetEntry(i);
    if (type==20 || type==21){
      if(vortex_ph->at(2) >= spTh){
	sca = (Double_t)ca - 0.5 - (r*(vortex_ph->at(2) - vortex_ph->at(0)) / (vortex_ph->at(2)+vortex_ph->at(0)))*cor;
	sra = (Double_t)ra + randD;
	stype = 22;
      }else if(vortex_ph->at(6) >= spTh){ 
	sca = (Double_t)ca + 0.5 - (r*(vortex_ph->at(0) - vortex_ph->at(6)) / (vortex_ph->at(0)+vortex_ph->at(6)))*cor;
	sra = (Double_t)ra + randD;
	stype = 22;
      }else if(vortex_ph->at(4) >= spTh){ 
	sca = (Double_t)ca + randD;
	sra = (Double_t)ra - 0.5 - (r*(vortex_ph->at(4) - vortex_ph->at(0)) / (vortex_ph->at(4)+vortex_ph->at(0)))*cor;
	stype = type;
      }else if(vortex_ph->at(8) >= spTh){ 
	sca = (Double_t)ca + randD;
	sra = (Double_t)ra + 0.5 - (r*(vortex_ph->at(0) - vortex_ph->at(8)) / (vortex_ph->at(0)+vortex_ph->at(8)))*cor;
	stype = type;
      }else{
	std::cout<<"???\n";
	continue;
      }
    }else if (type == 10 || type == 11){
      sca= ca + randSX;
      sra= ra + randSY;
      stype = type;
    }else{
      continue;
    }
      tr->Fill();     
    }
  tr->Write();
  sf->Close();
  return 0;
}

 void GetRandomRange(Double_t *r,Double_t *Rs,Double_t *Rw,int SINGLE,int DOUBLE)
{
  Double_t s = (Double_t)SINGLE;
  Double_t w = (Double_t)DOUBLE;
  Double_t spw=(Double_t)SINGLE/(Double_t)DOUBLE;
  *r = (s + w - sqrt(pow(s+w,2)-w*(2*s+w)))/(2*(2*s+w));
  *Rs = sqrt(spw*(*r)*(1-2*(*r))); // range of single pixel count region [-Rs,Rs]
  *Rw = (1-2*(*r))/2; // range of double pixel count region [-Rw,Rw]
}

 /* Using TTreeReader
// TTreeReader is able to use from ROOT v6.

int SubPixAna(const TString rootFile){

  const int SINGLE=1;
  const int DOUBLE=1;
  const int PH1=0;

  if(!rootFile){
    std::cerr<<"Usage: SubPixAna(root file name)"<<std::endl;
    return 0;
  }

  TFile *f=TFile::Open(rootFile);
  if(!f){
    std::cerr<<"Error: no such file." << std::endl;
    return 0;
  }

  TTreeReader evlReader("tree_evlist",f);
  TTreeReaderValue<UShort_t> ra(evlReader,"ra");
  TTreeReaderValue<UShort_t> ca(evlReader,"ca");
  TTreeReaderValue<UShort_t> type(evlReader,"type");
  TTreeReaderValue<Float_t> ph_merge(evlReader,"ph_merge");
  TTreeReaderValue<std::vector<float>> vortex_ph(evlReader,"vortex_ph");

  TTreeReader hdReader("tree_header",f);
  TTreeReaderValue<std::vector<string>> value(hdReader,"value");

  TFile *file = new TFile("subpix.root","recreate");
  std::cout<<"## make file subpix.root ##"<<std::endl;  
  TTree *tr = new TTree("tree_subpix","Sub Pixel Analyzed Tree");

  // variables for SubPixTree 
  Double_t sca,sra;
  Int_t stype;
  Double_t r,Rs,Rw;
  Double_t randSX,randSY; // Rondom single x & y direction
  Double_t randD;
  //  Double_t spTh=GetThreshold(value->at(5));

  hdReader.Next();

  Double_t spTh=std::stod(value->at(5));

  Double_t cor=TMath::Abs((PH1+spTh)/(PH1-spTh));// correction factor
  //Double_t cor=1;
  tr->Branch("sca",&sca);
  tr->Branch("sra",&sra);
  tr->Branch("stype",&stype);

  gRandom->SetSeed(time(NULL));
  GetRandomRange(&r,&Rs,&Rw,SINGLE,DOUBLE);//初期化
  randSX = gRandom->Uniform(-Rs,Rs);
  randSY = gRandom->Uniform(-Rs,Rs);
  randD = gRandom->Uniform(-Rw,Rw);

  while(evlReader.Next()){
    if (*type==20 || *type==21){
      if(vortex_ph->at(2) >= spTh){
	sca = (Double_t)*ca - 0.5 - (r*(vortex_ph->at(2) - vortex_ph->at(0)) / (vortex_ph->at(2)+vortex_ph->at(0)))*cor;
	sra = (Double_t)*ra + randD;
	stype = 22;
      }else if(vortex_ph->at(6) >= spTh){ 
	sca = (Double_t)*ca + 0.5 - (r*(vortex_ph->at(0) - vortex_ph->at(6)) / (vortex_ph->at(0)+vortex_ph->at(6)))*cor;
	sra = (Double_t)*ra + randD;
	stype = 22;
      }else if(vortex_ph->at(4) >= spTh){ 
	sca = (Double_t)*ca + randD;
	sra = (Double_t)*ra - 0.5 - (r*(vortex_ph->at(4) - vortex_ph->at(0)) / (vortex_ph->at(4)+vortex_ph->at(0)))*cor;
	stype = *type;
      }else if(vortex_ph->at(8) >= spTh){ 
	sca = (Double_t)*ca + randD;
	sra = (Double_t)*ra + 0.5 - (r*(vortex_ph->at(0) - vortex_ph->at(8)) / (vortex_ph->at(0)+vortex_ph->at(8)))*cor;
	stype = *type;
      }else{
	std::cout<<"???\n";
	continue;
      }
    }else if (*type == 10 || *type == 11){
      sca= *ca + randSX;
      sra= *ra + randSY;
      stype = *type;
    }else{
      continue;
    }
      tr->Fill();     
    }
  tr->Write();
  f->Close();
  return 0;
  }*/
