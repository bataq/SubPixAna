// ROOT Macro file
// make & ./MakeCanvas
// Auther KAWABATA Tomoki
// 2017-01-30

void MakeCanvas()
{   
  double rot_angle=30;
  int energy_m=50;
  int energy_M=100;
  
  TFile *rootfile=new TFile("subpix.root");  
  
  // Cut ////////////////////////////////////////
  // energy
  TCut ecut_m=Form("ph_merge>=%d",energy_m);
  TCut ecut_M=Form("ph_merge<=%d",energy_M);
  TCut type_select = "type == 10||type == 11||type == 20||type == 21||type == 22";
  // projection region
  Double_t rad = rot_angle * TMath::Pi()/ 180;
  Int_t p_xrm = 71*cos(rad)-71*sin(rad)+143*sin(rad)-60;
  Int_t p_xrM = 71*cos(rad)-71*sin(rad)+143*sin(rad)+60;
  Int_t p_yrm = 71*sin(rad)+71*cos(rad)-20;
  Int_t p_yrM = 71*sin(rad)+71*cos(rad)+20;
  TCut xxm = Form("xx>=%d",p_xrm);
  TCut xxM = Form("xx<=%d",p_xrM);
  TCut yym = Form("yy>=%d",p_yrm);
  TCut yyM = Form("yy<=%d",p_yrM);
  // /////////////////////////////////////////////

  // Draw Region & bin ///////////////////////////
  Int_t rM = 144*sin(rad)+144*cos(rad);
  Int_t nbin =10*rM;
  // /////////////////////////////////////////////
  TCanvas *canvas=new TCanvas("SubPixCanvas","subpix",700,700);
  
  gStyle->SetOptStat(0);

  
  canvas->Divide(2,2);
  canvas->cd(1);
  TH2D *rot_map = new TH2D("subpixmap","Subpix map; CA; RA",nbin,0,rM,nbin,0,rM); 
  subpix_tree->Draw("yy:xx>>subpixmap",ecut_m && ecut_M && type_select,"colz");
  
  canvas->cd(2);
  TH2D *reg_map = new TH2D("regionmap","Projection Region; CA; RA",rM,0,rM,rM,0,rM); 
  subpix_tree->Draw("yy:xx>>regionmap",xxm && xxM && yym && yyM,"scat");
  
  canvas->cd(3);
  TH2D *pro_map = new TH2D("projection","Rotation map; CA; RA",nbin/2,0,rM,nbin/2,0,rM);
  *pro_map = *rot_map;
  pro_map->GetXaxis()->SetRangeUser(p_xrm,p_xrM);
  TH1D *h_projection = new TH1D("projection","Projection; CA",90,10,100); 
  Projection(h_projection,pro_map,p_yrm,p_yrM);//(2D hist,y min,y max)
  
  canvas->cd(4);
  WriteConditions(ecut_m,ecut_M,xxm,xxM,yym,yyM);

  std::cout<<Form(" rotation angle = %f",rot_angle)<<std::endl;
  std::cout<<"#### Cut Conditions ####"<<std::endl;
  std::cout<<Form(" ph_merge>=%d",energy_m)<<"(PH)"<<std::endl;
  std::cout<<Form(" ph_merge<=%d",energy_M)<<"(PH)"<<std::endl;
  std::cout<<Form(" xx>=%d",p_xrm)<<std::endl;
  std::cout<<Form(" xx<=%d",p_xrM)<<std::endl;
  std::cout<<Form(" yy>=%d",p_yrm)<<std::endl;
  std::cout<<Form(" yy<=%d",p_yrM)<<std::endl;
  std::cout<<" type==single && double"<<std::endl;
  std::cout<<"########################"<<std::endl;

    TFile *file = new TFile("subpix.root","update");  
    std::cout<<std::endl;
    std::cout<<"######### update file subpix.root #########"<<std::endl;  
    rot_map->Write();
    //reg_map->Write();
    h_projection->Write();
    canvas->Write();
    file->Close();
}

void Projection(TH1D *h1,TH2D *h2,Int_t first_bin,Int_t last_bin)
{
  Int_t first = h2->GetXaxis()->FindBin(first_bin);
  Int_t last = h2->GetXaxis()->FindBin(last_bin);
  h1 = h2->ProjectionX("projection",first,last);
  h1->SetTitle("Projection");
  h1->GetXaxis()->SetRangeUser(40,160);
  h1->Draw();
}

void WriteConditions(TCut ecut_m,TCut ecut_M,TCut xxm,TCut xxM,TCut yym,TCut yyM){

  //  TText *filename=new TText(0.1,0.6,rootfile);

  TText *t0=new TText(0.09,0.551,"Cut Conditions");
  TText *t1=new TText(0.1,0.5,ecut_m);
  TText *t2=new TText(0.1,0.45,ecut_M);
  TText *t3=new TText(0.1,0.4,xxm);
  TText *t4=new TText(0.1,0.35,xxM);
  TText *t5=new TText(0.1,0.3,yym);
  TText *t6=new TText(0.1,0.25,yyM);
  TText *t7=new TText(0.1,0.2,"type==single&&double");
  
  t0->Draw();
  t1->Draw();
  t2->Draw();
  t3->Draw();
  t4->Draw();
  t5->Draw();
  t6->Draw();
  t7->Draw();
}
