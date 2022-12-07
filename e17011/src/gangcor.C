void gangcor() {

  TH1F *dis = new TH1F("dis","dis",1000,0,1000);

  TH1F *angdis = new TH1F("angdis","angdis",100,0,1000);


  Double_t angcor[180];
  Double_t angsum = 0;
  for(Int_t i=0;i<180;i++){
    Double_t ctheta = TMath::Cos(TMath::Pi()*(i/180.));
    angcor[i] = 1. + (1./8.)*ctheta*ctheta + (1./24.)*ctheta*ctheta*ctheta*ctheta;
    dis->Fill(i,angcor[i]);
    angsum += angcor[i];
    //G4cout << "ang " << ctheta << " " << angcor[deg] << " " << angsum << G4endl;
  }

  for(Int_t z=0;z< 100000;z++){

    if(z%100==0) cout << z << " done" << endl;

    Double_t prob = 1. * gRandom->Rndm();
    Double_t ang = 0;
    Double_t tempsum = 0;
    for(Int_t i=0;i<180;i++){
      tempsum += angcor[i]/angsum;
      //cout << "temps sum " << tempsum << " " << prob << " " << angcor[i]<<" " << angsum << endl;
      if(prob < tempsum) {
	ang = i;
	break;
      }
    } 
    //cout << "ang "<<ang << endl;
    angdis->Fill(ang); 
  }

  TF1 *fitf = new TF1("fitf","50*(1.+(1./8.)*(x)*(x) + (1./24.)*(x)*(x)*(x)*(x))",-1,1);
  //fitf->SetParameter(0,50);
  //fitf->Draw();
  //angdis->Fit("fitfunc");
  //angdis->Draw();
  //fitf->Draw("same");
  dis->Draw();
  angdis->Draw();
}
