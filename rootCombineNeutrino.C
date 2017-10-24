{
  gROOT.SetBatch(0);

  TChain* chainA = new TChain("tree");

  TString base ="/data/neutrino06/phallsjo";

  // TString folderA ="/multiNeutrino/";

  TString folderA ="/multiNeutrinoFix/";

  vector<vector<double> > valuesp;
  vector<double> currentp;

  //double num = 0;
  
  //for(double cnt = 0.1; cnt < 4.1; cnt+=0.1)
  for(double cnt = 0.1; cnt < 10.1; cnt+=0.1)
  //for(int i = 1; i < 21; i++)
    {
      //double temp =(double) i;
      //cout<<temp<<endl;
      
      //double cnt = double(i)/10.0;
      //cout<<std::setprecision(2)<<cnt<<endl;

      ostringstream convert;
      convert << cnt;
      TString fileName;

      fileName = convert.str()+".root";

      chainA->Add(base+"/out/rec_out/nd_mu-CCQE"+folderA+fileName);
    }

  

  //chainA->Add("/data/neutrino06/phallsjo/out/rec_out/nd_mu-CCQE/nd_mu-CCQE_1000.root");

  /*
 TCanvas* c = new TCanvas("c","c",600,400);

 TH1D* effA = new TH1D("effA",";True Momentum (GeV/c); Fractional Efficiency",100, 0, 10);
 
 TH1D* allA = new TH1D("allA",";True Momentum (GeV/c); Fractional Efficiency",100, 0, 10);
  

 chainA->Draw("abs(MCtr_Mom)*0.001>>allA",""); //Base

 chainA->Draw("abs(MCtr_Mom)*0.001>>effA","MCtr_ID && abs(traj_Mom)>10 && abs(traj_Mom)<20000"); //ChargeID

 //chainA->Draw("abs(MCtr_Mom)*0.001>>effA","traj_Fitted && abs(traj_Mom)>10 && abs(traj_Mom)<20000"); //Reconstruction
  effA->Sumw2();
 
  allA->Sumw2();

  TGraphAsymmErrors testA;// = TGraphAsymmErrors(eff);
  testA = TGraphAsymmErrors(effA,allA);
  testA.SetLineColor(1);
  testA.SetLineWidth(2);
  testA.SetLineWidth(2);

  //testA.SetTitle("Charge reconstruction efficiency;Momentum of Muon [GeV/c];Charge eff");
  testA.SetTitle("Reconstruction Efficiency;Momentum of Muon [GeV/c];Rec eff");

  testA.Draw();
  
  TLegend* leg = new TLegend(0.8,0.1,0.9,0.3);
  leg->SetHeader("Legend");
  leg->AddEntry(&testA,"#mu_{-}","l");
  leg->Draw();
  */
 // Neutrino energy from three body interaction.

 // L -> Charged lepton. Muon, mL = 150.658 MeV/c
 // PL -> 3-Momentum of charged lepton.
 // EL -> Energy of charged lepton.
 // thetaL -> angle between these two. No! Angle of muon wrt incoming neutrino.
 // M -> Mass of stationary (proton) 938.272 MeV/c// (neutron) target mass. nu to mu- , n to p .Neutrino 939.565 MeV/c

 // Erec = (M * EL - mL^2/2-M^2/2) / (M-EL + abs(PL)cos(thetaL))

  //old  Erec = (M * EL - mL^2/2) / (M-EL + abs(PL)cos(thetaL))

 // E^2 = p^2 +(mc^2)^2;

  //chainA->Draw("trajsAngle:MC_PrimaryEng","trajsAngle>0","COLZ");

  //chainA->Draw("(939.565*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)+abs(MCtr_Mom)*cos(trajsAngle)):MC_PrimaryEng","trajsAngle>0 && abs((939.565*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)+abs(MCtr_Mom)*cos(trajsAngle)))<5000 && (939.565*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)+abs(MCtr_Mom)*cos(trajsAngle))>0","COLZ");

  //chainA->Draw("(939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle)):MC_PrimaryEng>>a","abs((939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle)))<5000 && (939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle))>100","COLZ");


  // New from https://journals.aps.org/prd/pdf/10.1103/PhysRevD.85.113008
  chainA->Draw("(939.565*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)+abs(MCtr_Mom)*cos(trajsAngle)):MC_PrimaryEng>>a","abs((939.565*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)+abs(MCtr_Mom)*cos(trajsAngle)))<5000 && (939.565*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658)+abs(MCtr_Mom)*cos(trajsAngle))>100","COLZ");
  // Energy
  //chainA->Draw("(939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle)):MC_PrimaryEng>>a","abs((939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle)))<5000 && (939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle))>100","COLZ");

//Residual

//chainA->Draw("abs(939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle))-MC_PrimaryEng:MC_PrimaryEng>>a","abs((939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle)))<5000 && (939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle))>1000","COLZ");


  //chainA->Draw("(939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle))-MC_PrimaryEng:MC_PrimaryEng>>a","abs((939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle)))<5000 && (939.565*sqrt(traj_Mom*traj_Mom+150.658*150.658)-150.658*150.658/2-939.565*939.565/2)/(939.565-sqrt(traj_Mom*traj_Mom+150.658*150.658)+abs(traj_Mom)*cos(trajsAngle))>100","COLZ");
  
  //a->SetTitle("Neutrino energy reconstruction;Neutrino energy [MeV/c];Reconstructed energy [MeV/c]");
  a->SetTitle("Neutrino energy reconstruction Prec-Ptrue ;Neutrino energy [MeV/c];Residual");
  
  a->GetYaxis()->SetTitleOffset(1.4);

  a->SetStats(false);
  TLine trend(1000,1000,4000,4000);
  trend->SetLineWidth(8);
  trend->SetLineColor(1);
  trend->Draw("same");
  

  //chainA->Draw("939.565*abs((150.658*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) - 150.658*150.658/2 -939.565*939.565/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) + MCtr_Mom*cos(trajsAngle))):MC_PrimaryEng>>a","abs(939.565*abs((150.658*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) - 150.658*150.658/2 -939.565*939.565/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) + MCtr_Mom*cos(trajsAngle))))<20000","COLZ");


  //chainA->Draw("abs((939.565*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) - 150.658*150.658/2)/(939.565-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) + abs(MCtr_Mom)*cos(trajsAngle))):MC_PrimaryEng>>a","","COLZ");


  //chainA->Draw("abs((150.658*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) - 150.658*150.658/2)/(150.658-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) + abs(MCtr_Mom)*cos(trajsAngle))):MC_PrimaryEng>>a","abs((150.658*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) - 150.658*150.658/2)/(150.658-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) + abs(MCtr_Mom)*cos(trajsAngle)))<5000 && abs((150.658*sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) - 150.658*150.658/2)/(150.658-sqrt(MCtr_Mom*MCtr_Mom+150.658*150.658) + abs(MCtr_Mom)*cos(trajsAngle)))> 700 && abs(MCtr_Mom)>10","COLZ");

  //chainA->Draw("abs((150.658*sqrt(traj_Mom*traj_Mom+150.658*150.658) - 150.658*150.658/2)/(150.658-sqrt(traj_Mom*traj_Mom+150.658*150.658) + abs(traj_Mom)*cos(trajsAngle))):MC_PrimaryEng>>a","abs((150.658*sqrt(traj_Mom*traj_Mom+150.658*150.658) - 150.658*150.658/2)/(150.658-sqrt(traj_Mom*traj_Mom+150.658*150.658) + abs(traj_Mom)*cos(trajsAngle)))<5000 && abs((150.658*sqrt(traj_Mom*traj_Mom+150.658*150.658) - 150.658*150.658/2)/(150.658-sqrt(traj_Mom*traj_Mom+150.658*150.658) + abs(traj_Mom)*cos(trajsAngle)))> 700 && abs(traj_Mom)>10","COLZ");

 //gSystem->Exit(0);
}
