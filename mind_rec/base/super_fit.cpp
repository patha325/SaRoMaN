#include <super_fit.h>
#include <TMath.h>

using namespace bhep;///

//**********************************************************************
super_fit::super_fit() {
  //**********************************************************************
}

//***********************************************************************
super_fit::~super_fit() {
  //***********************************************************************
}


/***************************************************************************************/
double super_fit::RangeMomentum(double length,double nodeZ){
  /***************************************************************************************/
  // Get momentum depending on length. Did it go through each submodule? Get specific wFe
  // for each submodule.

  std::map<dict::Key,vector<double> > moduleDataMap = _supergeom.getModuleDataMap();
  double p = 0;

  for (std::map<dict::Key,vector<double> >::iterator it=moduleDataMap.begin();
       it!=moduleDataMap.end(); ++it)
    {
      double module_pos = it->second[0];
      double module_half_size = it->second[1];
      double wFe = it->second[2];

      //Sanity checking not outside range forward 
      if(!((module_pos-module_half_size)<(nodeZ+length)))
	{ 
	  continue;
	}

      // Through the whole module
      else if((nodeZ + length) > (module_pos + module_half_size))
	{
	   p += 2*module_half_size * wFe;
	}
      // Stop in the module
      else if((nodeZ + length) < (module_pos + module_half_size))
	{
	  p+=((length + nodeZ) - (module_pos - module_half_size))*wFe;
	}	 
    }

  return p;

}

/***************************************************************************************/
double super_fit::MomentumFromCurvature2(const Trajectory& traj, int startPoint,double minP){
  /***************************************************************************************/

  // Unsure which hit is first.

  bool startLow = false;

  if(traj.nodes()[startPoint]->measurement().vector()[2] < traj.nodes()[traj.size()-1]->measurement().vector()[2])
    {
      startLow = true;
    }

  // Find straight line between 2 and 3 point.
  // y = k*z+m

  int start = startPoint;
  int point = traj.size()-3;

  double scalar = dot(_supergeom.getRawBField(traj.nodes()[start]->measurement().vector()),
		     _supergeom.getRawBField(traj.nodes()[point+2]->measurement().vector()));

  while(scalar < 0 && start <= point && 0 <=point)
    {
      if(startLow)
	{
	  point--;
	}
      else
	{
	  start++;
	}
      scalar = dot(_supergeom.getRawBField(traj.nodes()[start]->measurement().vector()),
		     _supergeom.getRawBField(traj.nodes()[point+2]->measurement().vector()));
    }


  double k1 = 0;

  while(k1 == 0 && start <= point)
    {
      k1 = ((traj.nodes()[start]->measurement().vector()[1]-traj.nodes()[start+1]->measurement().vector()[1])/
	    (traj.nodes()[start]->measurement().vector()[2]-traj.nodes()[start+1]->measurement().vector()[2]) +
	    (traj.nodes()[start+1]->measurement().vector()[1]-traj.nodes()[start+2]->measurement().vector()[1])/
	    (traj.nodes()[start+1]->measurement().vector()[2]-traj.nodes()[start+2]->measurement().vector()[2]))/2;
      
      start++;
    }

  double k2 = 0;

  while(k2 == 0 & 0 <= point)
    {
      k2 = ((traj.nodes()[point]->measurement().vector()[1]-traj.nodes()[point+1]->measurement().vector()[1])/
	    (traj.nodes()[point]->measurement().vector()[2]-traj.nodes()[point+1]->measurement().vector()[2]) +
	    (traj.nodes()[point+1]->measurement().vector()[1]-traj.nodes()[point+2]->measurement().vector()[1])/
	    (traj.nodes()[point+1]->measurement().vector()[2]-traj.nodes()[point+2]->measurement().vector()[2]))/2;

      point--;
    }

  cout<<"start: "<<start<<endl;
  cout<<"point: "<<point<<endl;
  cout<<"size: "<<traj.size()<<endl;

  double m1 = traj.nodes()[start+1]->measurement().vector()[1] - k1* (traj.nodes()[start+1]->measurement().vector()[2]);  
  double m2 = traj.nodes()[point+1]->measurement().vector()[1] - k2* (traj.nodes()[point+1]->measurement().vector()[2]);
    
  // Find their orthogonal lines.

  double k3 = -1/k1;
  double k4 = -1/k2;

  double m3 = traj.nodes()[start+1]->measurement().vector()[1] - k3* (traj.nodes()[start+1]->measurement().vector()[2]); 
  double m4 = traj.nodes()[point+1]->measurement().vector()[1] - k4* (traj.nodes()[point+1]->measurement().vector()[2]); 

  // Find intersection

  double zc = (m3-m4)/(k4-k3);

  double yc = k3 *zc + m3;

  double r = sqrt((zc-traj.nodes()[1]->measurement().vector()[2])*(zc-traj.nodes()[1]->measurement().vector()[2])
		  + (yc-traj.nodes()[1]->measurement().vector()[1])*(yc-traj.nodes()[1]->measurement().vector()[1]));

  double r2 = sqrt((zc-traj.nodes()[point+1]->measurement().vector()[2])*(zc-traj.nodes()[point+1]->measurement().vector()[2])
		  + (yc-traj.nodes()[point+1]->measurement().vector()[1])*(yc-traj.nodes()[point+1]->measurement().vector()[1]));

  //EVector B = _supergeom.getBField(pos);

  double field = _supergeom.getRawBField(traj.nodes()[0]->measurement().vector())[0]*_supergeom.getBScaleAvr();

  double correction = 3.5;

  double finalP = 300*fabs(field) * correction *(r+r2)/2;

  cout<<"Momentum highPt final super"<<finalP<<endl;

  if(finalP<minP)
    {
      // Try to recursivly find a better value.
      if((start +1)< point)
	{
	  finalP =  MomentumFromCurvature(traj,start+1,minP);
	}
      else
	{
	  finalP = minP;
	}
    }
  // Need a better estimate of a maximum value for the detector.
  if(finalP > 14000)
    {
      // Try to recursivly find a better value.
      if((start +1)< point)
	{
	  finalP =  MomentumFromCurvature(traj,start+1,minP);
	}
      else
	{
	  finalP = 14000;
	}
    }

  return finalP;
}

//*****************************************************************************
double fitf3(Double_t *x,Double_t *par) { 
  //*****************************************************************************

  double z = x[0];

  double fitval = par[0] + par[1]*z+par[2]*z*z+par[3]*z*z*z+par[4]*z*z*z*z;

  return fitval;
}


//void fitter::ComputeMomFromParabola(const Trajectory& traj, int nplanes, int firsthit, EVector& V){
//MomentumFromCurvature(const Trajectory& startTrack, int startPoint,double minP){
//*****************************************************************************
double super_fit::MomentumFromCurvature(const Trajectory& traj, int startPoint,double minP){
  //*****************************************************************************

  //cout<<"in new MomentumFromCurvature"<<endl;

  //Some catchers for pointless returns.
  int fitcatch;
  //
  int nfit, sign;
  int fitRange[3];
  int nplanes = traj.size()-2;
  int firsthit = startPoint;
  const int fitpoints =  traj.size()-2;//nplanes - firsthit;
  
  double xpos[fitpoints], ypos[fitpoints], zpos[fitpoints];
  double upos[fitpoints], vpos[fitpoints];

  int pos = 0;

  EVector currentpos = EVector(3,0);
  EVector currentB   = EVector(3,0);
  EVector V = EVector(6,0);
  EVector z = EVector(3,0);
  z[2] = 1;
  double Bmean=0;
  for( int ipoint=firsthit; ipoint < nplanes; ipoint ++ ){
    xpos[pos] = traj.node(ipoint).measurement().position()[0];
    ypos[pos] = traj.node(ipoint).measurement().position()[1];
    zpos[pos] = traj.node(ipoint).measurement().position()[2]
      - traj.node(firsthit).measurement().position()[2];
    currentpos[0] = traj.node(ipoint).measurement().position()[0];
    currentpos[1] = traj.node(ipoint).measurement().position()[1];
    currentpos[2] = 0.;
    currentB = _supergeom.getBField(currentpos);
    //double field = _supergeom.getRawBField(track.nodes()[0]->measurement().vector())[0]*_supergeom.getBScaleAvr();
    upos[pos] = xpos[pos] > 0 ? asin(ypos[pos]/currentpos.norm())
      : -asin(ypos[pos]/currentpos.norm());
    vpos[pos] = dot(currentpos,crossprod(z,currentB))/currentB.norm();
    Bmean += currentB.norm();
    ++pos;
  }
  Bmean /= pos;
  Bmean /= tesla;
  
  if (fitpoints <= 15) { nfit = 1; fitRange[0] = fitpoints;}
  else if (fitpoints <= 40) { 
    nfit = 2;
    fitRange[0] = 15; fitRange[1] = (int)(0.7*fitpoints);
  }
  else if (fitpoints > 40) { 
    nfit = 3;
    fitRange[0] = 15; fitRange[1] = (int)(fitpoints/2); fitRange[2] = (int)(0.7*fitpoints);
  }
  for (int ifit = 0;ifit < nfit;ifit++) {
    //cout<<"for ifit "<<ifit<<endl;
    TGraph *trajFitXZ = new TGraph(fitRange[ifit],zpos, xpos);
    TGraph *trajFitYZ = new TGraph(fitRange[ifit],zpos, ypos);
    TGraph *trajFitUZ = new TGraph(fitRange[ifit],zpos, upos);
    TGraph *trajFitVZ = new TGraph(fitRange[ifit],zpos, vpos);
    
    TF1 *func = new TF1("fit",fitf3,-3,3,3);
    func->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func->SetParNames("a", "b", "c", "d", "e");
    
    TF1 *func2 = new TF1("fit2",fitf3,-3,3,3);
    func2->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func2->SetParNames("f", "g", "h", "i", "j");

    TF1 *func3 = new TF1("fit3",fitf3,-3,3,3);
    func3->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func3->SetParNames("a1", "b1", "c1", "d1", "e1");
    
    TF1 *func4 = new TF1("fit4",fitf3,-3,3,3);
    func4->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func4->SetParNames("f1", "g1", "h1", "i1", "j1");

    fitcatch = trajFitXZ->Fit("fit", "QN");
    fitcatch = trajFitYZ->Fit("fit2", "QN");
    fitcatch = trajFitUZ->Fit("fit3", "QN");
    fitcatch = trajFitVZ->Fit("fit4", "QN");
    
    double b = func->GetParameter(1);
    double c = func->GetParameter(2);
    double g = func2->GetParameter(1);
    /*double f = func2->GetParameter(0);
      double a = func->GetParameter(0);
      double h = func2->GetParameter(2);  
      double a1 = func3->GetParameter(0);
      double b1 = func3->GetParameter(1);
      double c1 = func3->GetParameter(2);  
      double f1 = func4->GetParameter(0);*////
    double g1 = func4->GetParameter(1);
    double h1 = func4->GetParameter(2);  
    
    //cout<<"after fitters"<<endl;

    if (ifit == 0) {

      V[4] = g;   //func2->GetParameter(1);
      V[3] = b;

      if (h1!=0) {
	V[5] = 1./(-0.3*Bmean*pow((1+g1*g1),3./2.)/
		   (2*h1)*0.01);
	V[5] /= GeV;
	sign = (int)( V[5]/fabs( V[5] ));
      } else V[5] = 0;
    } else {
      if ((int)(-c/fabs(c)) == sign) {
	V[4] = g;
	V[3] = b;
	V[5] = 1/(-0.3*Bmean*pow((1+g1*g1),3./2.)/(2*h1)*0.01);
	V[5] /= GeV;
      } else break;
    }
    //cout<<"deleter death"<<endl;
    delete trajFitXZ;
    delete trajFitYZ;
    delete trajFitUZ;
    delete trajFitVZ;
  
    delete func;
    delete func2;
    delete func3;
    delete func4;
  }
  //cout<<"before return"<<endl;
  //return fabs((1/V[5]/10.0));


  double meansign = _supergeom.getDetectorModel()->CalculateChargeMomentum();

  //return fabs(meansign);

  return meansign;

  //std::cout<<"Momentum guess from polynomial fit: p/q = "<<1./V[5]<<std::endl;
}




//***********************************************************************
double super_fit::MomentumFromCurvature3(const Trajectory& startTrack, int startPoint,double minP){
  //***********************************************************************

  int fitcatcher;
  int nMeas = startTrack.size();
  bool removeHits = false;

  Trajectory track = startTrack;

  track.sort_nodes(RP::z, -1);

  // Run this with different track starts! Will make really bad guesses when we change field.
  // Either from start or from the end.

  //Try both start from back, only 4 hits when pt to large?

  //if (nMeas > 4) nMeas = 4;

  //if(nMeas > 16)
  //removeHits = true;

  // If we pass through the detector, remove the last 2? Plane hits.

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;

  double firstNodeZ = track.nodes()[nMeas-1]->measurement().position()[2];

  int start = nMeas-2;
  int end = 0;

  if(removeHits)
    {
      double firstNodeZ = track.nodes()[nMeas-3]->measurement().position()[2];
      pos[0] = x[nMeas-3] = track.nodes()[nMeas-3]->measurement().vector()[0];
      pos[1] = y[nMeas-3] = track.nodes()[nMeas-3]->measurement().vector()[1];
      pos[2] = z[nMeas-3] = track.nodes()[nMeas-3]->measurement().vector()[2];
      start = nMeas-4;
      end = 2;
    }
  else
    {
      pos[0] = x[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[0];
      pos[1] = y[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[1];
      pos[2] = z[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[2];
    }
  


  EVector prevB(3,0);

  for (int iMeas = start;iMeas >= end;iMeas--){
    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().position()[2];

    // get the b-field from the previous step
    EVector B = _supergeom.getBField(pos);

    if(dot(prevB,B) < 0)
      {
	break;
      }
    prevB=B;

    pos[0] = x[iMeas];  pos[1] = y[iMeas];   pos[2] = z[iMeas];

    if(crossprod(Z,B).norm() ==0 )
      {
	u[iMeas] = 0;
      }
    else
      {
	u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
      }
  }

  TGraph *gr = new TGraph((const int)minindex, z, u);
  
  TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
  fun->SetParameters(0.,0.001,0.001);

  fitcatcher = gr->Fit("parfit2", "QN");

  // Positive particle -> Positive parameters -> Qtilde +
  // Negative Particle -> Negative parameters -> Qtilde -

  double field = _supergeom.getRawBField(track.nodes()[0]->measurement().vector())[0]*_supergeom.getBScaleAvr();

  double qtilde = 10.0/3.0 * 300*pow(1+pow(fun->GetParameter(1),2),3./2.)/
    (2*fun->GetParameter(2)) * field;


  if(fabs(qtilde) > 14000)
    qtilde = 14000;


  return fabs(qtilde);
}



//***********************************************************************
double super_fit::CalculateCharge(const Trajectory& startTrack) {
  //***********************************************************************
  /*
  int fitcatcher;
  int nMeas = startTrack.size();
  bool removeHits = false;

  Trajectory track = startTrack;

  track.sort_nodes(RP::z, -1);

  // Run this with different track starts! Will make really bad guesses when we change field.
  // Either from start or from the end.

  //Try both start from back, only 4 hits when pt to large?

  //if (nMeas > 4) nMeas = 4;

  //if(nMeas > 16)
  //removeHits = true;

  // If we pass through the detector, remove the last 2? Plane hits.

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas];
  int minindex = nMeas;
  double minR = 99999.999, pdR = 0.0, sumdq=0;

  double firstNodeZ = track.nodes()[nMeas-1]->measurement().position()[2];

  int start = nMeas-2;
  int end = 0;

  if(removeHits)
    {
      double firstNodeZ = track.nodes()[nMeas-3]->measurement().position()[2];
      pos[0] = x[nMeas-3] = track.nodes()[nMeas-3]->measurement().vector()[0];
      pos[1] = y[nMeas-3] = track.nodes()[nMeas-3]->measurement().vector()[1];
      pos[2] = z[nMeas-3] = track.nodes()[nMeas-3]->measurement().vector()[2];
      start = nMeas-4;
      end = 2;
    }
  else
    {
      pos[0] = x[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[0];
      pos[1] = y[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[1];
      pos[2] = z[nMeas-1] = track.nodes()[nMeas-1]->measurement().vector()[2];
    }
  


  EVector prevB(3,0);

  for (int iMeas = start;iMeas >= end;iMeas--){
    x[iMeas] = track.nodes()[iMeas]->measurement().vector()[0];
    y[iMeas] = track.nodes()[iMeas]->measurement().vector()[1];
    z[iMeas] = track.nodes()[iMeas]->measurement().position()[2];

    // get the b-field from the previous step
    EVector B = _supergeom.getBField(pos);

    if(dot(prevB,B) < 0)
      {
	break;
      }
    prevB=B;

    pos[0] = x[iMeas];  pos[1] = y[iMeas];   pos[2] = z[iMeas];

    if(crossprod(Z,B).norm() ==0 )
      {
	u[iMeas] = 0;
      }
    else
      {
	u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
      }
  }

  TGraph *gr = new TGraph((const int)minindex, z, u);
  
  TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-3,3);
  fun->SetParameters(0.,0.001,0.001);

  fitcatcher = gr->Fit("parfit2", "QN");

  // Positive particle -> Positive parameters -> Qtilde +
  // Negative Particle -> Negative parameters -> Qtilde -

  double qtilde = 0.3*pow(1+pow(fun->GetParameter(1),2),3./2.)/
    (2*fun->GetParameter(2));

  int meansign = (int)(qtilde/fabs(qtilde));
  

  // Better to be wrong than to return 0.
  if(meansign == 0)
    {
      meansign = 1;
    }
  */

  //double meansign = _supergeom.getDetectorModel()->CalculateCharge();

  double meansign = _supergeom.getDetectorModel()->CalculateChargeMomentum();

  meansign/=fabs(meansign);

  //double meansign =Helper(startTrack,0,0);
  
  return meansign;
}

//*****************************************************************************
double super_fit::Helper(const Trajectory& traj, int startPoint,double minP){
  //*****************************************************************************

  //cout<<"in new MomentumFromCurvature"<<endl;

  //Some catchers for pointless returns.
  int fitcatch;
  //
  int nfit, sign;
  int fitRange[3];
  int nplanes = traj.size()-2;
  int firsthit = startPoint;
  const int fitpoints =  traj.size()-2;//nplanes - firsthit;
  
  double xpos[fitpoints], ypos[fitpoints], zpos[fitpoints];
  double upos[fitpoints], vpos[fitpoints];

  int pos = 0;

  EVector currentpos = EVector(3,0);
  EVector currentB   = EVector(3,0);
  EVector V = EVector(6,0);
  EVector z = EVector(3,0);
  z[2] = 1;
  double Bmean=0;
  for( int ipoint=firsthit; ipoint < nplanes; ipoint ++ ){
    xpos[pos] = traj.node(ipoint).measurement().position()[0];
    ypos[pos] = traj.node(ipoint).measurement().position()[1];
    zpos[pos] = traj.node(ipoint).measurement().position()[2]
      - traj.node(firsthit).measurement().position()[2];
    currentpos[0] = traj.node(ipoint).measurement().position()[0];
    currentpos[1] = traj.node(ipoint).measurement().position()[1];
    currentpos[2] = 0.;
    currentB = _supergeom.getBField(currentpos);
    //double field = _supergeom.getRawBField(track.nodes()[0]->measurement().vector())[0]*_supergeom.getBScaleAvr();
    upos[pos] = xpos[pos] > 0 ? asin(ypos[pos]/currentpos.norm())
      : -asin(ypos[pos]/currentpos.norm());
    vpos[pos] = dot(currentpos,crossprod(z,currentB))/currentB.norm();
    Bmean += currentB.norm();
    ++pos;
  }
  Bmean /= pos;
  Bmean /= tesla;
  
  if (fitpoints <= 15) { nfit = 1; fitRange[0] = fitpoints;}
  else if (fitpoints <= 40) { 
    nfit = 2;
    fitRange[0] = 15; fitRange[1] = (int)(0.7*fitpoints);
  }
  else if (fitpoints > 40) { 
    nfit = 3;
    fitRange[0] = 15; fitRange[1] = (int)(fitpoints/2); fitRange[2] = (int)(0.7*fitpoints);
  }
  for (int ifit = 0;ifit < nfit;ifit++) {
    //cout<<"for ifit "<<ifit<<endl;
    TGraph *trajFitXZ = new TGraph(fitRange[ifit],zpos, xpos);
    TGraph *trajFitYZ = new TGraph(fitRange[ifit],zpos, ypos);
    TGraph *trajFitUZ = new TGraph(fitRange[ifit],zpos, upos);
    TGraph *trajFitVZ = new TGraph(fitRange[ifit],zpos, vpos);
    
    TF1 *func = new TF1("fit",fitf3,-3,3,3);
    func->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func->SetParNames("a", "b", "c", "d", "e");
    
    TF1 *func2 = new TF1("fit2",fitf3,-3,3,3);
    func2->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func2->SetParNames("f", "g", "h", "i", "j");

    TF1 *func3 = new TF1("fit3",fitf3,-3,3,3);
    func3->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func3->SetParNames("a1", "b1", "c1", "d1", "e1");
    
    TF1 *func4 = new TF1("fit4",fitf3,-3,3,3);
    func4->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func4->SetParNames("f1", "g1", "h1", "i1", "j1");

    fitcatch = trajFitXZ->Fit("fit", "QN");
    fitcatch = trajFitYZ->Fit("fit2", "QN");
    fitcatch = trajFitUZ->Fit("fit3", "QN");
    fitcatch = trajFitVZ->Fit("fit4", "QN");
    
    double b = func->GetParameter(1);
    double c = func->GetParameter(2);
    double g = func2->GetParameter(1);
    /*double f = func2->GetParameter(0);
      double a = func->GetParameter(0);
      double h = func2->GetParameter(2);  
      double a1 = func3->GetParameter(0);
      double b1 = func3->GetParameter(1);
      double c1 = func3->GetParameter(2);  
      double f1 = func4->GetParameter(0);*////
    double g1 = func4->GetParameter(1);
    double h1 = func4->GetParameter(2);  
    
    //cout<<"after fitters"<<endl;

    if (ifit == 0) {

      V[4] = g;   //func2->GetParameter(1);
      V[3] = b;

      if (h1!=0) {
	V[5] = 1./(-0.3*Bmean*pow((1+g1*g1),3./2.)/
		   (2*h1)*0.01);
	V[5] /= GeV;
	sign = (int)( V[5]/fabs( V[5] ));
      } else V[5] = 0;
    } else {
      if ((int)(-c/fabs(c)) == sign) {
	V[4] = g;
	V[3] = b;
	V[5] = 1/(-0.3*Bmean*pow((1+g1*g1),3./2.)/(2*h1)*0.01);
	V[5] /= GeV;
      } else break;
    }
    //cout<<"deleter death"<<endl;
    delete trajFitXZ;
    delete trajFitYZ;
    delete trajFitUZ;
    delete trajFitVZ;
  
    delete func;
    delete func2;
    delete func3;
    delete func4;
  }

  //cout<<"before return"<<endl;
  return ((1/V[5]/10.0)/fabs((1/V[5]/10.0)));
  //std::cout<<"Momentum guess from polynomial fit: p/q = "<<1./V[5]<<std::endl;
}
