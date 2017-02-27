#include <mind/detector.h>

class sortDetectorsByZ{
public:
  bool  operator()(SubDetector* det1, SubDetector* det2 ){
    if (det1->GetZ() < det2->GetZ()) return true;
    return false;
  } 
};

class sortPlanesByZ{
public:
  bool  operator()(plane_info* plane1, plane_info* plane2 ){
    if (plane1->GetZ() < plane2->GetZ()) return true;
    return false;
  } 
};

//*************************************************************
Detector::Detector() {
//*************************************************************
  
}


//*************************************************************
Detector::~Detector() {
//*************************************************************


}

//*************************************************************
SubDetector* Detector::GetSubDetector(dict::Key vol_name) {
//*************************************************************

  SubDetector* ret_val = NULL;

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      if(_subDetectorVec[i]->GetName()==vol_name)
	{
	  ret_val=_subDetectorVec[i];
	  break;
	}
    }

  return ret_val;
}

//*************************************************************
SubDetector* Detector::GetSubDetector(double zPosition) {
//*************************************************************

  SubDetector* ret_val = NULL;

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      //cout<<"Detectormin"<<_subDetectorVec[i]->GetZMin()<<endl;
      //cout<<(_subDetectorVec[i]->GetZMin()<zPosition)<<endl;
      //cout<<(_subDetectorVec[i]->GetZMax()>zPosition)<<endl;
      //cout<<"Detectormax"<<_subDetectorVec[i]->GetZMax()<<endl;
      //cout<<"zPos"<<zPosition<<endl;
      if((_subDetectorVec[i]->GetZMin()<zPosition) && 
	 (_subDetectorVec[i]->GetZMax()>zPosition))
	{
	  ret_val=_subDetectorVec[i];
	  //cout<<"Found subDetector for "<<zPosition<<endl;
	  break;
	}
    }

  return ret_val;
}

//*************************************************************
void Detector::Reset() {
//*************************************************************

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      _subDetectorVec[i]->GetPlanes()->clear();
    }
}

//*************************************************************
int Detector::GetNPlanes() {
//*************************************************************

  int retVal = 0;

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      retVal+=_subDetectorVec[i]->GetPlanes()->size();
    }
  return retVal;
}


//*************************************************************
double Detector::CalculateChargeMomentum(vector<double>& debug) {
//*************************************************************

// fill prev de_dx for all subdetectors.
// combine all of the guessess.

  //cout<<"hits in detector="<<GetNHits()<<endl;


  double retVal = 0;

  //bool used = false;

  sort( _subDetectorVec.begin(), _subDetectorVec.end(), sortDetectorsByZ() );
  
  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      // Fill de_dx properties and do sorting of planes.
      _subDetectorVec[i]->ClearEstimate();
      sort( _subDetectorVec[i]->GetPlanes()->begin(), _subDetectorVec[i]->GetPlanes()->end(), sortPlanesByZ() );
      if(i!=0)
	{
	  /*
	  _subDetectorVec[i]->SetPrevde_dx(_subDetectorVec[i-1]->GetLength()*
					   _subDetectorVec[i-1]->Getde_dx()+
					   _subDetectorVec[i]->GetLength()/2.0*
					   _subDetectorVec[i]->Getde_dx()+
					   _subDetectorVec[i-1]->GetPrevde_dx());
	  */
	  _subDetectorVec[i]->SetPrevde_dx(_subDetectorVec[i-1]->GetLength()*
					   _subDetectorVec[i-1]->Getde_dx()+
					   _subDetectorVec[i-1]->GetPrevde_dx());
	}
    }
  /*
    for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
    
    cout<<_subDetectorVec[i]->GetName()<<endl;
    cout<<_subDetectorVec[i]->Getde_dx()<<endl;
    cout<<_subDetectorVec[i]->GetLength()<<endl;
    cout<<_subDetectorVec[i]->GetPrevde_dx()<<endl;
    }
  */
  std::vector<cluster*> stackHits;

  std::vector<cluster*> hits;

  int nPlanes = 0;

  int leverArmFirst = 0;

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      for(unsigned int j=0;j<_subDetectorVec[i]->GetPlanes()->size();j++)
	{
	  if(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size()) nPlanes++;
	  for(unsigned int k=0;k<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size();k++)
	    {
	      hits.push_back(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
	    }
	}
    }


  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {

      //cout<<"name="<<_subDetectorVec[i]->GetName()<<endl;
      // Count number of hits in the stack to see if we can use helix.
      if(_subDetectorVec[i]->GetName() == "SFFFS0" ||
	 _subDetectorVec[i]->GetName() == "S0" ||
	 _subDetectorVec[i]->GetName() == "S1" ||
	 _subDetectorVec[i]->GetName() == "SFS" ||
	 _subDetectorVec[i]->GetName() == "SFFS" ||
	 _subDetectorVec[i]->GetName() == "SFFFS" ||
	 _subDetectorVec[i]->GetName() == "SFS0" ||
	 _subDetectorVec[i]->GetName() == "SFFS0" ||
	 //_subDetectorVec[i]->GetName() == "SFFFFS0" ||
	 _subDetectorVec[i]->GetName() == "SFFFFS1" ||
	 _subDetectorVec[i]->GetName() == "SFFS1" ||
	 _subDetectorVec[i]->GetName() == "SFFFS1"
	 )
	{
	  if(i+1 < _subDetectorVec.size() && i == 0)
	    {
	      leverArmFirst = i;
	    }

	  continue;
	}
      else if(_subDetectorVec[i]->GetName() == "TASD" 
	      || _subDetectorVec[i]->GetName() == "TASD0" 
	      || _subDetectorVec[i]->GetName() == "TASD1" 
	      || _subDetectorVec[i]->GetName() == "SF"
	      || _subDetectorVec[i]->GetName() == "SF0"
	      || _subDetectorVec[i]->GetName() == "SF1"
	      || _subDetectorVec[i]->GetName() == "SF2"
	      || _subDetectorVec[i]->GetName() == "SF3")
	{
	  continue;
	}
      else
	{
	  
	  for(unsigned int j=0;j<_subDetectorVec[i]->GetPlanes()->size();j++)
	    {
	      if(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size()) nPlanes++;
	      for(unsigned int k=0;k<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size();k++)
		{
		  stackHits.push_back(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
		}
	    }
	}
      
    }
  
  //cout<<"helixHits.size()="<<helixHits.size()<<endl;
  // If enough hits, use only helix and leverarm.
  // else? Lever arm and quadratic in one region?
  //if(helixHits>10)
  
  //double retVal;
  
  //vector<double> debug2;

  vector<double> temp;

  if(stackHits.size()>5)
    {
      //double helixVal =Helix(helixHits);
      retVal = Helix(stackHits);
      //cout<<"Helix returns="<<helixVal<<endl;
      //if(helixVal==0 || fabs(helixVal) >8000 || isnan(helixVal)) helixVal = 4000;
      
      debug.push_back(retVal);
      debug.push_back(Quadratic(stackHits)[0]);
      temp = LeverArm(leverArmFirst+1);
      debug.push_back(temp[0]);
      debug.push_back(temp[1]);
      temp.clear();
      temp = LeverArm(leverArmFirst+2);
      debug.push_back(temp[0]);
      debug.push_back(temp[1]);

      cout<<"Helix  used"<<endl;

      //return helixVal;
    }
  else
    {
      double posCharge = 1;
      double negCharge = 1;

      //cout<<"LeverArm"<<endl;
      debug.push_back(0);
      debug.push_back(0);

      //cout<<"Call="<<LeverArm(leverArmFirst+1)[0]<<endl;
      temp.clear();

      temp = LeverArm(leverArmFirst+2);

      retVal = temp[0];

      if(fabs(retVal)>1)
	{
	  posCharge *= temp[2];
	  negCharge *= temp[3];
	}

      //cout<<"Should be good here"<<retVal<<endl;

      debug.push_back(temp[0]);
      debug.push_back(temp[1]);


      temp.clear();
      temp = LeverArm(leverArmFirst+1);
      //cout<<"Second call"<<retVal<<endl;
      debug.push_back(temp[0]);
      debug.push_back(temp[1]);

      if(fabs(retVal) <10 || fabs(retVal) >1000 || isnan(retVal) ) retVal = temp[0];

      //cout<<"hits.size()="<<hits.size()<<endl;

      if((fabs(retVal) <10 || fabs(retVal) >1000 || isnan(retVal) ) && hits.size() > 3)
	{
	  retVal = Quadratic(hits)[0];
	  //cout<<"Lever arm does not work"<<endl;
	  //cout<<"Quadratic used"<<endl;
	}
      //else
      //{
	  //cout<<"Lever arm used"<<endl;
      //}




      if(fabs(retVal)>1 && temp[0])
	{
	  posCharge *= temp[2];
	  negCharge *= temp[3];
	}

      //cout<<"posCharge="<<posCharge<<endl;
      //cout<<"negCharge="<<negCharge<<endl;

      if(posCharge>negCharge) retVal = fabs(retVal);
      else retVal = -1 * fabs(retVal);


    }

      //retVal = LeverArmQuadratic(debug2);
      //cout<<"Debug size="<<debug2.size()<<endl;
    
  
  cout<<"in detector.cpp="<<retVal<<endl;

  //if(fabs(retVal)<10 || fabs(retVal) >8000 || isnan(retVal)) retVal = 2000;

  if(isnan(retVal)) retVal = -2000;
  
  //retVal = 1000;

  //cout<<"in detector.cpp="<<retVal<<endl;
  /*
  if(hits.size()>0)
    {
      cout<<"track length="<<fabs(fabs(hits[hits.size()-1]->vector()[2])-fabs(hits[0]->vector()[2]))<<endl;
      cout<<"de_dx="<<GetSubDetector(hits[hits.size()-1]->vector()[2])->GetPrevde_dx()<<endl;
    }
  else
    {
      cout<<"size0_in_detector"<<endl;
    }
  */
  return retVal;
}


//*************************************************************
int Detector::GetNHits() {
//*************************************************************

  int ret_val =0;

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      ret_val+=_subDetectorVec[i]->GetNHits();
    }
  return ret_val;
}


//*************************************************************
double Detector::CalculateCharge() {
//*************************************************************

/*
Check if first planes, then leaver-arm.
Else in main detector or end.
If in main detector need to find subdetector and next plane in detector.

Sort subdetectors by z position and planes. Have function to find next subdetector and next plane.
*/

  //CalculateChargeMomentum();

  std::vector<cluster*> hits;

  //cout<<"hits in detector="<<GetNHits()<<endl;


  sort( _subDetectorVec.begin(), _subDetectorVec.end(), sortDetectorsByZ() );

  

  double retVal = 0;
  bool checkNext = false;
  unsigned int detectorNum = 0;
  /*
  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      cout<<"Detector="<<_subDetectorVec[i]->GetName()<<endl;
    }
  */

  //for(unsigned int i=0;i<_subDetectorVec.size();i++)
  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      sort( _subDetectorVec[i]->GetPlanes()->begin(), _subDetectorVec[i]->GetPlanes()->end(), sortPlanesByZ() );

      if(_subDetectorVec[i]->GetName() == "S0") continue;
      if(_subDetectorVec[i]->GetName() == "SFFFS0") continue;
      if(_subDetectorVec[i]->GetName() == "SFS") continue;
      if(_subDetectorVec[i]->GetName() == "SFFFFS1") continue;
      if(_subDetectorVec[i]->GetName() == "S1") continue;

      //cout<<"Detector="<<_subDetectorVec[i]->GetName()<<endl;


      for(unsigned int j=0;j<_subDetectorVec[i]->GetPlanes()->size();j++)
	{
	  for(unsigned int k=0;k<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size();k++)
	    {
	      //cout<<"Z="<<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]->position()[2]<<" Y="<<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]->position()[1]<<endl;
	      hits.push_back(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
	    }
	}
    }
  //cout<<"Nhits="<<hits.size()<<endl;
  if(hits.size())
    {
      retVal = ChargeOne(hits,_subDetectorVec[0]);
    }
  /*
  else if(_subDetectorVec[0]->GetNHits()>0 && _subDetectorVec[1]->GetNHits()>0 && _subDetectorVec[2]->GetNHits()>0 && _subDetectorVec[3]->GetNHits()>0)
    {
      double theta1 =  atan((_subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[1] -_subDetectorVec[0]->GetPlanes()->at(0)->GetHits()[0]->position()[1])/
			    (_subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[2] - _subDetectorVec[0]->GetPlanes()->at(0)->GetHits()[0]->position()[2]));

      double theta2 =  atan((_subDetectorVec[1]->GetPlanes()->at(1)->GetHits()[0]->position()[1] -_subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[1])/
			    (_subDetectorVec[1]->GetPlanes()->at(1)->GetHits()[0]->position()[2] - _subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[2]));

      double theta3 =  atan((_subDetectorVec[2]->GetPlanes()->at(0)->GetHits()[0]->position()[1] -_subDetectorVec[1]->GetPlanes()->at(1)->GetHits()[0]->position()[1])/
			    (_subDetectorVec[2]->GetPlanes()->at(0)->GetHits()[0]->position()[2] - _subDetectorVec[1]->GetPlanes()->at(1)->GetHits()[0]->position()[2]));

      retVal = (theta3-theta2)/fabs(theta3-theta2) *(theta2-theta1)/fabs(theta2-theta1);
      
    }
  else if(_subDetectorVec[0]->GetNHits()>0 && _subDetectorVec[1]->GetNHits()>0 && _subDetectorVec[2]->GetNHits()>0)
    {
      double theta1 =  atan((_subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[1] -_subDetectorVec[0]->GetPlanes()->at(0)->GetHits()[0]->position()[1])/
			    (_subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[2] - _subDetectorVec[0]->GetPlanes()->at(0)->GetHits()[0]->position()[2]));

      double theta2 =  atan((_subDetectorVec[1]->GetPlanes()->at(1)->GetHits()[0]->position()[1] -_subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[1])/
			    (_subDetectorVec[1]->GetPlanes()->at(1)->GetHits()[0]->position()[2] - _subDetectorVec[1]->GetPlanes()->at(0)->GetHits()[0]->position()[2]));




      retVal=(theta2-theta1)/fabs(theta2-theta1);
    }
  */
  else
    {
      retVal = 1;
    }
      

    
  //return 1;
  //cout<<retVal <<endl;
  hits.clear();
  
  //retVal = ChargeOne(hits,_subDetectorVec[detectorNum]);
  //cout<<retVal <<endl;
  //hits.clear();


  return retVal;
}

//*************************************************************
double Detector::ChargeOne(std::vector<cluster*>& hits,SubDetector* detector) {
//*************************************************************

  //cout<<"In ChargeOne"<<endl;

  int nMeas = hits.size();
  int minindex = nMeas;  


  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas];
  
  double firstNodeZ = hits[0]->position()[2];
  
  
  pos[0] = x[0] = hits[0]->vector()[0];
  pos[1] = y[0] = hits[0]->vector()[1];
  pos[2] = z[0] = hits[0]->vector()[2];
  
  
  //EVector BfieldVector = detector->GetBField()->vector(pos);
  
  EVector prevB(3,0);
  
  for (int iMeas = 1; iMeas<nMeas ;iMeas++){
    x[iMeas] = hits[iMeas]->vector()[0];
    y[iMeas] = hits[iMeas]->vector()[1];
    z[iMeas] = hits[iMeas]->vector()[2];

    // get the b-field from the previous step
    //EVector B = detector->GetBField()->vector(pos);
    EVector B = GetSubDetector(hits[iMeas]->vector()[2])->GetBField()->vector(pos);

    //cout<<"Is b null"<<endl;

    //cout<<"FIELD="<<B[0]<<" "<<B[1]<<" "<<B[2]<<endl;
    //cout<<"Is b no"<<endl;

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

  (void) gr->Fit("parfit2", "QN");

  // Positive particle -> Positive parameters -> Qtilde +
  // Negative Particle -> Negative parameters -> Qtilde -

  double qtilde = 1000 * 0.3* prevB[0]*pow(1+pow(fun->GetParameter(1),2),3./2.)/
    (2*fun->GetParameter(2)*0.01);

  int meansign = (int)(qtilde/fabs(qtilde));

  if(meansign == 0) meansign = 1;

  delete gr;
  delete fun;


  return meansign;

}


//*************************************************************
//double Detector::Quadratic(std::vector<cluster*>& hits) {
vector<double> Detector::Quadratic(std::vector<cluster*>& hits) {
//*************************************************************

  int nMeas = hits.size();
  int minindex = nMeas;  

  double erZ = 10; //mm
  double errY = 15; //mm

  for (int iMeas = 0; iMeas<nMeas ;iMeas++){
    /*
    cout<< hits[iMeas]->vector()[0]<<" "
	<< hits[iMeas]->vector()[1]<<" "
	<< hits[iMeas]->vector()[2]<<endl;
    */
  }


  vector<double> retVec;

  EVector pos(3,0);
  EVector Z(3,0); Z[2] = 1;
  double x[(const int)nMeas], y[(const int)nMeas], 
    z[(const int)nMeas], u[(const int)nMeas], errZ[(const int)nMeas], errU[(const int)nMeas];
  
  double firstNodeZ = hits[0]->position()[2];
  
  
  pos[0] = x[0] = hits[0]->vector()[0];
  pos[1] = y[0] = hits[0]->vector()[1];
  pos[2] = z[0] = hits[0]->vector()[2];

  errZ[0] = erZ;
  
  
  //EVector BfieldVector = detector->GetBField()->vector(pos);
  
  EVector prevB(3,0);
  
  for (int iMeas = 1; iMeas<nMeas ;iMeas++){
    x[iMeas] = hits[iMeas]->vector()[0]; 
    y[iMeas] = hits[iMeas]->vector()[1];
    z[iMeas] = hits[iMeas]->vector()[2];

    errZ[iMeas] = erZ;

    // get the b-field from the previous step
    //EVector B = detector->GetBField()->vector(pos);
    EVector B = GetSubDetector(hits[iMeas]->vector()[2])->GetBField()->vector(pos);

    //cout<<"FIELD="<<B[0]<<" "<<B[1]<<" "<<B[2]<<endl;

    if(dot(prevB,B) < 0)
      {
	break;
      }
    prevB=B;

    pos[0] =  x[iMeas];  pos[1] =  y[iMeas];   pos[2] = z[iMeas];

    if(crossprod(Z,B).norm() ==0 )
      {
	u[iMeas] = 0;
      }
    else
      {
	//u[iMeas] = dot(pos, crossprod(Z,B))/crossprod(Z,B).norm();
	u[iMeas] = dot(pos, crossprod(Z,B))/B.norm();
	//u[iMeas] = pos[1];
      }
    errU[iMeas] = errY;
  }

  //TGraphErrors *gr = new TGraphErrors ((const int)minindex, z, u, errZ, errU); 

  TGraph *gr = new TGraph((const int)minindex, z, u);
  
  TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x",-300,300);
  //TF1 *fun = new TF1("parfit2","[0]+[1]*x+[2]*x*x+[3]*x*x",-300,300);
  //TF1 *fun = new TF1("parfit2","pol3",-3000,3000);
  //TF1 *fun = new TF1("parfit2","pol2",-3000,3000);
  fun->SetParameters(0.,0.001,0.001);

  //fun->SetParameters(100.0,0.1,0.00001);

  //(void) gr->Fit("parfit2", "QNEM");
  //(void) gr->Fit("parfit2", "QNEF");
  (void) gr->Fit("parfit2", "QN");
  // (void) gr->Fit("parfit2", "QN");
  /*
  cout<<"Param0="<<fun->GetParameter(0)<<endl;
  cout<<"Param1="<<fun->GetParameter(1)<<endl;
  cout<<"Param2="<<fun->GetParameter(2)<<endl;
  //cout<<"Param3="<<fun->GetParameter(3)<<endl;

  cout<<"Error0="<<fun->GetParError(0)<<endl;
  cout<<"Error1="<<fun->GetParError(1)<<endl;
  cout<<"Error2="<<fun->GetParError(2)<<endl;
  //cout<<"Error3="<<fun->GetParError(3)<<endl;
  */
  //cout<<"x2="<<fun->GetChisquare()<<endl;

  //cout<<"Px2="<<fun->GetProb()<<endl;

  // Positive particle -> Positive parameters -> Qtilde +
  // Negative Particle -> Negative parameters -> Qtilde -


  double charge = 1.0;

  if(y[nMeas-1]>y[0]) charge = -1.0;
 

  //double qtilde = charge * 1000* 0.3*prevB[0]*pow(1+pow(fun->GetParameter(1),2),3./2.)/
  //(2*fun->GetParameter(2));

  SubDetector* sub = GetSubDetector(hits[nMeas-1]->vector()[2]);  

  double fieldVal = sub->GetBField()->vector(sub->GetPlanes()->at(0)->GetHits()[0]->position())[0];

  charge*=fieldVal/fabs(fieldVal);


  //double p = (sub->GetPrevde_dx() + 
  //(sub->GetMinPlane()->GetHits()[0]->position()[2] - sub->GetZMin())*sub->Getde_dx());

  double p = RangeMomentum(sub);
  
  double qtilde = charge * p;


  //qtilde = fabs(qtilde);

  //qtilde *=- fun->GetParameter(1)/fabs(fun->GetParameter(1));

  //double eqtilde = 1000* 0.3*prevB[0]*pow(1+pow(fabs(fun->GetParameter(1))+fun->GetParError(1),2),3./2.)/
  //(2*(fabs(fun->GetParameter(2))-fun->GetParError(2)));


  //cout<<"errorCalc="<<fabs(fabs(qtilde)-fabs(eqtilde))<<endl;


  //double qerr = sqrt(pow(3.0/2.0 * 2 * fun->GetParError(1)/ fun->GetParameter(1),2) + pow(fun->GetParError(2)/ fun->GetParameter(2),2));


  // sigma^2 = (dp/da)^2 sigma_a^2 + (dp/db)^2 sigma_b^2

  /*
  double qerr = (qtilde / fun->GetParameter(2))* (qtilde / fun->GetParameter(2)) * fun->GetParError(2) * fun->GetParError(2) +
    qtilde * 3* fun->GetParameter(1)/ (1+fun->GetParameter(1)*fun->GetParameter(1)) *  fun->GetParError(1) *
    qtilde * 3* fun->GetParameter(1)/ (1+fun->GetParameter(1)*fun->GetParameter(1)) * fun->GetParError(1);
  */  
  
  //qerr =sqrt(qerr);
  
  //cout<<"errorCalc="<<qerr<<endl;
  
  //int meansign = (int)(qtilde/fabs(qtilde));

  //if(meansign == 0) meansign = 1;

  cout<<"quadratic="<<qtilde<<endl;

  //return qtilde;

  //if(qtilde==0 || fabs(qtilde) >8000 || isnan(qtilde)) qtilde = 400;

  //if(isnan(qtilde)) qtilde = 0;

  //if(qtilde != 0) qtilde = qtilde/fabs(qtilde) * (fabs(qtilde) + GetSubDetector(hits[0]->vector()[2])->Getde_dx());

  retVec.push_back(qtilde);
  //retVec.push_back(qerr);

  delete gr;

  delete fun;

  return retVec;

}

//*************************************************************
vector<double> Detector::ThreePointCircle(std::vector<cluster*>& hits) {
//*************************************************************

//Using for instance the theory from:
//http://paulbourke.net/geometry/circlesphere/

  double z0 = hits[0]->vector()[2];
  double z1 = hits[1]->vector()[2];
  double z2 = hits[2]->vector()[2];

  double y0 = hits[0]->vector()[1];
  double y1 = hits[1]->vector()[1];
  double y2 = hits[2]->vector()[1];


  vector<double> retVec;

  double ma = (y1-y0)/
    (z1-z0);

  double mb = (y2-y1)/
    (z2-z1);

  double zc = (ma*mb*(y0-y2)+mb*(z0+z1)-ma*(z1+z2))/(2*(mb-ma));

  double yc = -1/mb * (zc -(z1+z2)/2) + (y1+y2)/2;


  //double yc = -1/mb *(zc-(hits[1]->vector()[2]+hits[2]->vector()[2])/2) +
  //(hits[1]->vector()[2]+hits[2]->vector()[1])/2;

  double z = z1-zc;

  double y = y1-yc;

  double r = sqrt(z*z+y*y);

  EVector B = GetSubDetector(z1)->GetBField()->vector(hits[1]->vector());

  int charge = 0;

  if(y>0)
    {
      charge = +1;
    }
  else
    {
      charge = -1;
    }
  

  double p = charge* 0.3* B[0] * r * 1000;

  // error calc.
  // read from saroman.
  double errY = 15; //mm
  double errZ = 10; //mm

  double errMa = ma*ma * (errY*errY+errY*errY)/
    ((y1-y0)*(y1-y0)) +
    ma*ma * (errZ*errZ+errZ*errZ)/
    ((z1-z0)*(z1-z0));

  errMa = sqrt(errMa);
  
 double errMb = mb*mb * (errY*errY+errY*errY)/
    ((y2-y1)*(y2-y1)) +
    ma*ma * (errZ*errZ+errZ*errZ)/
    ((z2-z1)*(z2-z1));

  errMb = sqrt(errMb);
  
  double mba= mb-ma;

  double ema = -zc/(mba) - (zc-mb*(z0+z1))/ma;

  double emb = zc/(mba) + (zc+ma*(z1+z2))/mb;

  double ey2 = -ma*mb/(2*mba);
  double ey0 = ma*mb/(2*mba);
  double ez2 = -ma/(2*mba);
  double ez1 = 1.0/2.0;
  double ez0 = mb/(2*mba);

  double ezc = ema*ema*errMa*errMa+emb*emb*errMb*errMb+
    (ey2*ey2+ey0*ey0)*errY*errY +
    (ez2*ez2+ez1*ez1+ez0*ez0)*errZ*errZ;

  ezc = sqrt(ezc);

  emb = (yc-(y1+y2)/2)/(-mb);

  double eyc = emb*emb*mb*mb + zc*zc /(mb*mb) + (z1*z1+z2*z2)/(2*mb*2*mb) + (y1*y1+y2*y2)/(2*2+2*2);

  eyc = sqrt(eyc);
  

  double errR = z*z/(r*r)*(errZ*errZ+ezc*ezc)+y*y/(r*r)*(errY*errY+eyc*eyc);
  errR = sqrt(errR);
  
  //double p = charge* 0.3* B[0] * r * 1000;

  double err = charge*0.3*1000*B[0] *errR;

  err=fabs(err);

  //cout<<"threepoint="<<p<<endl;
  //cout<<err<<endl;

  retVec.push_back(p);
  retVec.push_back(err);

  //return p;
  return retVec;

}


//*****************************************************************************
double fitfunc(Double_t *x,Double_t *par) { 
  //*****************************************************************************

  double z = x[0];

  double fitval = par[0] + par[1]*z+par[2]*z*z+par[3]*z*z*z+par[4]*z*z*z*z;

  return fitval;
}


//*****************************************************************************
double Detector::Helix(std::vector<cluster*>& hits){
  //*****************************************************************************
  /*
  cout<<"in Helix"<<endl;

  for(unsigned int cnt = 0; cnt<hits.size();cnt++) 
    {
      cout<< hits[cnt]->vector()[0]<<" "
	  <<hits[cnt]->vector()[1]<<" "
	  <<hits[cnt]->vector()[2]<<endl;
    }
  
  */

  //Some catchers for pointless returns.
  int fitcatch;
  //
  int nfit, sign;
  int fitRange[3];
  int nplanes = hits.size()-1;
  int firsthit = 0;
  const int fitpoints =  hits.size()-1;//nplanes - firsthit;
  
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

    xpos[pos] = hits[ipoint]->vector()[0];
    ypos[pos] = hits[ipoint]->vector()[1];
    zpos[pos] = hits[ipoint]->vector()[2]
      - hits[firsthit]->vector()[2];
    currentpos[0] = hits[ipoint]->vector()[0];
    currentpos[1] = hits[ipoint]->vector()[1];
    currentpos[2] = hits[ipoint]->vector()[2];
    currentB = GetSubDetector(hits[ipoint]->vector()[2])->GetBField()->vector(currentpos);
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
    
    TF1 *func = new TF1("fit",fitfunc,-3,3,3);
    func->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func->SetParNames("a", "b", "c", "d", "e");
    
    TF1 *func2 = new TF1("fit2",fitfunc,-3,3,3);
    func2->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func2->SetParNames("f", "g", "h", "i", "j");

    TF1 *func3 = new TF1("fit3",fitfunc,-3,3,3);
    func3->SetParameters(0.,0.,0.001,0.0001,0.0001);
    func3->SetParNames("a1", "b1", "c1", "d1", "e1");
    
    TF1 *func4 = new TF1("fit4",fitfunc,-3,3,3);
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
	V[5] = 1./(0.3*Bmean*pow((1+g1*g1),3./2.)/
		   (2*h1)*0.001);
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


  double retVal = 1/V[5];

  //if(retVal==0 || fabs(retVal) >8000 || isnan(retVal)) retVal = 400;

  if(isnan(retVal)) retVal = 0;

  if(retVal > 8000) retVal = 8000;
  
  if(retVal < -8000) retVal = -8000;

  if(retVal != 0) retVal = retVal/fabs(retVal) * (fabs(retVal) + GetSubDetector(hits[0]->vector()[2])->Getde_dx());
  
  return retVal;

  //std::cout<<"Momentum guess from polynomial fit: p/q = "<<1./V[5]<<std::endl;
}

//*****************************************************************************
double Detector::LeverArmQuadratic(vector<double>& debug){
//*****************************************************************************


  double retVal = 0;

  bool used = false;

  for(unsigned int i=1;i<_subDetectorVec.size()-2;i++)
    {
      //cout<<i<<endl;
      if(_subDetectorVec[i]->GetName() == "SFFFS0" ||
	 _subDetectorVec[i]->GetName() == "SFS" ||
	 _subDetectorVec[i]->GetName() == "SFFS" ||
	 _subDetectorVec[i]->GetName() == "SFFFS" ||
	 _subDetectorVec[i]->GetName() == "SFS0" ||
	 _subDetectorVec[i]->GetName() == "SFFS0" ||
	 //_subDetectorVec[i]->GetName() == "SFFFFS0" ||
	 _subDetectorVec[i]->GetName() == "SFFFFS1" ||
	 _subDetectorVec[i]->GetName() == "SFFS1" ||
	 _subDetectorVec[i]->GetName() == "SFFFS1"
	 )
	{
	  
	  if( _subDetectorVec[i-1]->GetMaxPlane() &&_subDetectorVec[i]->GetMinPlane() && 
	      _subDetectorVec[i]->GetMaxPlane() &&  _subDetectorVec[i+1]->GetMinPlane())
	    {
	      if(_subDetectorVec[i-1]->GetMaxPlane()->GetHits().size()>0 &&
		 _subDetectorVec[i]->GetMinPlane()->GetHits().size()>0 &&
		 _subDetectorVec[i]->GetMinPlane() !=  _subDetectorVec[i]->GetMaxPlane() &&
		 _subDetectorVec[i]->GetMaxPlane()->GetHits().size()>0 &&
		 _subDetectorVec[i+1]->GetMinPlane()->GetHits().size()>0)
		{
		  //debug.push_back(LeverArm(i)[0]);

		} // End if enough hits
	    } // End if good pointers
	  //}
	  
	}
      else // Do not use leverarm. (In stack) but can not use helix.
	{
	  //cout<<"In stack"<<endl;
	  
	  if(_subDetectorVec[i]->GetName() == "SFFFFS0" && used) continue;
	  
	  std::vector<cluster*> hits;
	  
	  if(_subDetectorVec[i]->GetName() == "SFFFFModule" &&
	     _subDetectorVec[i+1]->GetName() == "SFFFFS0")
	    {
	      for(unsigned int j=0;j<_subDetectorVec[i]->GetPlanes()->size();j++)
		{
		  for(unsigned int k=0;k<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size();k++)
		    {
		      hits.push_back(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
		    }
		}
	      
	      for(unsigned int j=0;j<_subDetectorVec[i+1]->GetPlanes()->size();j++)
		{
		  for(unsigned int k=0;k<_subDetectorVec[i+1]->GetPlanes()->at(j)->GetHits().size();k++)
		    {
		      hits.push_back(_subDetectorVec[i+1]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
		    }
		}
	      used = true;
	    }
	  else
	    {	
	      
	      for(unsigned int j=0;j<_subDetectorVec[i]->GetPlanes()->size();j++)
		{
		  for(unsigned int k=0;k<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size();k++)
		    {
		      hits.push_back(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
		    }
		}
	      if(_subDetectorVec[i+1]->GetPlanes()->size())
		{
		  if(_subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits().size())
		    {
		      hits.push_back(_subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits()[0]);
		    }
		}
	    }		
	  
	  //cout<<_subDetectorVec[i]->GetName()<<" "<<hits.size()<<endl;
	  
	  if(hits.size()>2)
	    {
	      if(hits.size()>3 && _subDetectorVec[i]->GetNPlanes()> 2)
		{
		  vector<double> p;
		  //cout<<"Before quadratic "<<_subDetectorVec[i]->GetName()<<endl;
		  p = Quadratic(hits);
		  
		  _subDetectorVec[i]->SetCharge(p[0]/fabs(p[0]));
		  _subDetectorVec[i]->SetMomentum(fabs(p[0]));
		  //_subDetectorVec[i]->SetChargeProb(chargeProb);
		  _subDetectorVec[i]->SetError(p[1]);
		  
		  // Set to crazy value, needs to be calculated after mean p;
		  
		  _subDetectorVec[i]->SetChargeProb(1.0);		  
		  _subDetectorVec[i]->SetNegChargeProb(1.0);
		  
		}
	      else
		{
		  /*
		    vector<double> p;
		    cout<<"Using threePoint"<<endl;
		    p = ThreePointCircle(hits);
		    
		    _subDetectorVec[i]->SetCharge(p[0]/fabs(p[0]));
		    _subDetectorVec[i]->SetMomentum(fabs(p[0]));
		    _subDetectorVec[i]->SetError(p[1]);
		    cout<<"Returns="<<p[0]<<" "<<p[1]<<endl;
		    
		    // Set to crazy value, needs to be calculated after mean p;
		    _subDetectorVec[i]->SetChargeProb(1.0);
		    _subDetectorVec[i]->SetNegChargeProb(1.0);
		  */
		}
	    }
	  //else
	  
	}
      
      
    }//end for.

  //Combine the guesses.
  
  double tempP =0;
  double tempErr =0;
  
  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      //cout<<"subDetector="<<_subDetectorVec[i]->GetMomentum()<<endl;
      
      
      if(_subDetectorVec[i]->GetMomentum() == 0 || isinf(_subDetectorVec[i]->GetMomentum()) ||
	 isnan(_subDetectorVec[i]->GetMomentum()) || isnan(_subDetectorVec[i]->GetError()) ||
	 _subDetectorVec[i]->GetError() == 0 || isinf(_subDetectorVec[i]->GetError()))
	{
	  continue;
	}
      else
	{
	  /*
	    cout<<"Using detector="<<_subDetectorVec[i]->GetName()<<endl;
	    cout<<_subDetectorVec[i]->GetMomentum()<<endl;
	    cout<<_subDetectorVec[i]->GetError()<<endl;
	    cout<<_subDetectorVec[i]->GetChargeProb()<<endl;
	    cout<< _subDetectorVec[i]->GetNegChargeProb()<<endl;
	  */
	  tempP += (_subDetectorVec[i]->GetMomentum()+ _subDetectorVec[i]->GetPrevde_dx())
	    /(_subDetectorVec[i]->GetError() *_subDetectorVec[i]->GetError() );
	  //+ _subDetectorVec[i]->GetPrevde_dx();
	  
	  tempErr += 1/(_subDetectorVec[i]->GetError()*_subDetectorVec[i]->GetError());
	}
      
    }

//cout<<"TempP="<<tempP<<endl;
//cout<<"TempErr="<<tempErr<<endl;

 tempP=tempP/tempErr;

 //cout<<"Final="<<tempP<<endl;


 double chargeProb = 1;
 double negChargeProb = 1;
 for(unsigned int i=0;i<_subDetectorVec.size();i++)
   {
     //cout<<"chargeProb="<<chargeProb<<endl;
     //cout<<"NchargeProb="<<negChargeProb<<endl;
     
     if(_subDetectorVec[i]->GetMomentum() == 0 || isinf(_subDetectorVec[i]->GetMomentum()) ||
	isnan(_subDetectorVec[i]->GetMomentum()) || isnan(_subDetectorVec[i]->GetError()) ||
	_subDetectorVec[i]->GetError() == 0 || isinf(_subDetectorVec[i]->GetError()))
       {
	 continue;
       }
     else
       {
	 double tempChargeProb;
	 double tempNegChargeProb;
	 
	 if(_subDetectorVec[i]->GetChargeProb() == 1 &&  _subDetectorVec[i]->GetNegChargeProb() == 1)
	   {
	     // Need to calculate prob from meanP;
	     
	     double subP = _subDetectorVec[i]->GetCharge() *_subDetectorVec[i]->GetMomentum();
	     //double subC = _subDetectorVec[i]->GetCharge();
	     double subE = _subDetectorVec[i]->GetError();
	     
	     
	     tempChargeProb=  1/(sqrt(2* M_PI) * fabs(subE)) * exp(-(subP-tempP)*(subP-tempP)/(2*subE*subE));
	     tempNegChargeProb= 1/(sqrt(2* M_PI) * fabs(subE)) * exp(-(-subP-tempP)*(-subP-tempP)/(2*subE*subE));
	     
	     if(!tempChargeProb || !tempNegChargeProb) 
	       {
		 //cout<<"Prob is 0 ?! Mom="<<subP<<"+-"<<subE<<endl;
		 continue;
	       }
	     else
	       {
		 chargeProb *= tempChargeProb;
		 //chargeProb *= tempNegChargeProb;
		 negChargeProb *= tempNegChargeProb;
		 //negChargeProb *= tempChargeProb;
	       }
	     
	   }
	 else
	   {
	     tempChargeProb =_subDetectorVec[i]->GetChargeProb();
	     tempNegChargeProb = _subDetectorVec[i]->GetNegChargeProb();
	     
	     chargeProb *= tempChargeProb;
	     negChargeProb *= tempNegChargeProb;
	     //chargeProb*=_subDetectorVec[i]->GetChargeProb();
	     //negChargeProb*= _subDetectorVec[i]->GetNegChargeProb();
	   }
	 /*	
	   cout<<_subDetectorVec[i]->GetName()<<" "
	   <<_subDetectorVec[i]->GetMomentum()*_subDetectorVec[i]->GetCharge()<<" "
	   <<_subDetectorVec[i]->GetMomentum()+_subDetectorVec[i]->GetPrevde_dx()<<" "
	   <<_subDetectorVec[i]->GetError()<<" "
	   <<tempChargeProb<<" "
	   <<tempNegChargeProb<<endl;
	 */
       }
     
   }//end for

 //cout<<"final chargeProb="<<chargeProb<<endl;
 //cout<<"final NchargeProb="<<negChargeProb<<endl;
 
 if(chargeProb<negChargeProb)
   {
     retVal = -1* tempP;
   }
 else
   {
     retVal = tempP;
   }
 
 
 
 if(retVal==0 || fabs(retVal) >8000 || isnan(retVal)) retVal = 4000;
 
 
 //cout<<"retVal="<<retVal<<endl;
 
 //retVal = 4500;
 
 return retVal;
}

//*****************************************************************************
vector<double> Detector::LeverArm(unsigned int i){
//*****************************************************************************

  double retVal;
  double angle = 0;
  double negChargeProb =0;
  double chargeProb = 0;
  double badGuess = 1000;
  vector<double> retValVec;
  /*
  cout<<_subDetectorVec[i-1]->GetName()<<endl;
  cout<<_subDetectorVec[i]->GetName()<<endl;
  cout<<_subDetectorVec[i+1]->GetName()<<endl;
  */
 
  //cout<<_subDetectorVec[i-1]->GetMaxPlane()<<endl;
  //cout<<_subDetectorVec[i]->GetMinPlane()<<endl;
  //cout<<_subDetectorVec[i]->GetMaxPlane()<<endl;
  //cout<<_subDetectorVec[i+1]->GetMinPlane()<<endl;

  
  if( _subDetectorVec[i-1]->GetMaxPlane() &&_subDetectorVec[i]->GetMinPlane() && 
      _subDetectorVec[i]->GetMaxPlane() &&  _subDetectorVec[i+1]->GetMinPlane())
    {

      //cout<<(_subDetectorVec[i-1]->GetMaxPlane()->GetHits().size()>0)<<endl;
      //cout<<(_subDetectorVec[i]->GetMinPlane()->GetHits().size()>0)<<endl;
      //cout<<(_subDetectorVec[i]->GetMinPlane() !=  _subDetectorVec[i]->GetMaxPlane())<<endl;
      //cout<<(_subDetectorVec[i]->GetMaxPlane()->GetHits().size()>0)<<endl;
      //cout<<(_subDetectorVec[i+1]->GetMinPlane()->GetHits().size()>0)<<endl;

      
      if(_subDetectorVec[i-1]->GetMaxPlane()->GetHits().size()>0 &&
	 _subDetectorVec[i]->GetMinPlane()->GetHits().size()>0 &&
	 _subDetectorVec[i]->GetMinPlane() !=  _subDetectorVec[i]->GetMaxPlane() &&
	 _subDetectorVec[i]->GetMaxPlane()->GetHits().size()>0 &&
	 _subDetectorVec[i+1]->GetMinPlane()->GetHits().size()>0)
	{


  double dy = (_subDetectorVec[i]->GetMinPlane()->GetHits()[0]->position()[1] -
	       _subDetectorVec[i-1]->GetMaxPlane()->GetHits()[0]->position()[1]);
  double dz = (_subDetectorVec[i]->GetMinPlane()->GetHits()[0]->position()[2] - 
	       _subDetectorVec[i-1]->GetMaxPlane()->GetHits()[0]->position()[2]);
  
  double theta1 = atan(dy/dz);
  
  //double errZ = 10; //mm
  double errZ = 0; //mm
  double errY = 15; //mm
  
  //double errTheta1 = 1.0/(1.0+(dy/dz * dy/dz)) * 1.0/(1.0+(dy/dz * dy/dz)) * (errY*errY+errY*errY/(errZ*errZ*errZ*errZ));//(2*errZ*errZ+2*errY*errY);
  
  //double errTheta1 = 1.0/(1.0+(dy/dz * dy/dz)) * 1.0/(1.0+(dy/dz * dy/dz)) * errY*errY /(dz*dz) + errZ*errZ*dy*dy/(dz*dz*dz*dz);  
  
  double errTheta1 = 1.0/(dy*dy+dz*dz) * 1.0/(dy*dy+dz*dz) * (dz*dz*errY*errY + dy*dy*errZ*errZ);
  
  errTheta1 = sqrt(errTheta1);
  
  dy = (_subDetectorVec[i+1]->GetMinPlane()->GetHits()[0]->position()[1] -_subDetectorVec[i]->GetMaxPlane()->GetHits()[0]->position()[1]);
  
  dz = (_subDetectorVec[i+1]->GetMinPlane()->GetHits()[0]->position()[2] - _subDetectorVec[i]->GetMaxPlane()->GetHits()[0]->position()[2]);
  
  double theta2 =atan(dy/dz);
  
  //double errTheta2 = 1.0/(1.0+(dy/dz * dy/dz)) * 1.0/(1.0+(dy/dz * dy/dz)) * (errY*errY+errY*errY/(errZ*errZ*errZ*errZ));
  
  //double errTheta2 = 1.0/(1.0+(dy/dz * dy/dz)) * 1.0/(1.0+(dy/dz * dy/dz)) * errY*errY /(dz*dz) + errZ*errZ*dy*dy/(dz*dz*dz*dz);  
  
  double errTheta2 = 1.0/(dy*dy+dz*dz) * 1.0/(dy*dy+dz*dz) * (dz*dz*errY*errY + dy*dy*errZ*errZ);
  
  errTheta2 = sqrt(errTheta2);
  
  double delta1 = theta2-theta1;

  angle = delta1;
  
  //double delta1 = theta1-theta2;
  
  double errDelta1 = errTheta2*errTheta2 + errTheta1*errTheta1;
  
  errDelta1 = sqrt(errDelta1);
  
  //cout<<"before l"<<endl;
  double l = (_subDetectorVec[i]->GetPlanes()->at(1)->GetHits()[0]->position()[2]  - 
	      _subDetectorVec[i]->GetPlanes()->at(0)->GetHits()[0]->position()[2])*1000;
  //cout<<"before p"<<endl;
  //double p = 0.3* 
  //_subDetectorVec[i]->GetBField()->vector(_subDetectorVec[i]->GetPlanes()->at(0)->GetHits()[0]->position())[0] *
  //l/ delta1;
  /*
  double p = 0.3* 
    _subDetectorVec[i]->GetBField()->vector(_subDetectorVec[i]->GetMinPlane()->GetHits()[0]->position())[0]*
    l/ delta1;
  */

  double p;
  double fieldVal = _subDetectorVec[i]->GetBField()->vector(_subDetectorVec[i]->GetPlanes()->at(0)->GetHits()[0]->position())[0];

  if(delta1)
    {
      //p= fieldVal/fabs(fieldVal) *delta1/fabs(delta1) * (_subDetectorVec[i+1]->GetPrevde_dx() + 
      //			(_subDetectorVec[i+1]->GetMinPlane()->GetHits()[0]->position()[2] - _subDetectorVec[i+1]->GetZMin())*_subDetectorVec[i+1]->Getde_dx());
    
      p= fieldVal/fabs(fieldVal) *delta1/fabs(delta1) * RangeMomentum(_subDetectorVec[i+1]);
}
  else
    {
      //p = badGuess * fieldVal/fabs(fieldVal);
      p=0;
    }

  
  //cout<<"Momentum="<<p<<endl;
  
  _subDetectorVec[i]->SetCharge(p/fabs(p));
  _subDetectorVec[i]->SetMomentum(fabs(p));
  _subDetectorVec[i]->SetError(fabs(p*errDelta1/delta1));
  
  double X0 = _subDetectorVec[i]->GetX0Eff();
  
  double chargeWidth = 13.6/(fabs(p)*(fabs(p)/sqrt(p*p+105.65*105.65))) *
    sqrt(l/X0) *
    (1+0.0036*log(l/X0));
  
  chargeWidth = chargeWidth*chargeWidth + errDelta1*errDelta1;
  chargeWidth= sqrt(chargeWidth);
  
  //double angularResolution = 1.0/
  
  // B field negative, changes places of the probabilities.
		  	  
  //double chargeProb =  1/(sqrt(2* M_PI) * chargeWidth) * exp(-(delta1-chargeWidth)*(delta1-chargeWidth)/(2*chargeWidth*chargeWidth));
  
  //double negChargeProb =  1/(sqrt(2* M_PI) * chargeWidth) * exp(-(-delta1-chargeWidth)*(-delta1-chargeWidth)/(2*chargeWidth*chargeWidth));
  
  negChargeProb =  1/(sqrt(2* M_PI) * chargeWidth) * exp(-(delta1-chargeWidth)*(delta1-chargeWidth)/(2*chargeWidth*chargeWidth));
  
  chargeProb =  1/(sqrt(2* M_PI) * chargeWidth) * exp(-(-delta1-chargeWidth)*(-delta1-chargeWidth)/(2*chargeWidth*chargeWidth));
  
  _subDetectorVec[i]->SetChargeProb(chargeProb);
  
  _subDetectorVec[i]->SetNegChargeProb(negChargeProb);
    
  retVal = p;

  //retVal = delta1;
	}
      else
	{
	  cout<<"Lever arm fail 1"<<endl;
	}
    }
  else
    {
      cout<<"Lever arm fail 2"<<endl;
    }

  //if(retVal==0 || fabs(retVal) >8000 || isnan(retVal)) retVal = 400;
  if(isnan(retVal)) retVal = 0;

  //if(retVal==0) retVal = badGuess;

  //if(retVal==0) retVal = badGuess;

  //if(retVal != 0) retVal = retVal/fabs(retVal) * (fabs(retVal) + _subDetectorVec[i]->Getde_dx());

  retValVec.push_back(retVal);
  retValVec.push_back(angle);
  retValVec.push_back(chargeProb);
  retValVec.push_back(negChargeProb);



  cout<<"lever_arm function="<<retVal<<endl;

  return retValVec;
}

//*****************************************************************************
double Detector::RangeMomentum(SubDetector* subDet){
//*****************************************************************************
// Calculates the momentum through a range calculation.
// Requires the subdetector object for the last (highest Z) hit.
// Basically using a continuous slowing-down approximation (CSDA)


  double ke = subDet->GetPrevde_dx() + 
    (fabs(subDet->GetMinPlane()->GetHits()[0]->position()[2]) - 
     fabs(subDet->GetZMin()))*subDet->Getde_dx(); // Kinetic energy in MeV.

  double mass = 105.658; //mev/c^2 // taking muon mass

  double p = 1.5 * sqrt(ke*(ke+2.0*mass));
  
    
  return p;
}

