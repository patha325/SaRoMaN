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
//double Detector::CalculateCharge() {
//*************************************************************

/*
Check if first planes, then leaver-arm.
Else in main detector or end.
If in main detector need to find subdetector and next plane in detector.

Sort subdetectors by z position and planes. Have function to find next subdetector and next plane.
*/
/*
  std::vector<cluster*> hits;


  sort( _subDetectorVec.begin(), _subDetectorVec.end(), sortDetectorsByZ() );

  

  double retVal = 0;
  bool checkNext = false;
  unsigned int detectorNum = 0;

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      cout<<"Detector="<<_subDetectorVec[i]->GetName()<<endl;
    }
  

  //for(unsigned int i=0;i<_subDetectorVec.size();i++)
  for(unsigned int i=3;i<_subDetectorVec.size()-2;i++)
    {
      sort( _subDetectorVec[i]->GetPlanes()->begin(), _subDetectorVec[i]->GetPlanes()->end(), sortPlanesByZ() );

      //if(_subDetectorVec[i]->GetNHits()<3) continue;
      
      //if(!checkNext && _subDetectorVec[i]->GetNPlanes()<3) continue; 

      if(!checkNext && _subDetectorVec[i]->GetNPlanes()<2) continue; 
      

      //detectorNum = i -1;
      
      //cout<<"Detector="<<_subDetectorVec[i]->GetName()<<endl;
      //cout<<"Detector="<<_subDetectorVec[i]->GetwSc()<<endl;
      //cout<<"Detector="<<_subDetectorVec[i]->GetwFe()<<endl;
      //cout<<"Detector="<<_subDetectorVec[i]->GetX0Eff()<<endl;
      //cout<<"Detector="<<_subDetectorVec[i]->Getde_dx()<<endl;
      //cout<<"Detector="<<_subDetectorVec[i]->GetFieldScale()<<endl;
      //cout<<"Detector="<<_subDetectorVec[i]->GetLength()<<endl;
      //cout<<_subDetectorVec[i]->GetNHits()<<endl;
      
      for(unsigned int j=0;j<_subDetectorVec[i]->GetPlanes()->size();j++)
	{
	  for(unsigned int k=0;k<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size();k++)
	    {
	      cout<<"Z="<<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]->position()[2]<<" Y="<<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]->position()[1]<<endl;
	      hits.push_back(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
	    }

	  //if(checkNext) break;
	}

      for(unsigned int j=0;j<_subDetectorVec[i+1]->GetPlanes()->size();j++)
	{
	  //for(unsigned int k=0;k<_subDetectorVec[i+1]->GetPlanes()->at(j)->GetHits().size();k++)
	  //{
	  unsigned int k=0;
	  cout<<"Z="<<_subDetectorVec[i+1]->GetPlanes()->at(j)->GetHits()[k]->position()[2]<<" Y="<<_subDetectorVec[i+1]->GetPlanes()->at(j)->GetHits()[k]->position()[1]<<endl;
	  hits.push_back(_subDetectorVec[i+1]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
	      //}
	  
	  //if(checkNext) break;
	}

      //checkNext = true;
      //hits.push_back(_subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits()[0]);
      //cout<<"Z="<<_subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits()[0]->position()[2]<<" Y="<<_subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits()[0]->position()[1]<<endl;
	      
      //cout<<hits.size()<<endl;
 
      cout<<"Nhits="<<hits.size()<<endl;
      retVal = ChargeOne(hits,_subDetectorVec[i]);
      cout<<retVal <<endl;
      hits.clear();
      //checkNext = false;
      //break;

    }
  //retVal = ChargeOne(hits,_subDetectorVec[detectorNum]);
  //cout<<retVal <<endl;
  //hits.clear();


  return retVal;

}
*/

//*************************************************************
double Detector::CalculateChargeMomentum() {
//*************************************************************

// fill prev de_dx for all subdetectors.
// combine all of the guessess.



  double retVal = 0;

  sort( _subDetectorVec.begin(), _subDetectorVec.end(), sortDetectorsByZ() );

    for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      _subDetectorVec[i]->ClearEstimate();
      sort( _subDetectorVec[i]->GetPlanes()->begin(), _subDetectorVec[i]->GetPlanes()->end(), sortPlanesByZ() );
      if(i!=0)
	{
	  _subDetectorVec[i]->SetPrevde_dx(_subDetectorVec[i-1]->GetLength()*_subDetectorVec[i-1]->Getde_dx());
	}
    }

    for(unsigned int i=1;i<_subDetectorVec.size()-2;i++)
      {
	cout<<i<<endl;
	if(_subDetectorVec[i]->GetName() == "SFFFS0" ||
	   _subDetectorVec[i]->GetName() == "SFS" ||
	   _subDetectorVec[i]->GetName() == "SFFFFS1" ||
	   _subDetectorVec[i]->GetName() == "SFFS1" ||
	   _subDetectorVec[i]->GetName() == "SFFFS1"
	   )
	  {
	    /*
	    if(_subDetectorVec[i-1]->GetNHits()>0 && 
	       _subDetectorVec[i]->GetNHits()>0 && 
	       _subDetectorVec[i+1]->GetNHits()>0)

	    */
	    //if(_subDetectorVec[i-1]->GetNPlanes()==1 &&
	    // _subDetectorVec[i]->GetNPlanes()>1 &&
	    // _subDetectorVec[i+1]->GetNPlanes()>1)
	    //{
	    /*
		if(
	       _subDetectorVec[i]->GetPlanes()->at(0)->GetHits().size()>0 &&
	       _subDetectorVec[i-1]->GetPlanes()->at(0)->GetHits().size()>0 &&
	       _subDetectorVec[i]->GetPlanes()->at(1)->GetHits().size()>0 &&
	       _subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits().size()>0)
	    */

	    if( _subDetectorVec[i-1]->GetMaxPlane() &&_subDetectorVec[i]->GetMinPlane() && 
		_subDetectorVec[i]->GetMaxPlane() &&  _subDetectorVec[i+1]->GetMinPlane())
	      {
	    if(_subDetectorVec[i-1]->GetMaxPlane()->GetHits().size()>0 &&
	       _subDetectorVec[i]->GetMinPlane()->GetHits().size()>0 &&
	       _subDetectorVec[i]->GetMinPlane() !=  _subDetectorVec[i]->GetMaxPlane() &&
	       _subDetectorVec[i]->GetMaxPlane()->GetHits().size()>0 &&
	       _subDetectorVec[i+1]->GetMinPlane()->GetHits().size()>0)
	      {
		//cout<<"before theta1"<<endl;
		//double theta1 =  atan((_subDetectorVec[i]->GetPlanes()->at(0)->GetHits()[0]->position()[1] -_subDetectorVec[i-1]->GetPlanes()->at(0)->GetHits()[0]->position()[1])/
		//	      (_subDetectorVec[i]->GetPlanes()->at(0)->GetHits()[0]->position()[2] - _subDetectorVec[i-1]->GetPlanes()->at(0)->GetHits()[0]->position()[2]));
		
		
		double dy = (_subDetectorVec[i]->GetMinPlane()->GetHits()[0]->position()[1] -_subDetectorVec[i-1]->GetMaxPlane()->GetHits()[0]->position()[1]);
		double dz = (_subDetectorVec[i]->GetMinPlane()->GetHits()[0]->position()[2] - _subDetectorVec[i-1]->GetMaxPlane()->GetHits()[0]->position()[2]);

		



		//double theta1 =  atan((_subDetectorVec[i]->GetMinPlane()->GetHits()[0]->position()[1] -_subDetectorVec[i-1]->GetMaxPlane()->GetHits()[0]->position()[1])/
		//	      (_subDetectorVec[i]->GetMinPlane()->GetHits()[0]->position()[2] - _subDetectorVec[i-1]->GetMaxPlane()->GetHits()[0]->position()[2]));
		
		double theta1 = atan(dy/dz);

		double errZ = 10; //mm
		double errY = 150; //mm

		double errTheta1 = 1.0/(1.0+(dy/dz * dy/dz)) * 1.0/(1.0+(dy/dz * dy/dz)) * (2*errZ*errZ+2*errY*errY);

		errTheta1 = sqrt(errTheta1);


		//cout<<"before theta2"<<endl;
		//double theta2 =  atan((_subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits()[0]->position()[1] -_subDetectorVec[i]->GetPlanes()->at(1)->GetHits()[0]->position()[1])/
		//	      (_subDetectorVec[i+1]->GetPlanes()->at(0)->GetHits()[0]->position()[2] - _subDetectorVec[i]->GetPlanes()->at(1)->GetHits()[0]->position()[2]));
		

		dy = (_subDetectorVec[i+1]->GetMinPlane()->GetHits()[0]->position()[1] -_subDetectorVec[i]->GetMaxPlane()->GetHits()[0]->position()[1]);

		dz = (_subDetectorVec[i+1]->GetMinPlane()->GetHits()[0]->position()[2] - _subDetectorVec[i]->GetMaxPlane()->GetHits()[0]->position()[2]);

		double theta2 =atan(dy/dz);

		double errTheta2 = 1.0/(1.0+(dy/dz * dy/dz)) * 1.0/(1.0+(dy/dz * dy/dz)) * (2*errZ*errZ+2*errY*errY);

		errTheta2 = sqrt(errTheta2);

		//double theta2 =  atan((_subDetectorVec[i+1]->GetMinPlane()->GetHits()[0]->position()[1] -_subDetectorVec[i]->GetMaxPlane()->GetHits()[0]->position()[1])/
		//	      (_subDetectorVec[i+1]->GetMinPlane()->GetHits()[0]->position()[2] - _subDetectorVec[i]->GetMaxPlane()->GetHits()[0]->position()[2]));
		

		//cout<<"before delta1"<<endl;


		//double delta1 = theta2-theta1;

		double delta1 = theta1-theta2;

		double errDelta1 = errTheta2*errTheta2 + errTheta1*errTheta1;

		errDelta1 = sqrt(errDelta1);

		//double delta1 = theta1-theta2;


		//cout<<"before l"<<endl;
		double l = (_subDetectorVec[i]->GetPlanes()->at(1)->GetHits()[0]->position()[2]  - 
			    _subDetectorVec[i]->GetPlanes()->at(0)->GetHits()[0]->position()[2])*1000;
		//cout<<"before p"<<endl;
		double p = 0.3* 
		  _subDetectorVec[i]->GetBField()->vector(_subDetectorVec[i]->GetPlanes()->at(0)->GetHits()[0]->position())[0] *
		  l/ delta1;


		cout<<"Momentum="<<p<<endl;
		
		_subDetectorVec[i]->SetCharge(p/fabs(p));
		_subDetectorVec[i]->SetMomentum(fabs(p));
		_subDetectorVec[i]->SetError(p*errDelta1/delta1);

		double X0 = _subDetectorVec[i]->GetX0Eff();

		double chargeWidth = 13.6/(fabs(p)*(fabs(p)/sqrt(p*p+105.65*105.65))) *
					  sqrt(l/X0) *
					  (1+0.0036*log(l/X0));

		//double angularResolution = 1.0/


		double chargeProb =  1/sqrt(2* M_PI * chargeWidth) * exp(-(delta1-chargeWidth)*(delta1-chargeWidth)/(2*chargeWidth*chargeWidth));

		double negChargeProb =  1/sqrt(2* M_PI * chargeWidth) * exp(-(-delta1-chargeWidth)*(-delta1-chargeWidth)/(2*chargeWidth*chargeWidth));



		_subDetectorVec[i]->SetChargeProb(chargeProb);

		_subDetectorVec[i]->SetNegChargeProb(negChargeProb);


	      } // End if enough hits
	      } // End if good pointers
		//}


	  }
	else // Do not use leverarm. (In stack)
	  {
	    cout<<"In stack"<<endl;

	    std::vector<cluster*> hits;
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

	    cout<<_subDetectorVec[i]->GetName()<<" "<<hits.size()<<endl;

	    if(hits.size()>2)
	      {
		if(hits.size()>3){
		  vector<double> p;
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
		    vector<double> p;
		    cout<<"Using threePoint"<<endl;
		    p = ThreePointCircle(hits);
		     _subDetectorVec[i]->SetCharge(p[0]/fabs(p[0]));
		     _subDetectorVec[i]->SetMomentum(fabs(p[0]));
		     _subDetectorVec[i]->SetError(p[1]);
		     //cout<<"Returns="<<p<<endl;

		     // Set to crazy value, needs to be calculated after mean p;
		     _subDetectorVec[i]->SetChargeProb(1.0);
		     _subDetectorVec[i]->SetNegChargeProb(1.0);
		
		  }
	      }
	    //else
		 
	  }


      }//end for.


    double tempP =0;
    double tempErr =0;

for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      cout<<"subDetector="<<_subDetectorVec[i]->GetMomentum()<<endl;


      if(_subDetectorVec[i]->GetMomentum() == 0 || isinf(_subDetectorVec[i]->GetMomentum()) ||
	 isnan(_subDetectorVec[i]->GetMomentum()) || isnan(_subDetectorVec[i]->GetError()) ||
	 _subDetectorVec[i]->GetError() == 0 || isinf(_subDetectorVec[i]->GetError()))
	{
	  continue;
	}
      else
	{
	  cout<<"Using detector="<<_subDetectorVec[i]->GetName()<<endl;
	  tempP += _subDetectorVec[i]->GetMomentum()
	    /(_subDetectorVec[i]->GetError() *_subDetectorVec[i]->GetError() )
	    + _subDetectorVec[i]->GetPrevde_dx();

	  tempErr += 1/(_subDetectorVec[i]->GetError()*_subDetectorVec[i]->GetError());
	}

    }

 cout<<"TempP="<<tempP<<endl;
 cout<<"TempErr="<<tempErr<<endl;

 tempP=tempP/tempErr;

 cout<<"Final="<<tempP<<endl;


 double chargeProb = 1;
 double negChargeProb = 1;
for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      cout<<"chargeProb="<<chargeProb<<endl;
      cout<<"NchargeProb="<<negChargeProb<<endl;

      if(_subDetectorVec[i]->GetMomentum() == 0 || isinf(_subDetectorVec[i]->GetMomentum()) ||
	 isnan(_subDetectorVec[i]->GetMomentum()) || isnan(_subDetectorVec[i]->GetError()) ||
	 _subDetectorVec[i]->GetError() == 0 || isinf(_subDetectorVec[i]->GetError()))
	{
	  continue;
	}
      else
	{
	  if(_subDetectorVec[i]->GetChargeProb() == 1 &&  _subDetectorVec[i]->GetNegChargeProb() == 1)
	    {
	      // Need to calculate prob from meanP;

	      double subP = _subDetectorVec[i]->GetMomentum();
	      double subE = _subDetectorVec[i]->GetError();

	      double tempChargeProb=  1/sqrt(2* M_PI * fabs(subP)) * exp(-(subP-tempP)*(subP-tempP)/(2*subE*subE));
	      double tempNegChargeProb= 1/sqrt(2* M_PI * fabs(subP)) * exp(-(-subP-tempP)*(-subP-tempP)/(2*subE*subE));

	      if(!tempChargeProb || !tempNegChargeProb) 
		{
		  cout<<"Prob is 0 ?! Mom="<<subP<<"+-"<<subE<<endl;
		  continue;
		}
	      else
		{
		  chargeProb *= tempChargeProb;
		  negChargeProb *= tempNegChargeProb;
		}

	    }
	  else
	    {
	      chargeProb*=_subDetectorVec[i]->GetChargeProb();
	      negChargeProb*= _subDetectorVec[i]->GetNegChargeProb();
	    }
	
      cout<<_subDetectorVec[i]->GetName()<<" "
	  <<_subDetectorVec[i]->GetChargeProb()<<" "
	  <<_subDetectorVec[i]->GetNegChargeProb()<<endl;
	}

    }

 cout<<"chargeProb="<<chargeProb<<endl;
 cout<<"NchargeProb="<<negChargeProb<<endl;
 
 if(chargeProb>negChargeProb)
   {
     retVal = tempP;
   }
 else
   {
     retVal = -1*tempP;
   }



 if(retVal==0 || retVal >8000 || isnan(retVal)) retVal = 4000;

  return retVal;
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


  sort( _subDetectorVec.begin(), _subDetectorVec.end(), sortDetectorsByZ() );

  

  double retVal = 0;
  bool checkNext = false;
  unsigned int detectorNum = 0;

  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      cout<<"Detector="<<_subDetectorVec[i]->GetName()<<endl;
    }
  

  //for(unsigned int i=0;i<_subDetectorVec.size();i++)
  for(unsigned int i=0;i<_subDetectorVec.size();i++)
    {
      sort( _subDetectorVec[i]->GetPlanes()->begin(), _subDetectorVec[i]->GetPlanes()->end(), sortPlanesByZ() );

      if(_subDetectorVec[i]->GetName() == "S0") continue;
      if(_subDetectorVec[i]->GetName() == "SFFFS0") continue;
      if(_subDetectorVec[i]->GetName() == "SFS") continue;
      if(_subDetectorVec[i]->GetName() == "SFFFFS1") continue;
      if(_subDetectorVec[i]->GetName() == "S1") continue;

      cout<<"Detector="<<_subDetectorVec[i]->GetName()<<endl;


      for(unsigned int j=0;j<_subDetectorVec[i]->GetPlanes()->size();j++)
	{
	  for(unsigned int k=0;k<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits().size();k++)
	    {
	      cout<<"Z="<<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]->position()[2]<<" Y="<<_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]->position()[1]<<endl;
	      hits.push_back(_subDetectorVec[i]->GetPlanes()->at(j)->GetHits()[k]); //Fill planes with mult occupancy...
	    }
	}
    }
  cout<<"Nhits="<<hits.size()<<endl;
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
  cout<<retVal <<endl;
  hits.clear();
  
  //retVal = ChargeOne(hits,_subDetectorVec[detectorNum]);
  //cout<<retVal <<endl;
  //hits.clear();


  return retVal;
}

//*************************************************************
double Detector::ChargeOne(std::vector<cluster*>& hits,SubDetector* detector) {
//*************************************************************

  cout<<"In ChargeOne"<<endl;

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

    cout<<"FIELD="<<B[0]<<" "<<B[1]<<" "<<B[2]<<endl;
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

  return meansign;

}


//*************************************************************
//double Detector::Quadratic(std::vector<cluster*>& hits) {
vector<double> Detector::Quadratic(std::vector<cluster*>& hits) {
//*************************************************************

  int nMeas = hits.size();
  int minindex = nMeas;  

  vector<double> retVec;

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
    x[iMeas] = hits[iMeas]->vector()[0]; //mm to meter
    y[iMeas] = hits[iMeas]->vector()[1];
    z[iMeas] = hits[iMeas]->vector()[2];

    // get the b-field from the previous step
    //EVector B = detector->GetBField()->vector(pos);
    EVector B = GetSubDetector(hits[iMeas]->vector()[2])->GetBField()->vector(pos);

    cout<<"FIELD="<<B[0]<<" "<<B[1]<<" "<<B[2]<<endl;

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

  (void) gr->Fit("parfit2", "QNEM");

  cout<<"Param0="<<fun->GetParameter(0)<<endl;
  cout<<"Param1="<<fun->GetParameter(1)<<endl;
  cout<<"Param2="<<fun->GetParameter(2)<<endl;

  cout<<"Error0="<<fun->GetParError(0)<<endl;
  cout<<"Error1="<<fun->GetParError(1)<<endl;
  cout<<"Error2="<<fun->GetParError(2)<<endl;

  cout<<"x2="<<fun->GetChisquare()<<endl;

  cout<<"Px2="<<fun->GetProb()<<endl;

  // Positive particle -> Positive parameters -> Qtilde +
  // Negative Particle -> Negative parameters -> Qtilde -

  double qtilde = 1000* 0.3*prevB[0]*pow(1+pow(fun->GetParameter(1),2),3./2.)/
    (2*fun->GetParameter(2));

  //double eqtilde = 1000* 0.3*prevB[0]*pow(1+pow(fabs(fun->GetParameter(1))+fun->GetParError(1),2),3./2.)/
  //(2*(fabs(fun->GetParameter(2))-fun->GetParError(2)));


  //cout<<"errorCalc="<<fabs(fabs(qtilde)-fabs(eqtilde))<<endl;


  //double qerr = sqrt(pow(3.0/2.0 * 2 * fun->GetParError(1)/ fun->GetParameter(1),2) + pow(fun->GetParError(2)/ fun->GetParameter(2),2));


  // sigma^2 = (dp/da)^2 sigma_a^2 + (dp/db)^2 sigma_b^2

  double qerr = (qtilde / fun->GetParameter(2))* (qtilde / fun->GetParameter(2)) * fun->GetParError(2) * fun->GetParError(2) +
    qtilde * 3* fun->GetParameter(1)/ (1+fun->GetParameter(1)*fun->GetParameter(1)) *  fun->GetParError(1) *
    qtilde * 3* fun->GetParameter(1)/ (1+fun->GetParameter(1)*fun->GetParameter(1)) * fun->GetParError(1);
    
  
  qerr =sqrt(qerr);
  
  cout<<"errorCalc="<<qerr<<endl;
  
  //int meansign = (int)(qtilde/fabs(qtilde));

  //if(meansign == 0) meansign = 1;

  cout<<"quadratic="<<qtilde<<endl;

  //return qtilde;

  retVec.push_back(qtilde);
  retVec.push_back(qerr);

  return retVec;

}

//*************************************************************
vector<double> Detector::ThreePointCircle(std::vector<cluster*>& hits) {
//*************************************************************

//Using for instance the theory from:
//http://paulbourke.net/geometry/circlesphere/

  vector<double> retVec;

  double ma = (hits[1]->vector()[1]-hits[0]->vector()[1])/
    (hits[1]->vector()[2]-hits[0]->vector()[2]);

  double mb = (hits[2]->vector()[1]-hits[1]->vector()[1])/
    (hits[2]->vector()[2]-hits[1]->vector()[2]);

  double zc = (ma*mb*(hits[0]->vector()[1]-hits[2]->vector()[1])
	       +mb*(hits[0]->vector()[2]+hits[1]->vector()[2])
	       -ma*(hits[1]->vector()[2]+hits[2]->vector()[2]))/
    (2*(mb-ma));

  double yc = -1/mb *(zc-(hits[1]->vector()[2]+hits[2]->vector()[2])/2) +
    (hits[1]->vector()[2]+hits[2]->vector()[1])/2;

  double z = hits[1]->vector()[2]-zc;

  double y = hits[1]->vector()[1]-yc;

  double r = sqrt(z*z+y*y);

  EVector B = GetSubDetector(hits[1]->vector()[2])->GetBField()->vector(hits[1]->vector());

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

  double err = 1;
  
  cout<<"threepoint="<<p<<endl;

  retVec.push_back(p);
  retVec.push_back(err);

  //return p;
  return retVec;

}
