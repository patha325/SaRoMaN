#include <plane_info.h>
#include <recpack/stc_tools.h>
#include <TMath.h>



using namespace bhep;


//**********************************************************************
plane_info::plane_info() {
  //**********************************************************************
}


//**********************************************************************
plane_info::plane_info(int no, double zpos, const bhep::gstore& pstore) {
  //**********************************************************************

  _planeNo = no;
  _zPos = zpos;
  _engPlane=0;
  _store = pstore;
  
  /*detX = _store.fetch_dstore("MIND_x") * bhep::m;
   detY = _store.fetch_dstore("MIND_y") * bhep::m;

   if(_store.find_dstore("WLSatten"))
     WLSAtten = _store.fetch_dstore("WLSatten");
   else						
   WLSAtten = 5000.;*/
}

//**********************************************************************
plane_info::plane_info(int no, double zpos, double rpos, const bhep::gstore& pstore) {
  //**********************************************************************

  _planeNo = no;
  _zPos = zpos;
  _engPlane=0;
  _rPos= rpos;    
  _store = pstore;
}


//***********************************************************************
plane_info::~plane_info() {
  //***********************************************************************
  //  for(int i=0; i<(int)_hits.size(); i++)
  // delete _hits[i];
     
    _hits.clear();
  

}


//**********************************************************************
void plane_info::AddHit(cluster* hit) {
  //**********************************************************************


  _hits.push_back(hit);

  _engPlane += correctEdep(hit->get_eng() * MeV, hit);
  

}


//***********************************************************************
void plane_info::reset(){
  //***********************************************************************
  
  _zPos = 0;
  _planeNo = 0;
}





//***********************************************************************
double plane_info::correctEdep(double edep, cluster* hit)
  //***********************************************************************
{

  double sum1 = 0;
  

  double X = hit->position()[0];
  double Y = hit->position()[1];
  double Z = hit->position()[2];

  double detX, detY, detZ, vdetX, vdetY, vdetZ, WLSAtten, octGeom;


  /// retrieve the values from parameter file
  // _store = pstore;
  
  detX = _store.fetch_dstore("MIND_x") * m;
  detY = _store.fetch_dstore("MIND_y") * m;
  detZ = _store.fetch_dstore("MIND_z") * m;
  vdetX = _store.find_dstore("vertex_x")?
    _store.fetch_dstore("vertex_x") * m : 0.0;
  vdetY = _store.find_dstore("vertex_y")?
    _store.fetch_dstore("vertex_y") * m : 0.0;;
  vdetZ = _store.find_dstore("vertexDepth")?
    _store.fetch_dstore("vertexDepth") * m : 0.0;

  if(_store.find_dstore("WLSatten"))
    WLSAtten = _store.fetch_dstore("WLSatten");
  else						
    WLSAtten = 5000.;

    
  octGeom  = _store.find_istore("IsOctagonal") ?
    _store.fetch_istore("IsOctagonal"): 1;

  ///for the time being provided manually
  /* detX = 14* m;
  detY = 14* m;
  WLSAtten = 5* m;
  int octGeom =1;
  */
  
  double slope = octGeom==1 ? (detY - detX*tan(atan(1)/2.))/
    (detY*tan(atan(1)/2.) - detX) : -1.0;
  if(octGeom==0 && Z > (vdetZ - detZ)/2.){ // Assume a cylindrical geometry
    double xedge = detX * sqrt(1 - 4.* pow(Y/detY, 2))/2.;
    double yedge = detY * sqrt(1 - 4.* pow(X/detX, 2))/2.;
    sum1 =  exp( -(xedge - fabs(X))/WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/WLSAtten );
  }
  else if(octGeom==0 && Z < (vdetZ - detZ)/2.){ // Assume a cylindrical geometry
    double xedge = detX/2.;
    double yedge = detY/2.;
    sum1 =  exp( -(xedge - fabs(X))/WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/WLSAtten );
  }
  else if(octGeom==1) {
    //need to take into account drift distance to closest and furthest edge.
    if((fabs(X) < detY*tan(atan(1)/2.)/2 &&
	fabs(Y) < detX*tan(atan(1)/2.)/2 ) ){
      sum1 =  exp( -(detX/2 - fabs(X))/WLSAtten );
      sum1 += exp( -(detX/2 + fabs(X))/WLSAtten );
      sum1 += exp( -(detY/2 - fabs(Y))/WLSAtten );
      sum1 += exp( -(detY/2 + fabs(Y))/WLSAtten );
    }
    else if(fabs(X) > detY*tan(atan(1)/2.)/2 &&
	    fabs(Y) < detX*tan(atan(1)/2.)/2  ){
      double xedge = detX/2 + (fabs(Y) - 
				detX/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge  - fabs(X))/WLSAtten );
      sum1 += exp( -(xedge  + fabs(X))/WLSAtten );
      sum1 += exp( -(detY/2 - fabs(Y))/WLSAtten );
      sum1 += exp( -(detY/2 + fabs(Y))/WLSAtten );
    }
    else if(fabs(X) < detY*tan(atan(1)/2.)/2 &&
	    fabs(Y) > detX*tan(atan(1)/2.)/2  ){
      double yedge = detY/2 + (fabs(X) - 
				detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(detX/2 - fabs(X))/WLSAtten );
      sum1 += exp( -(detX/2 + fabs(X))/WLSAtten );
      sum1 += exp( -(yedge  - fabs(Y))/WLSAtten );
      sum1 += exp( -(yedge  + fabs(Y))/WLSAtten );
    }
    else if(fabs(X) > detY*tan(atan(1)/2.)/2 &&
	    fabs(Y) > detX*tan(atan(1)/2.)/2  ){
      double xedge = detX/2 + (fabs(Y) - 
				detX/2.*tan(atan(1)/2.))*slope;
      double yedge = detY/2 + (fabs(X) - 
				detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge - fabs(X))/WLSAtten );
      sum1 += exp( -(xedge + fabs(X))/WLSAtten );
      sum1 += exp( -(yedge - fabs(Y))/WLSAtten );
      sum1 += exp( -(yedge + fabs(Y))/WLSAtten );
    }
  }
  else if(octGeom==2) {
    //need to take into account drift distance to closest and furthest edge.
    if((fabs(X) < detY*tan(atan(1)/2.)/2. + detX/4.&&
	fabs(Y) < detY*tan(atan(1)/2.)/2. ) ){
      sum1 =  exp( -(detX/2 - fabs(X))/WLSAtten );
      sum1 += exp( -(detX/2 + fabs(X))/WLSAtten );
      sum1 += exp( -(detY/2 - fabs(Y))/WLSAtten );
      sum1 += exp( -(detY/2 + fabs(Y))/WLSAtten );
    }
    else if(fabs(X) > detY*tan(atan(1)/2.)/2. + detX/4. &&
	    fabs(Y) < detY*tan(atan(1)/2.)/2.  ){
      double xedge = detX/2 + (fabs(Y) - 
				detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge  - fabs(X))/WLSAtten );
      sum1 += exp( -(xedge  + fabs(X))/WLSAtten );
      sum1 += exp( -(detY/2 - fabs(Y))/WLSAtten );
      sum1 += exp( -(detY/2 + fabs(Y))/WLSAtten );
    }
    else if(fabs(X) < detY*tan(atan(1)/2.)/2 + detX/4. &&
	    fabs(Y) > detY*tan(atan(1)/2.)/2  ){
      double yedge = detX/2 + (fabs(X) - 
				detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(detX/2 - fabs(X))/WLSAtten );
      sum1 += exp( -(detX/2 + fabs(X))/WLSAtten );
      sum1 += exp( -(yedge  - fabs(Y))/WLSAtten );
      sum1 += exp( -(yedge  + fabs(Y))/WLSAtten );
    }
    else if(fabs(X) > detY*tan(atan(1)/2.)/2 + detX/4. &&
	    fabs(Y) > detY*tan(atan(1)/2.)/2  ){
      double xedge = detX/2 + (fabs(Y) - 
				detX/2.*tan(atan(1)/2.))*slope;
      double yedge = detX/2 + (fabs(X) - 
				detY/2.*tan(atan(1)/2.))*slope;
      sum1 =  exp( -(xedge - fabs(X))/WLSAtten );
      sum1 += exp( -(xedge + fabs(X))/WLSAtten );
      sum1 += exp( -(yedge - fabs(Y))/WLSAtten );
      sum1 += exp( -(yedge + fabs(Y))/WLSAtten );
    }
  }
  else if(octGeom == -2 || octGeom == 3){ // Assume a square cross-section
    double xedge = detY/2.;
    double yedge = detY/2.;
    sum1 =  exp( -(xedge - fabs(X))/WLSAtten );
    sum1 += exp( -(xedge + fabs(X))/WLSAtten );
    sum1 += exp( -(yedge - fabs(Y))/WLSAtten );
    sum1 += exp( -(yedge + fabs(Y))/WLSAtten );
  }

  double corrEng = 4*edep/sum1;
  return corrEng;
}

//*********************************************************
int plane_info::mostSeparatedHitIndex(){
  //*******************************
  if (_hits.size() <= 2) return 0;
  else if (_hits.size() >  2){
    double maxSumRDiff = -1.;
    int iMaxSumRDiff = -1;
    for ( int i=0; i< _hits.size(); i++){
      // select a single hit and difine a sum of differncese in position in the plane
      double sumRDiff=0.;
      double x0 = _hits[i]->position()[0];
      double y0 = _hits[i]->position()[1];
      
      for ( int j=0; j< _hits.size(); j++){
	if ( i == j ) continue;
	double x1 = _hits[j]->position()[0];
	double y1 = _hits[j]->position()[1];
	sumRDiff += sqrt((x0 - x1)*(x0 - x1) + (y0 - y1)*(y0 - y1));
      }
      
      if ( sumRDiff > maxSumRDiff ) { // update maximum
	maxSumRDiff = sumRDiff;
	iMaxSumRDiff = i;
      }
    }
    return iMaxSumRDiff;
  }
}
