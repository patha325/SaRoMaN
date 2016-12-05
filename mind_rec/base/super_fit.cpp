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

//void fitter::ComputeMomFromParabola(const Trajectory& traj, int nplanes, int firsthit, EVector& V){
//MomentumFromCurvature(const Trajectory& startTrack, int startPoint,double minP){
//*****************************************************************************
double super_fit::MomentumFromCurvature(const Trajectory& traj, int startPoint,double minP, vector<double>& debug){
  //*****************************************************************************
 
  double meansign = _supergeom.getDetectorModel()->CalculateChargeMomentum(debug);

  //return fabs(meansign);

  return meansign;

  //std::cout<<"Momentum guess from polynomial fit: p/q = "<<1./V[5]<<std::endl;
}

//***********************************************************************
double super_fit::CalculateCharge(const Trajectory& startTrack) {
  //***********************************************************************

  vector<double> debug;

  //double meansign = _supergeom.getDetectorModel()->CalculateCharge();

  double meansign = _supergeom.getDetectorModel()->CalculateChargeMomentum(debug);

  cout<<"returned in super_fit::meansign="<<meansign<<endl;

  meansign/=fabs(meansign);

  //double meansign =Helper(startTrack,0,0);
  
  return meansign;
}

