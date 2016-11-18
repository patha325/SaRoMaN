#ifndef _detector___
#define _detector___

#include <recpack/RecPackManager.h>
//#include <mind/SetupSk.h>
//#include <mind/MINDfieldMapReader.h>
//#include <mind/DeDxMap.h>

#include <mind/subDetector.h>

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>

#include <math.h>

#include <vector>


using namespace Recpack;

class Detector{

public:
    
  Detector();
  virtual ~Detector();

  int GetNSubDetectors(){return _subDetectorVec.size();};
  std::vector<SubDetector*>* GetSubDetectorVec(){return &_subDetectorVec;};

  SubDetector* GetSubDetector(dict::Key vol_name);
  SubDetector* GetSubDetector(double zPosition);

  double Helix(std::vector<cluster*>& hits);
  
  double LeverArmQuadratic(vector<double>& debug);
  //double LeverArm(unsigned int i);

  vector<double> LeverArm(unsigned int i);

  double CalculateChargeMomentum(vector<double>& debug);

  int GetNPlanes();

  int GetNHits();

  void Reset();

  double CalculateCharge();

  double ChargeOne(std::vector<cluster*>& hits, SubDetector* detector);

  //double Quadratic(std::vector<cluster*>& hits);
  vector<double> Quadratic(std::vector<cluster*>& hits);

  vector<double> ThreePointCircle(std::vector<cluster*>& hits);

 private:

  double _posX;
  double _posY;
  double _posZ;
  double _sizeX;
  double _sizeY;
  double _sizeZ; //length
  std::vector<SubDetector*> _subDetectorVec;


};



#endif 
