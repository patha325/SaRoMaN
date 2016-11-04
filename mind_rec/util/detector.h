#ifndef _detector___
#define _detector___

#include <recpack/RecPackManager.h>
//#include <mind/SetupSk.h>
//#include <mind/MINDfieldMapReader.h>
//#include <mind/DeDxMap.h>

#include <mind/subDetector.h>

#include <vector>


using namespace Recpack;

class Detector{

public:
    
  Detector();
  virtual ~Detector();

  int GetNSubDetectors(){return _subDetectorVec.size();};
  std::vector<SubDetector*>* GetSubDetectorVec(){return &_subDetectorVec;};

  SubDetector* GetSubDetector(dict::Key vol_name);

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
