#ifndef _subDetector___
#define _subDetector___

#include <recpack/RecPackManager.h>
#include <mind/MINDfieldMapReader.h>
#include <mind/DeDxMap.h>

//#include <mind/plane_info.h>

#include <vector>

using namespace Recpack;

class SubDetector{

public:
    
  SubDetector();

  SubDetector(dict::Key vol_name, EVector pos, std::vector<double> size);

  virtual ~SubDetector();

  dict::Key GetName(){return _name;};

  void SetProperties(MINDfieldMapReader* field,DeDxMap* map){_BFieldMap = field; _DeDxMap = map;};

  //std::vector<plane_info*>* GetPlanes() {return &_planes;};
  //int GetNPlanes() {return _planes.size();};

 private:
  double _posX;
  double _posY;
  double _posZ;
  double _sizeX;
  double _sizeY;
  double _sizeZ; //length

  EVector _pos;
  std::vector<double> _size;

  dict::Key _name;

  MINDfieldMapReader* _BFieldMap;

  DeDxMap* _DeDxMap;
  //std::vector<plane_info*> _planes;


};



#endif 
