#ifndef _subDetector___
#define _subDetector___

#include <recpack/RecPackManager.h>
#include <mind/MINDfieldMapReader.h>
#include <mind/DeDxMap.h>

#include <mind/plane_info.h>

#include <vector>

using namespace Recpack;

class SubDetector{

public:
    
  SubDetector();

  SubDetector(dict::Key vol_name, EVector pos, std::vector<double> size);

  virtual ~SubDetector();

  dict::Key GetName(){return _name;};

  void SetProperties(MINDfieldMapReader* field,DeDxMap* map){_BFieldMap = field; _DeDxMap = map;};

  void SetSecondaryProperties(double wSc,
			      double wFe,
			      double X0Eff,
			      double de_dx,
			      double fieldScale,
			      double length){
    _wSc=wSc; _wFe=wFe; _X0Eff=X0Eff; _de_dx=de_dx; _fieldScale=fieldScale; _length=length;};


  void ClearEstimate()
  {
    _momentum = 0;
    _charge = 0;
    _chargeProb=0;
    _error=0;
  };


  double GetMomentum(){return _momentum;};
  void SetMomentum(double momentum){_momentum=momentum;};

  double GetCharge(){return _charge;};
  void SetCharge(double charge){_charge=charge;};


  double GetChargeProb(){return _chargeProb;};
  void SetChargeProb(double chargeProb){_chargeProb=chargeProb;};
  double GetNegChargeProb(){return _negChargeProb;};
  void SetNegChargeProb(double negChargeProb){_negChargeProb=negChargeProb;};

  double GetError(){return _error;};
  void SetError(double error){_error=error;};


 
double GetwSc(){return _wSc;};
double GetwFe(){return _wFe;};
double GetX0Eff(){return _X0Eff;};
double Getde_dx(){return _de_dx;};
double GetFieldScale(){return _fieldScale;};
double GetLength(){return _length;};


 void SetPrevde_dx(double de_dx){_previousde_dx = de_dx;};
 double GetPrevde_dx(){return _previousde_dx;};

  std::vector<plane_info*>* GetPlanes() {return &_planes;};

  //plane_info* GetMinPlane() {return _planes[0];};

  //plane_info* GetMaxPlane() {return _planes[ _planes.size()-1];};

  plane_info* GetMinPlane();
  plane_info* GetMaxPlane();

  int GetNPlanes() {return _planes.size();};

  double GetZMax(){return _pos[2]+_size[2]/2.0;};
  double GetZ(){return _pos[2];};
  double GetZMin(){return _pos[2]-_size[2]/2.0;};

  int GetNHits();

  MINDfieldMapReader* GetBField(){return _BFieldMap;};

 private:
  //double _posX;
  //double _posY;
  //double _posZ;
  //double _sizeX;
  //double _sizeY;
  //double _sizeZ; //length

  double _wSc;
  double _wFe;
  double _X0Eff;
  double _de_dx;
  double _previousde_dx;
  double _fieldScale;
  double _length;

  double _momentum;
  double _charge;
  double _chargeProb;
  double _negChargeProb;
  double _error;
  

  EVector _pos;
  std::vector<double> _size;

  dict::Key _name;

  MINDfieldMapReader* _BFieldMap;

  DeDxMap* _DeDxMap;
  std::vector<plane_info*> _planes;


};



#endif 
