#include <mind/detector.h>

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
