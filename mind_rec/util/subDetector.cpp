#include <mind/subDetector.h>

//*************************************************************
SubDetector::SubDetector() {
//*************************************************************
  
}

//*************************************************************
SubDetector::SubDetector(dict::Key vol_name, EVector pos, std::vector<double> size) {
//*************************************************************
  _name = vol_name;
  _pos = pos;
  _size = size;
}


//*************************************************************
SubDetector::~SubDetector() {
//*************************************************************


}

//*************************************************************
int SubDetector::GetNHits() {
//*************************************************************

  int ret_val =0;

  for(unsigned int i=0;i<_planes.size();i++)
    {
      ret_val+=_planes[i]->GetNHits();
    }
  return ret_val;
}

//*************************************************************
 plane_info*  SubDetector::GetMaxPlane() {
//*************************************************************

   plane_info* ret_val;

   if(_planes.size())
     {
       ret_val = _planes[_planes.size()-1];
     }
   else
     {
       ret_val = NULL;
     }

  return ret_val;
}

//*************************************************************
 plane_info*  SubDetector::GetMinPlane() {
//*************************************************************

   plane_info* ret_val;

   if(_planes.size())
     {
       ret_val = _planes[0];
     }
   else
     {
       ret_val = NULL;
     }

  return ret_val;
}
