#include <mind/line.h>

//*************************************************************

Line::Line(cluster* downStream, cluster* upStream) { 
  

  _mx= downStream->position()[0];
  _my= downStream->position()[1];
  _mz= downStream->position()[2];

  _kx = upStream->position()[0]- _mx;
  _ky = upStream->position()[1]- _my;
  _kz = upStream->position()[2]- _mz;

  //hits.push_back(downStream);
  //hits.push_back(upStream);

}


Line::Line(Node* downStream, Node* upStream) { 
  

  _mx= downStream->measurement().position()[0];
  _my= downStream->measurement().position()[1];
  _mz= downStream->measurement().position()[2];

  _kx = upStream->measurement().position()[0]- _mx;
  _ky = upStream->measurement().position()[1]- _my;
  _kz = upStream->measurement().position()[2]- _mz;

  //hits.push_back(downStream);
  //hits.push_back(upStream);

}


//*************************************************************
Line::~Line() {
//*************************************************************


}

//*************************************************************
double Line::CalculateR(cluster* hit) {
//double Line::CalculateR(Node* hit) {
//*************************************************************

  //double dx = _mx - hit->position()[0] - (_mx - hit->position()[0])*_kx*_kx;
  //double dy = _my - hit->position()[1] - (_my - hit->position()[1])*_ky*_ky;
  //double dz = _mz - hit->position()[2] - (_mz - hit->position()[2])*_kz*_kz;

  //return sqrt(dx*dx+dy*dy+dz*dz);

  double param = (hit->position()[2] - _mz)/_kz;
  //double param = (hit->measurement().position()[2] - _mz)/_kz;

  std::cout<<"param="<<param<<endl;

  double dx = _mx + _kx * param - hit->position()[0];
  //double dx = _mx + _kx * param - hit->measurement().position()[0];

  double dy = _my + _ky * param - hit->position()[1];

  //double dy = _my + _ky * param - hit->measurement()position()[1];

  std::cout<<"dx="<<dx<<endl;
  std::cout<<"dy="<<dy<<endl;

  return sqrt(dx*dx+dy*dy);
  


}
//*************************************************************
void Line::Equation() {
//*************************************************************

  std::cout<<" x = "<<_kx<<" * t + "<<_mx<<endl;
  std::cout<<" y = "<<_ky<<" * t + "<<_my<<endl;
  std::cout<<" z = "<<_kz<<" * t + "<<_mz<<endl;
}

