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
  // line x1, x2. point x0
  // Distance abs (x2-x1) x (x1-x0) / abs (x2-x1)
  // rewrite as abs(a x b) / abs(a) 
  
  double ax = _kx;
  double ay = _ky;
  double az = _kz;

  double bx = _mx - hit->position()[0];
  double by = _my - hit->position()[1];
  double bz = _mz - hit->position()[2];
  // cross-product
  double cx = ay*bz - az*by;
  double cy = az*bx - ax*bz;
  double cz = ax*by - ay*bx;

  double d = sqrt(cx*cx+cy*cy+cz*cz);

  d/= sqrt(ax*ax+ay*ay+az*az);
  

  //double d = sqrt(CalculateRX(hit)*CalculateRX(hit)+
  //		  CalculateRY(hit)*CalculateRY(hit));

  return d;
}
//*************************************************************
double Line::CalculateRX(cluster* hit) {
//double Line::CalculateR(Node* hit) {
//*************************************************************
  double t = ( hit->position()[2]-_mz)/_kz; 
  double x = _kx*t+_mx;
  double d = fabs(x-hit->position()[0]);

  return d;
}
//*************************************************************
double Line::CalculateRY(cluster* hit) {
//double Line::CalculateR(Node* hit) {
//*************************************************************
  double t = ( hit->position()[2]-_mz)/_kz; 
  double y = _ky*t+_my;
  double d = fabs(y-hit->position()[1]);

  return d;
}
//*************************************************************
void Line::Equation() {
//*************************************************************

  std::cout<<" x = "<<_kx<<" * t + "<<_mx<<endl;
  std::cout<<" y = "<<_ky<<" * t + "<<_my<<endl;
  std::cout<<" z = "<<_kz<<" * t + "<<_mz<<endl;
}

//*************************************************************
void Line::PointAtZ(double z) {
//*************************************************************
  double t = (z-_mz)/_kz;
  
  std::cout<<" x = "<<_kx*t+_mx<<endl;
  std::cout<<" y = "<<_ky*t+_my<<endl;
  std::cout<<" z = "<<_kz*t+_mz<<endl;
}

