#include <cluster.h>

cluster::cluster() : Measurement()
{
  _eng = 0;
  _muProp = 0;
  _nhit = 0;
  _time = 0;
}

cluster::~cluster()
{
}

void cluster::add_hit(bhep::hit* dep)
{

  _voxes.push_back( dep );

  _eng += dep->ddata( "TotalEng" );

   //testbeam _muProp += dep->ddata( "MuonProportion" );
  
  _nVox = _voxes.size();

  _nhit += dep->idata( "NoPoints" );

  _time = dep->ddata( "HitTime" );

  _mother_particle = dep->vsdata("mother_particle");

  _IsTASD = dep->idata( "IsTASD" );


  set_hv("energy", HyperVector(_eng,0));
  //testbeam set_hv("MuonProp", HyperVector(_muProp/_nhit,0));
  set_hv("hittime", HyperVector(_time,0));
}
