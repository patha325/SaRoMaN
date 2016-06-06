#include "MindCryInterface.h"
#include "MindConfigService.h"
// #include "MindDetectorGeometry.h"
#include "SciNearDetectorGeometry.h"
#include "MindDetectorConstruction.h"
#include "MindParamStore.h"
#include "MindLookupTable.h"
#include <G4ParticleTable.hh>
#include <G4Event.hh>
#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <Randomize.hh>
#include <G4Material.hh>
#include <G4PrimaryParticle.hh>

MindCryInterface::MindCryInterface()
{
  Initialize();
}

MindCryInterface::~MindCryInterface()
{
}

void MindCryInterface::Initialize()
{
  
  // Read the cry input
  // std::string = config.GetSParam("CRYString");

  char len[6];
  double detwidth = MindConfigService::Instance().Geometry().GetDParam("length") * mm / m;
  int n = sprintf ( len, "%d", int(detwidth));
  std::string setupString = "returnNeutrons 1 ";
  setupString += "returnProtons 1 ";
  setupString += "returnGammas 1 ";
  setupString += "returnPions 1 ";
  setupString += "returnKaons 1 ";
  setupString += "date 11-11-2012 ";
  setupString += "latitude 41.83 ";
  setupString += "altitude 0 ";
  setupString += "subboxLength ";
  setupString += len;

  // CRY setup and generator initialization
  // std::string setupString = MindConfigService::Instance().Generation().GetSParam("crystring");
  std::cout<<setupString<<std::endl;
  std::string datastring = MindConfigService::Instance().Generation().GetSParam("crydata");
  _surface = MindConfigService::Instance().Generation().PeekIParam("crysurface") ? 
    MindConfigService::Instance().Generation().GetIParam("crysurface"):1;
  CRYSetup *setup=new CRYSetup(setupString, datastring);
			       
  gen = new CRYGenerator(setup);
  
  vect=new std::vector<CRYParticle*>;

  particle_table = G4ParticleTable::GetParticleTable();

  _fspdg = bhep::vdouble(6);
  for(int i=0; i<6; i++) _fspdg[i] = 0;
    
}
void MindCryInterface::GeneratePrimaryVertex(G4Event* event)
{

  // particle is generated at time zero
  G4double time = 0.;

  bhep::event& bevt = bhep::bhep_svc::instance()->get_event();

  MindDetectorConstruction* detConstr = (MindDetectorConstruction*)
    G4RunManager::GetRunManager()->GetUserDetectorConstruction();
  bevt.add_property( "IntType", "CRY");
  G4String particleName;
  vect->clear();
  gen->genEvent(vect);
  bevt.add_property("Charm", 0);
  bevt.add_property("CharmHad", 0);
  bevt.add_property("Q2",0);
  bevt.add_property("EngTrans",0);
  bool fillbevt=false;
  double totEnergy = 0;
  
  double detwidth = MindConfigService::Instance().Geometry().GetDParam("width") * mm  / 2.;
  double detheight = MindConfigService::Instance().Geometry().GetDParam("height") * mm / 2.;
  double detlength = MindConfigService::Instance().Geometry().GetDParam("length") * mm / 2.;
  for ( unsigned j=0; j<vect->size(); j++) {
    particleName=CRYUtils::partName((*vect)[j]->id());
    
    // MIND coordinate system is rotated with respect to that used for CRY
    G4ThreeVector position; // initial position co-ordinate depends on the face of the detector in question
    if(_surface==1){ // top of detector
      position = G4ThreeVector( (*vect)[j]->y()*m, (*vect)[j]->z()*m + detlength, (*vect)[j]->x()*m);
      bevt.set_vertex( (*vect)[j]->y()*m, (*vect)[j]->z()*m + detlength, (*vect)[j]->x()*m );
    } else if(_surface==2){ // right face of detector
      position = G4ThreeVector( (*vect)[j]->z()*m + detlength, (*vect)[j]->x()*m, (*vect)[j]->y()*m);
      bevt.set_vertex( (*vect)[j]->z()*m + detlength, (*vect)[j]->x()*m, (*vect)[j]->y()*m );
    } else if(_surface==3){ // left face of detector
      position = G4ThreeVector( (*vect)[j]->z()*m - detlength, (*vect)[j]->x()*m, (*vect)[j]->y()*m);
      bevt.set_vertex( (*vect)[j]->z()*m - detlength, (*vect)[j]->x()*m, (*vect)[j]->y()*m );
    } else if(_surface==4){ // upstream face of detector
      position = G4ThreeVector( (*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->z()*m - detlength);
      bevt.set_vertex( (*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->z()*m - detlength );
    } else if(_surface==5){ // downstream face of detector
      position = G4ThreeVector( (*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->z()*m + detlength);
      bevt.set_vertex( (*vect)[j]->x()*m, (*vect)[j]->y()*m, (*vect)[j]->x()*m + detlength);
    } else if(_surface==6){ // bottom face of detector
      position = G4ThreeVector( (*vect)[j]->y()*m, (*vect)[j]->z()*m - detlength, (*vect)[j]->x()*m);
      bevt.set_vertex( (*vect)[j]->y()*m, (*vect)[j]->z()*m - detlength, (*vect)[j]->x()*m );
    } 
    G4PrimaryVertex* vertex = new G4PrimaryVertex(position, (*vect)[j]->t());
    
    
    G4ParticleDefinition* particle_definition = particle_table->FindParticle((*vect)[j]->PDGid());
    G4double mass   = particle_definition->GetPDGMass();
    G4double charge    = particle_definition->GetPDGCharge();
    G4double energy = (*vect)[j]->ke()*MeV + mass;
    G4double pmom = std::sqrt(energy*energy - mass*mass);
    G4double px = pmom * (*vect)[j]->v();
    G4double py = pmom * (*vect)[j]->w();
    G4double pz = pmom * (*vect)[j]->u();

    totEnergy += energy;
    std::cout<<"Particle "<<j<<" has PID "<<(*vect)[j]->PDGid()<<" at Energy "<<energy<<std::endl;
    G4PrimaryParticle* particle = new G4PrimaryParticle(particle_definition, px, py, pz);
    particle->SetMass(mass);
    particle->SetCharge(charge);
    particle->SetPolarization(0.,0.,0.);
    vertex->SetPrimary(particle);
    if (!fillbevt)
      if ( vect->size() == 1 || vect->size() > 1 && fabs((*vect)[j]->PDGid())==13){
	bevt.add_property("partEnergy", energy);
	bevt.add_property("intpart", (*vect)[j]->PDGid());
	fillbevt=true;
      }
    event->AddPrimaryVertex(vertex);
  }

  
  for( int j=0; j < (vect->size() <=6? vect->size():6); j++)
    _fspdg[j] = (*vect)[j]->PDGid();
  
  bevt.add_property("nuEnergy", totEnergy);
  bevt.add_property("nuType", (*vect)[0]->id());

  bhep::vdouble had4P = bhep::vdouble(4);
  
  for (bhep::vdouble::iterator it1=had4P.begin();it1 != had4P.end();it1++)
    (*it1) = 0.;

  bevt.add_property("had4vec", had4P);
  bevt.add_property("npip",0);
  bevt.add_property("npi0",0);
  bevt.add_property("npim",0);
  bevt.add_property("nkp",0);
  bevt.add_property("nk0",0);
  bevt.add_property("nkm",0);
  bevt.add_property("hadEInNucleus",0.0);
  bevt.add_property("fspdg", _fspdg);
}

