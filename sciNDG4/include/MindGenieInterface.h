// ----------------------------------------------------------------------------
///  \file   MindGenieInterface.h
///  \brief
///  
///  \author  A Laing <a.laing@physics.gla.ac.uk>
///  \date    27 September 2010
///  \version $Id: MindGenieInterface.h 309 2009-06-22 08:30:02Z alaing $
///
///  Copyright (c) 2009 - IFIC Neutrino Group
// ----------------------------------------------------------------------------

#ifndef __GENIE_INTERFACE__
#define __GENIE_INTERFACE__

#include <G4VPrimaryGenerator.hh>
#include <globals.hh>
//#include <fstream>

#include <TSystem.h>
#include <TFile.h>
#include <TTree.h>
#include <TIterator.h>

#include "EVGCore/EventRecord.h"
#include "GHEP/GHepParticle.h"
#include "Ntuple/NtpMCFormat.h"
#include "Ntuple/NtpMCTreeHeader.h"
#include "Ntuple/NtpMCEventRecord.h"
#include "Messenger/Messenger.h"
#include "PDG/PDGCodes.h"
#include "Utils/CmdLnArgParser.h"

#include "MindConfigService.h"
#include "MindParamStore.h"
#include "MindDetectorConstruction.h"
#include "SciNearDetectorGeometry.h"
#include "MindLookupTable.h"

#include <G4RunManager.hh>
#include <G4VUserDetectorConstruction.hh>
#include <Randomize.hh>
#include <G4Material.hh>
#include <G4PrimaryParticle.hh>
// #include "Utils/CmdLineArgParserUtils.h"
// #include "Utils/CmdLineArgParserException.h"

#include <bhep/bhep_svc.h>

class G4ParticleDefinition;
class G4PrimaryParticle;
class G4Event;

using namespace genie;
using namespace CLHEP;

class MindGenieInterface: public G4VPrimaryGenerator
{
 public:
  /// Constructor
  MindGenieInterface();
  /// Destructor
  /* ~MindGenieInterface() {} */
  ~MindGenieInterface();
  
  void GeneratePrimaryVertex(G4Event*);

 private:
  void Initialize();

  void Finalize();

  G4PrimaryParticle*
    CreatePrimaryParticle(GHepParticle& part, G4int PDG, bool primLep=false);
  /* G4PrimaryParticle* */
/*     CreatePrimaryParticle(GHepParticle& part, G4int PDG); */

  G4String SelectVertexRegion();
  G4double GetTargetProb();
  /* G4String InteractionType(G4int code); */

private:
  TFile* _genieFiles[2]; ///< Genie input data files
  TTree* _genieData[2]; ///< Data trees
  /* NtpMCTreeHeader* _genieHeader[2]; ///<genie headers. */
  NtpMCEventRecord* _mcrec[2];
  
  G4ParticleDefinition* _particle_definition;
  G4double _mass; ///< Particle mass
  G4double _charge; ///< Particle charge
  G4int _vert_mat;
  G4double _FeTargetProp;
  G4String _vtx_location;
  G4int _fsl_select;

  G4int _evCount[2]; //< Event counter.

  G4int _tasd; //< a tasd detector?

  //Vertex in the case a fixed point is requested.
  G4int _fvecReg;
  G4ThreeVector _fvec;

  bhep::vdouble _had4P;

  bhep::vdouble _fspdg;
  G4double _hadEInNucleus; // total energy of hadron inside nucleus
};

#endif
