// ----------------------------------------------------------------------------
//  $Id: MindDetectorConstruction.cc 543 2014-11-01 23:03:10Z  $
//
//  Author : J Martin-Albo <jmalbos@ific.uv.es>
//  Created: 14 Jun 2008
//
//  Copyright (c) 2008, 2009 -- IFIC Neutrino Group
// ----------------------------------------------------------------------------

#include "MindDetectorConstruction.h"
#include "SciNearDetectorGeometry.h"
#include "MindMaterialsList.h"
#include "MindConfigService.h"
#include "MindParamStore.h"
#include "MindBarSD.h"

#include <G4Box.hh>
#include <G4LogicalVolume.hh>
#include <G4PVPlacement.hh>
#include <G4Material.hh>
#include <G4VisAttributes.hh>
#include <Randomize.hh>

//ROOT
/*
#include "TROOT.h"
#include "TFile.h"
#include "TSystem.h"
#include "Cintex/Cintex.h"

#include "MindGeoTree.h"
*/
#include "bhep/bhep_svc.h"

MindDetectorConstruction::MindDetectorConstruction():
  _detector(0)
{
  const MindParamStore& config = 
    MindConfigService::Instance().Geometry();
  _use_gdml = config.PeekIParam("useGDML") ? 
    ( config.GetIParam("useGDML")==0 ? false : true ) : false;
  _write_gdml = config.PeekIParam("writeGDML") ?
    ( config.GetIParam("writeGDML")==0 ? false : true ) : false;
  _gdml_file_name = config.PeekSParam("GDMLFileName") ? 
    config.GetSParam("GDMLFileName") : "";

  _myfile.open (config.GetSParam("xml_parsed"));
}



MindDetectorConstruction::~MindDetectorConstruction()
{
  _myfile.close();
  delete _detector;
}



G4VPhysicalVolume* MindDetectorConstruction::Construct()
{
  G4PVPlacement* world_physi;
  if( _write_gdml || !_use_gdml ){

    MindMaterialsList::DefineMaterials();
    
    
    // WORLD volumes definition .......................................
    // The WORLD is an empty (vacuum-filled) cube of 100m side.
    
    G4double size = 200.*m;
    G4Box* world_solid = new G4Box("WORLD", size/2., size/2., size/2.);
    
    G4LogicalVolume* world_logic = 
      new G4LogicalVolume(world_solid, G4Material::GetMaterial("G4_Galactic"),
			  "WORLD", 0, 0, 0, true);
    
    world_physi = 
      new G4PVPlacement(0, G4ThreeVector(), world_logic, "WORLD", 0, false, 0);
    
    world_logic->SetVisAttributes(G4VisAttributes::Invisible);
    
    
    // DETECTOR volumes definition ....................................
    _detector = new SciNearDetectorGeometry();
    
    
    G4LogicalVolume* detector_logic = _detector->GetLogicalVolume();
    
    
    // //SetMinimum Kinetic energy for tracking.
    //   G4double minKin = 100 * MeV;
    //   detector_logic->SetUserLimits(new G4UserLimits(DBL_MAX,DBL_MAX,DBL_MAX, minKin));
    
    G4PVPlacement* detector_physi = 
      new G4PVPlacement(0, G4ThreeVector(0,0,0), detector_logic, "DETECTOR", 
			world_logic, false, 0);
    
    // initiallize ROOT
    /****
	 TSystem ts;
	 
	 gSystem->Load("libMindClassesDict");
	 ROOT::Cintex::Cintex::Enable();
	 
	 const G4ElementTable* eltab = G4Element::GetElementTable();
	 const G4MaterialTable* mattab = G4Material::getMaterialTable();
	 
	 MindGeoTree* geotree = new MindGeoTree(world_physi, eltab, mattab);
	 
	 TFile* fo = new TFile("MindGeo.root","RECREATE");
	 
	 fo.WriteObject(geotree, "my_geo");
    ****/
    if ( _write_gdml ) _gdml.Write(_gdml_file_name , world_physi);
  }
  else if( _use_gdml ){
    _gdml.Read(_gdml_file_name);
    world_physi = dynamic_cast<G4PVPlacement*>(_gdml.GetWorldVolume());
    // SetNullField(*world_physi->GetLogicalVolume());
    regionmass["TASD"] = 0;
    regionmass["ACTIVE"] = 0;
    regionmass["PASSIVE"] = 0;
    G4String detectorName = "MIND/";
    _detectormass = world_physi->GetLogicalVolume()->GetMass();
    SetVolumeInformation(world_physi->GetLogicalVolume(), detectorName);
  }

  //myfile.open ("/data/neutrino05/phallsjo/copy/SaRoMan/example.txt");
  for ( int outerI = 0; outerI < world_physi->GetLogicalVolume()->GetNoDaughters(); outerI++ ) {
    
    G4int nDaughters = world_physi->GetLogicalVolume()->GetDaughter(outerI)->GetLogicalVolume()->GetNoDaughters();
    
    for ( int i = 0; i < nDaughters; i++ ) {
      G4VPhysicalVolume* outerDaughter = world_physi->GetLogicalVolume()->GetDaughter(outerI);
      G4VPhysicalVolume* daughter = world_physi->GetLogicalVolume()->GetDaughter(outerI)->GetLogicalVolume()->GetDaughter(i);
      
      G4LogicalVolume* myvol = daughter->GetLogicalVolume();
      G4Box* solidBox = (G4Box*) myvol->GetSolid();
      
      cerr<<outerDaughter->GetLogicalVolume()->GetName()<<endl;
      cerr<<outerDaughter->GetFrameTranslation()[2]<<endl;
      cerr<<myvol->GetName()<<endl;
      cerr<<daughter->GetFrameTranslation()[2]<<endl;

      _myfile<<myvol->GetName()<<" "
	//<<myvol->GetMass()/solidBox->GetCubicVolume()<<" "
	     <<outerDaughter->GetFrameTranslation()[0] + daughter->GetFrameTranslation()[0]<<" "
	     <<outerDaughter->GetFrameTranslation()[1] + daughter->GetFrameTranslation()[1]<<" "
	     <<outerDaughter->GetFrameTranslation()[2] + daughter->GetFrameTranslation()[2]<<" "
	     <<solidBox->GetXHalfLength()*2.0<<" "
	     <<solidBox->GetYHalfLength()*2.0<<" "
	     <<solidBox->GetZHalfLength()*2.0<<endl;
    }
  }
  
  return world_physi;
}

void MindDetectorConstruction::SetVolumeInformation(G4LogicalVolume* base, 
						    G4String detectorName) {

  // G4VPhysicalVolume* world = _parser.GetWorldVolume();
  G4int nDaughters = base->GetNoDaughters();
 
  for ( int i = 0; i < nDaughters; i++ ) {
    G4VPhysicalVolume* daughter = base->GetDaughter(i);
    G4LogicalVolume* myvol = daughter->GetLogicalVolume();
    G4String volumename = detectorName;
    volumename += myvol->GetName();
    volumename += "/";
    // first check if there are auxiliary objects
    std::cout<<volumename.c_str()<<std::endl;
    // If this name corresponds to one of the three regions add them to the map.
    if (myvol->GetName().contains("TASD")){
      regions["TASD"].push_back(daughter);
      regionmass["TASD"] += myvol->GetMass();
    } else if (myvol->GetName().contains("ACTIVE")){
      regions["ACTIVE"].push_back(daughter);
      regionmass["ACTIVE"] += myvol->GetMass();
    } else if (myvol->GetName().contains("PASSIVE")){
      regions["PASSIVE"].push_back(daughter);
      regionmass["PASSIVE"] += myvol->GetMass();
    } //Ignore otherwise
    
    const G4GDMLAuxListType auxlist = 
      _gdml.GetVolumeAuxiliaryInformation(myvol);

    // Squeak::mout(Squeak::info) << "Found volume " << myvol->GetName()
    // << " with "<< auxlist.size() <<" auxiliary elements." << std::endl;
    if (auxlist.size() > 0) {
      // Set auxiliary information
      SetAuxInformation(volumename, myvol, auxlist);
    }
    // 	 else {
    // 	 _detector->GetUserLimits().push_back(new G4UserLimits(_stepMax, _trackMax,
    // 	 _timeMax, _keThreshold));
    // 	 myvol->SetUserLimits(_detector->GetUserLimits().back());
    // 	 }
    if ( myvol->GetNoDaughters() > 0 ) {
      // Consider adding information to the daughter volumes
      SetVolumeInformation(myvol, volumename);
    }
  }
}

void MindDetectorConstruction::SetAuxInformation(G4String basename, 
						 G4LogicalVolume* myvol, 
						 const G4GDMLAuxListType auxlist){
  G4String sensdetname = basename;
  G4GDMLAuxListType::const_iterator vit = auxlist.begin();
  G4SDManager* SDMgr = G4SDManager::GetSDMpointer();
  do {
    try {
      if ((*vit).type.contains("SD")){
          G4VSensitiveDetector* mydet = SDMgr->FindSensitiveDetector((*vit).value);
          if ( mydet ){
              myvol->SetSensitiveDetector(mydet);
          } else {
              G4cout << sensdetname << " detector not found. Defining detector." << G4endl;
              sensDetList.push_back( new MindBarSD((*vit).value) );
              SDMgr->AddNewDetector(sensDetList.back());
              G4VSensitiveDetector* mydet = SDMgr->FindSensitiveDetector(sensdetname);
              myvol->SetSensitiveDetector(mydet);
          }
    }
    else if((*vit).type.contains("Color")){
	if((*vit).value.contains("Red"))
	  myvol->SetVisAttributes(new G4VisAttributes(G4Color(1, 0, 0)));

	else if((*vit).value.contains("Orange"))
	  myvol->SetVisAttributes(new G4VisAttributes(G4Color(1, 0.65, 0)));

	else if((*vit).value.contains("Blue"))
	  myvol->SetVisAttributes(new G4VisAttributes(G4Color(0.2, 0.2, 1)));

	else if((*vit).value.contains("Brown"))
	  myvol->SetVisAttributes(new G4VisAttributes(G4Color(0.61, 0.4, 0.12)));
      }
      else if ((*vit).type.contains("Invisible")) {
	myvol->SetVisAttributes(new G4VisAttributes(false));
      }
      else if ((*vit).type.contains("EMField")) {
	if ((*vit).value.contains("MagMap")){
	  SetMagneticField(*myvol);
	} else if((*vit).value.contains("NULL")){
	  SetNullField(*myvol);
	}
      }
    } catch (...) { 
      continue;
    }
    vit++;
  } while (vit != auxlist.end());	  
}


void MindDetectorConstruction::SetNullField(G4LogicalVolume& detector_logic)
{
  // apply a global uniform magnetic field along Y axis
  G4FieldManager* fieldMgr = new G4FieldManager();
  G4MagneticField* magField = 0;
  fieldMgr->SetDetectorField(magField);
  fieldMgr->CreateChordFinder(magField);
  fieldMgr->GetChordFinder()->SetDeltaChord(0.1*cm);
  
  detector_logic.SetFieldManager( fieldMgr, false );
}

void MindDetectorConstruction::SetMagneticField(G4LogicalVolume& vol) {
  
  G4FieldManager* fieldMgr = new G4FieldManager(); 
 
  const MindParamStore& config = 
    MindConfigService::Instance().Geometry();
   
  std::vector<G4double> fieldValue;  
  if ( config.PeekVParam("field") )
    fieldValue = config.GetVParam("field");

  double fieldScaling = +1.0;
  if(config.PeekDParam("FieldScaling"))
    fieldScaling = config.GetDParam("FieldScaling");
  if (config.PeekSParam("FieldMap")){
      // A field map has been provided.  Should add capability to
      // scale field -- but would such a scaling be physical???  As it
      // is now the "field" parameter is just a switch to turn the
      // field on and off.
      G4String Bmap = config.GetSParam("FieldMap");
      // Declaration of the magnetic field map object
      MindFieldMapR* magField = new MindFieldMapR(Bmap, fieldScaling, 30.0, 4, 30.0);
      // Now to embed the field in the detector geometry
      fieldMgr->SetDetectorField(magField);
      fieldMgr->CreateChordFinder(magField);
      fieldMgr->GetChordFinder()->SetDeltaChord(0.1*cm);
  } /*else {
      std::cout<<"Using uniform magnetic field."<<std::endl;
      G4UniformMagField* magField = 
	new G4UniformMagField(G4ThreeVector(fieldScaling*fieldValue[0]*tesla, 
					    fieldScaling*fieldValue[1]*tesla, 
					    fieldScaling*fieldValue[2]*tesla));
      fieldMgr->SetDetectorField(magField);
      fieldMgr->CreateChordFinder(magField);
      fieldMgr->GetChordFinder()->SetDeltaChord(0.1*cm);
      }*/
  vol.SetFieldManager( fieldMgr, true );
}

G4ThreeVector MindDetectorConstruction::GetVertex(G4String region_name){
  
  G4int nvols = regions[region_name].size();
  G4Box* solidBox;
  G4ThreeVector offset;
  if ( nvols > 1 ) {
    G4int volselect = int(floor(G4UniformRand() * double(nvols)));
    // Get the solid definition
    solidBox = (G4Box*) regions[region_name][volselect]->GetLogicalVolume()->GetSolid();
    offset = regions[region_name][volselect]->GetFrameTranslation();
  } else if (nvols == 1){
    solidBox = (G4Box*) regions[region_name][0]->GetLogicalVolume()->GetSolid();
    offset = regions[region_name][0]->GetFrameTranslation();
  } else {
    // give up
    return G4ThreeVector(0., 0., 0.);
  }
  double x = (2*G4UniformRand() - 1)*solidBox->GetXHalfLength();
  double y = (2*G4UniformRand() - 1)*solidBox->GetYHalfLength();
  double z = (2*G4UniformRand() - 1)*solidBox->GetZHalfLength();
  return G4ThreeVector(x + offset[0], y + offset[1], z + offset[2]);
  
}

G4String MindDetectorConstruction::GetRegion(G4ThreeVector vertex){
  bool inTASD=false;
  bool inACTIVE=false;
  bool inPASSIVE=false;
  for(int i=0; i<regions["TASD"].size(); i++){
    if (regions["TASD"][i]->GetLogicalVolume()->GetSolid()->Inside(vertex) == kInside)
      inTASD=true;
  }
  for(int i=0; i<regions["ACTIVE"].size(); i++){
    if (regions["ACTIVE"][i]->GetLogicalVolume()->GetSolid()->Inside(vertex) == kInside)
      inACTIVE=true;
  }
  for(int i=0; i<regions["PASSIVE"].size(); i++){
    if (regions["PASSIVE"][i]->GetLogicalVolume()->GetSolid()->Inside(vertex) == kInside)
      inPASSIVE=true;
  }
  G4String region="PASSIVE";
  if (inTASD){
    region = "TASD";
  }
  else if(inACTIVE){
    region = "ACTIVE";
  }
  else if(inPASSIVE){
    region = "PASSIVE";
  }
  return region;
}

/*
void MindDetectorConstruction::FindVolume(G4LogicalVolume* vol, G4String region_name,
					  std::vector<G4PhysicalVolume*> region_physi){
  
  for(G4int i=0; i < vol->GetNoDaughters(); i++){
    G4VPhysicalVolume* daughter = vol->GetDaughter(i);
    myvol = daughter->GetLogicalVolume();
    G4String volumename = myvol->GetName();
    if ( volumename.contains(region_name) ){
      // This is the volume we are looking for
      region_vol.push_back(daughter);
    }
    // Otherwise check if there are daughter volumes
    else if ( myvol->GetNoDaughters() > 0 ){
      // Use this function for the nex level down
      FindVolume(myvol, region_name, region_physi);
    }
  }
}
*/



/*
void MindDetectorConstruction::ConstructSDanField(){

  G4String SD_name = "/mind/ACTIVE/";
  MindSD* SD = new MindSD(SD_name);
  SetSensitiveDetector("ACTIVE", detector_logic);

}
*/

