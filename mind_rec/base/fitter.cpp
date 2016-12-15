#include <fitter.h>
#include <TMath.h>

using namespace bhep;///

//
class sortTrajByLength{
public:
  bool  operator()(const Trajectory* t1,const Trajectory* t2 ){
    if (t1->length() > t2->length()) return true;
    return false;
  }

};



//*************************************************************
fitter::fitter(const bhep::gstore& pstore,bhep::prlevel vlevel){
  //*************************************************************
  
  _level = vlevel;
  
  _store = pstore;
  
  _m = bhep::messenger(_level);

  _m.message("++fitter Messenger generated++",bhep::VERBOSE);


}

//*************************************************************
fitter::~fitter() {
  //************************************************************
}

//*************************************************************
void fitter::Initialize() {
  //*************************************************************
    
  _m.message("+++ fitter init  function ++++",bhep::VERBOSE);
  
  // initialize the vector<vector>
  _hadUnit.push_back(EVector(3,0)); _hadUnit.push_back(EVector(3,0));
  _hadUnit[0][2] = 1.; _hadUnit[1][2] = 1.;
  _detect = _store.fetch_sstore("detect");
  _testBeam = _store.fetch_istore("testBeam");

  // read parameters
  ReadParam();
  
  // initialize geometry
  _geom.init(_store, _level);

  //Instantiate recpack manager.
  MINDfitman::instance().set_man_parameters( _store, _geom.setup() );
  
  if (_X0 == 0){
    man().model_svc().enable_noiser(_model, RP::ms, false);
  }

  //If required make the clustering object.
  if ( _doClust )
    _clusters = new hit_clusterer( _store );

  ///initialize classifier
  get_classifier().Initialize( _store, _level, _geom.get_Fe_prop(), &_geom);

  _supergeom = _geom;

  _m.message("+++ End of init function ++++",bhep::VERBOSE);
}








//*************************************************************
bool fitter::Execute(bhep::particle& part,int evNo){
  //*************************************************************
  
  _m.message("+++ fitter execute function ++++",bhep::VERBOSE);
  
  bool ok = true;/// 
  _reseed_ok = false;
  _reseed_called = false;
  _fitted = false;
  _pr_count = 0;
  _muonindex[0] = 0; _muonindex[1] = 0;// The default muon track selection
  // _nonMuonEdep = 0 ;
  // _nonMuonHits = 0 ;
  
  ///create clusters or fill measurement vector
  ok = CreateMeasurements(part);
  
  //cout<<"What is patternRec: "<<_patternRec<<endl;

  ///if pattern recognition 
  if (_patternRec) {
    if ((int)_meas.size() < _min_seed_hits) {
      _failEvent = 1; 
      ok = false; 
    }
  } 
  else{ 
    ///for single traj
    Trajectory* straj = new Trajectory();
    ok = CreateSingleTrajectory(*straj);
    _trajs.push_back(straj);
  }
  
  if (!ok){
    _m.message("CreateMeasurements not ok", bhep::VERBOSE);
    _failEvent = 2; 
    return true;
  }
  
  //Sort in increasing z here when classifier up and running.!!!
  sort( _meas.begin(), _meas.end(), forwardSorter() );

  if(_testBeam)
    {
      const vector<bhep::hit*> hits = part.hits( _detect );

      //bool event_classif::Execute(const vector<cluster*>& hits,
      //		    vector<Trajectory*> &vtrajs, vector<cluster*>& hads) {
      
      //get_classifier().Execute( _meas, _trajs, _hadmeas);

      // fit the _meas clusters or 
      // const vector<bhep::hit*> hits = p.hits( _detect ); p= bhep::particle& part
      // What ever is the best. do line fit. add these hits to a track and return it.
      
      int fitcatcher;
      
      Trajectory* traj = new Trajectory();
      
     
      int nMeas = hits.size();
      double x[(const int)nMeas], y[(const int)nMeas], 
	z[(const int)nMeas], u[(const int)nMeas];
      int minindex = nMeas;
      
      for (int iMeas = 0;iMeas < nMeas;iMeas++){
	x[iMeas] = hits[iMeas]->x()[0];
	y[iMeas] = hits[iMeas]->x()[1];
	z[iMeas] = hits[iMeas]->x()[2];
      }
	
      TGraph *grY = new TGraph((const int)minindex, z, y);
      
      TF1 *funY = new TF1("liney","[0]+[1]*x",-4000,4000);
      funY->SetParameters(0.,0.001);
      
      fitcatcher = grY->Fit("liney", "QN");
      //cout<<"ConstantY="<<funY->GetParameter(0)<<endl;
      //cout<<"kY="<<funY->GetParameter(1)<<endl;
      //cout<<"x2Y="<<funY->GetChisquare()<<endl;

      TGraph *grX = new TGraph((const int)minindex, z, x);
      
      TF1 *funX = new TF1("linex","[0]+[1]*x",-4000,4000);
      funX->SetParameters(0.,0.001);
      
      fitcatcher = grX->Fit("linex", "QN");
      //cout<<"ConstantX="<<funX->GetParameter(0)<<endl;
      //cout<<"kX="<<funX->GetParameter(1)<<endl;
      //cout<<"x2X="<<funX->GetChisquare()<<endl;

      _xDir.clear();
      _yDir.clear();
      _x0.clear();
      _y0.clear();
      _xchi.clear();
      _ychi.clear();

      _xDir.push_back(funX->GetParameter(1));
      _yDir.push_back(funY->GetParameter(1));
      _x0.push_back(funX->GetParameter(0));
      _y0.push_back(funY->GetParameter(0));
      _xchi.push_back(funX->GetChisquare());
      _ychi.push_back(funY->GetChisquare());


      //Gets plane occupancies and total plane energies.
      //Needs hits in increasing z order.
      
      //   std::cout<<"+++++I am in get_plane_occupancy"<<std::endl;
      _m.message("++++ Calculating plane energies and occupancies ++++",bhep::VERBOSE);
      bool ok = true;
      
      std::vector<plane_info*> _planes;
      double _meanOcc;
      
      /// size of vector<cluster> hits
      size_t nHits = _meas.size();
      //cout<<"get_plane_occ ()  :: nHits="<<nHits<<endl;
      int single_count = 0;
      double testZ, testX, testY, curZ;
      size_t hits_used = 0, imeas = 0;
      int planeIndex=0 ;
      
      double _tolerance = _store.fetch_dstore("pos_resZ") * cm;
      
      /// loop over hits to calculate plane occupancy  
      do {
	
	testX = _meas[imeas]->position()[0];
	testY = _meas[imeas]->position()[1];
	testZ = _meas[imeas]->position()[2];
	hits_used++;
	// Avoid a hit with an undefined z position
	//if ( fabs(testZ) >  _detLength ) 
	//continue;
	// If nan.
	if(testX != testX) continue;
	if(testY != testY) continue;
	if(testZ != testZ) continue;
	
	///create plane info
	plane_info* plane = new plane_info(planeIndex, testZ, _store);
	plane->AddHit(_meas[imeas]);
	//cout<<"testZ: "<<testZ<<endl;
	
	///calculate the z position which is the current z for hits 1 -> total no of hits in the cluster 
	for (size_t i = hits_used;i <nHits;i++) {
	  curZ = _meas[i]->position()[2];
	  //cout<<"curZ: "<<curZ<<" testZ: "<<testZ<<" _tolerance: "<<_tolerance<<endl;
	  if (curZ <= testZ + _tolerance) {
	    
	    // add the hit to the same plane
	    plane->AddHit(_meas[i]);
	    //cout<<"Added hit"<<endl;
	    testZ = _meas[i]->position()[2];
	    hits_used++;
	  } else break;
	  
	}
	_m.message(" get plane info =",plane->GetZ()," Occ=",plane->GetNHits()," PlaneNo=",plane->GetPlaneNo(), bhep::VERBOSE);
	//cout<<"get plane info = "<<plane->GetZ()<<" Occ= "<<plane->GetNHits()<<" PlaneNo= "<<plane->GetPlaneNo()<<endl;
	///fill the plane_info vector
	_planes.push_back(plane);
	///increase the plane index
	planeIndex++;
	
	_meanOcc += (double) plane->GetNHits();
	
	
      } while (hits_used != nHits);
      
      ///total no of planes
      double _nplanes = (int)_planes.size();
      
      _meanOcc /= (double)_nplanes;

      // hits.size(); / _planes.size();

      //cout<<"_nplanes="<<_nplanes<<endl;

      //cout<<"_meanOcc="<<_meanOcc<<endl;

      _hitsPerPlanes.push_back(hits.size() / _planes.size());

      _avrHitsPerUsedPlanes.push_back(_meanOcc);

      /*
      std::vector<cluster*> meas;

      _clusters->execute( hits, meas ); 

      for(unsigned int cnt=0;cnt<meas.size();cnt++)
	{
	  RecObject* ro = dynamic_cast<RecObject*>(meas[cnt]);
	  traj->add_node(Node(*ro));
	}
      */

      //get_classifier().Execute( _meas, _trajs, _hadmeas);

      cout<<"planes="<<get_classifier().get_plane_info().size()<<endl;

      cout<<"free planes="<<get_classifier().get_free_planes()<<endl;

      //Root > Double_t chi2 = fit->GetChisquare();
      //Root > Double_t p1 = fit->GetParameter(1);
      //Root > Double_t e1 = fit->GetParError(1);
      
      //double qtilde = 0.3*pow(1+pow(fun->GetParameter(1),2),3./2.)/
      //(2*fun->GetParameter(2));	
      //_trajs.push_back(traj);
    } //end if(_testBeam)
  else
    {
  
      ///if pattern recognition and recTrajectory is ok
      if (_patternRec){
	
	/// execute event classification
	get_classifier().Execute( _meas, _trajs, _hadmeas);

	//cout<<"returned from event_classif.cpp _trajs.size ="<<_trajs.size()<<endl;
	/*
	  for(int i=0;i<_trajs.size();i++)
	  {
	  cout<<"after class trajs[i]->size() "<<_trajs[i]->size()<<endl;
	  cout<<"lowpt "<<_trajs[i]->quality("lowPt")<<endl;
	  }
	*/
	///sort the hadrons
	sort( _hadmeas.begin(), _hadmeas.end(), reverseSorter() );
	
	///PR seeds for all the trajectories from classifier
	_vPR_seed = get_classifier().get_patRec_seed_vector();
	
	//cout<<"Printing _vPR_seed"<<endl;
	//cout<<_vPR_seed[0]<<endl;
	
      }
    }

  //return true;

  /// for non PR track need to set tracks infos separately
  /*else if(ok) _trajs.push_back(_traj);*/
  double maxlength = -99999;
  double maxPlanes = -99999;

  /// loop over trajectories 
  for (unsigned int i=0; i<_trajs.size(); i++){ 
 
    _m.message("inside traj loop::if classifier ok, traj no =",i,"*********",bhep::DETAILED); 
    /*
    cout<<"traj no="<<i<<endl;
    cout<<"lowPt="<<_trajs[i]->quality("lowPt")<<endl;
    cout<<"trajs[i] size="<<_trajs[i]->size()<<endl;
    */
    // If the track is fitted
    //if(_trajs[i]->nodes()[0]->status("fitted"))
    // if(_trajs[i]->quality("fitted"))
    //if((_trajs[i]->first_fitted_node()).status("fitted"))
    
    //if(_trajs[i]->nodes()[_trajs[i]->first_fitted_node()]->status("fitted"))
    if(_trajs[i]->quality("lowPt") == 1)
      {
	//EVector vertex = EVector(3,0);
	//vertex = _trajs[i]->state(_trajs[i]->first_fitted_node()).vector();
	
	//cout<<"We have found a fitted track"<<endl;
	_fitted = true;
	ok = false; // No need to fit the track if it is fitted!
	_fitCheck = 1;// ?? Can not set this, thus the track is not added to fitted nodes. Must check this again.
	_reseed_ok = false;
	_reseed_called = false;

	_traj = *(_trajs[i]);

	//_traj.set_quality("lastIso", (int) _traj.size());
	_traj.set_quality("fitcheck", _fitCheck);
	_traj.set_quality("fitted",_fitted);
	_traj.set_quality("lowPt",1);
	_traj.set_quality("reseed",0);
	//cout<<"initialqP: "<<(double)(_traj.quality("initialqP"))<<endl;
	//cout<<"Momentum: "<<_traj.node(_traj.first_fitted_node()).state().vector()[5]<<endl;

	*(_trajs[i]) = _traj;
      }
    else
      {
    
	///
	_fitted = false;
	_reseed_ok = false;
	_reseed_called = false;
	_failType = 0;
	_intType = 0; //set to 'success' before run to avoid faults in value.
	ok = true;///track finded by PR or CA ??
	
	/// Get the trajectory
	_traj = *(_trajs[i]);
	//_traj.set_quality("lastIso", (int) _traj.size());
	_m.message("fitter::vector_PR size = ", _vPR_seed.size()," & trajno=",i,"  nmeas =",_traj.size(),bhep::DETAILED);
	
	//get traj informations
	int nplanes = 0, freeplanes = 0;
	double xtent = 0, vertZ =0;

	//Get traj info before fitting
	nplanes = (int)(_traj.quality("nplanes"));
	freeplanes = (int)(_traj.quality("freeplanes"));
	_intType = (int)(_traj.quality("intType"));
	_failType = (int)(_traj.quality("failType"));
	xtent = (double)(_traj.quality("xtent"));
	vertZ = (double)(_traj.quality("vertZ"));
	
	// _m.message("in fitter:: from classifier the traj =",_traj,bhep::DETAILED);
	_m.message("in fitter: for traj no =",i,"  intType =",_intType,"  failType=",_failType,bhep::VERBOSE);
	
	///sort the nodes in increasing Z (event when PR is not running) 
	_traj.sort_nodes(RP::z, 1);
	
	///if the traj finding fails during event_classification CA/PR anyone 
	if(_failType==4 || _failType==5 || _failType==6) ok = false;
	
	///track found by finder (CA or PR not failed)
	State seedState;
	if (ok) {
	  //SetFit mode to manager.
	  MINDfitman::instance().fit_mode();
	  
	  ok = CheckValidTraj(_traj);
	}
	
	///fit the trajectory 
	if (ok && _vPR_seed.size() > 0) {
	  /// seed for Fit 
	  ComputeSeed(_traj,seedState);
	  
	  _m.message("- traj node0=",*(_traj.nodes()[0]),bhep::DETAILED);

	  _m.message("- traj size=",_traj.nodes().size(),bhep::DETAILED);

	  // if (ok)cout<<"if classifier3="<<endl; 

	  cout<<"before fitTrajectory"<<endl;
	  for(int j = 0; j<_traj.size();j++)
	    {
	      if(_traj.nodes()[j]->status("fitted"))
		{
		  cout<<"momentum: "<<1.0/_traj.nodes()[j]->state().hv().vector()[5]<<endl;
		  cout<<_traj.nodes()[j]->measurement().position()[2]<<endl;
		  cout<<_geom.getBField(_traj.nodes()[j]->measurement().position())[0]<<endl;
		}
	      else break;
	    }

	  _fitted = FitTrajectory(seedState,i);
	  
	  _m.message("- traj node0=",*(_traj.nodes()[0]),bhep::DETAILED);
	  
	  _m.message("- copied trajectory =", _traj,bhep::DETAILED);
	  
	  //
	  //if(ok && _fitted && abs(_length) > maxlength){
	  if(_fitted && abs(_length) > maxlength){
	    maxlength = abs(_length);
	    _muonindex[0] = i;
	  }
	  if(nplanes > maxPlanes){
	    maxPlanes = double(nplanes);
	    _muonindex[1] = i;
	  }
	}

	cout<<"End of fitter.cpp Fitted: "<<_fitted<<endl;
	cout<<"End of fitter.cpp forwardFitcheck: "<<_forwardFitCheck<<endl;
	cout<<"End of fitter.cpp reseedFitcheck: "<<_reseedFitCheck<<endl;
	cout<<"End of fitter.cpp initialqP: "<<_initialqP<<endl;

    for(unsigned int cnt=0;cnt<_trajs.size();cnt++)
      {
	cout<<"end trajs[i]->size() "<<_trajs[cnt]->size()<<endl;
	cout<<"lowpt "<<_trajs[cnt]->quality("lowPt")<<endl;
      }

	//cout<<"End of fitter.cpp momentum: "<<_traj.state().hv().vector()[5]<<endl;
	if(_fitted && _fitCheck){

	  for(int j = 0; j<_traj.size();j++)
	    {
	      if(_traj.nodes()[j]->status("fitted"))
		{
		  cout<<"momentum: "<<1.0/_traj.nodes()[j]->state().hv().vector()[5]<<endl;
		  cout<<_traj.nodes()[j]->measurement().position()[2]<<endl;
		  cout<<_geom.getBField(_traj.nodes()[j]->measurement().position())[0]<<endl;
		}
	      else break;
	    }


	  
	  cout<<"End of fitter.cpp momentum: "<<1.0/_traj.node(_traj.first_fitted_node()).state().hv().vector()[5]<<endl;
	  /*
	  cout<<_geom.getBField(_traj.node(_traj.first_fitted_node()).measurement().position())[0]<<endl;
	  cout<<"End of fitter.cpp momentum: "<<1.0/_traj.node(_traj.last_fitted_node()).state().hv().vector()[5]<<endl;
	  cout<<_geom.getBField(_traj.node(_traj.last_fitted_node()).measurement().position())[0]<<endl;
	  */
	}
	///assign quality for each trajectory
	_traj.set_quality("failType",_failType);
	_traj.set_quality("intType",_intType);
	_traj.set_quality("nplanes",nplanes);
	_traj.set_quality("freeplanes",freeplanes);
	_traj.set_quality("reseed",_reseed_ok);
	_traj.set_quality("xtent",xtent);
	_traj.set_quality("initialqP",_initialqP);
	_traj.set_quality("fitted",_fitted);
	_traj.set_quality("vertZ", vertZ);
	_traj.set_quality("fitcheck", _fitCheck);
	_traj.set_quality("lowPt",0);
	
	*(_trajs[i]) = _traj;
	// cout<<"traj ="<<i<<" failType ="<<_failType<<" fitted ="<<(int)traj.quality("fitted")<<endl;;  
	
	// std::cout<<"*********************************************"<<std::endl;
	//std::cout<<" ++++++ Trajectory "<<i<<"  fitted+++++++++++"<<std::endl;
	//std::cout<<"*********************************************\n"<<std::endl;
      } // End loop over trajectories.
    
    /*
    for(int i=0;i<_trajs.size();i++)
      {
	cout<<"final trajs[i]->size() "<<_trajs[i]->size()<<endl;
	cout<<"lowpt "<<_trajs[i]->quality("lowPt")<<endl;
      }
    */
  }
    /// for hadron shower
    if((int)_trajs.size() != 0)
      for(int j=0; j<=(int)_trajs.size(); j++)
	rec_had_edep(j);
    else 
      rec_had_edep(0);
    // rec_had_energy();

    /*
    for(int i=0;i<_trajs.size();i++)
      {
	cout<<"after rec_had_edep trajs[i]->size() "<<_trajs[i]->size()<<endl;
	cout<<"lowpt "<<_trajs[i]->quality("lowPt")<<endl;
      }
    */

    //}  
  
    //std::cout<<"Final trajectory =" << _traj<<std::endl;
  
  _m.message(" I am ******************fitter end",bhep::VERBOSE); 
  
  //return _fitted;
  return true;///signifies fitter executed
}

//*************************************************************
void fitter::Reset() {
  //*************************************************************
  
  _m.message("+++ Reset function +++",bhep::VERBOSE);

  //Reset trajectory   
  _hadmeas.clear();
  _failEvent = 0;///
  _pr_count = 0;///
  _nonMuonEdep.clear();
  _traj.clear();
  _traj2.clear();
  _traj3.clear();
  _showerDir.clear();
  _showerVertZ.clear();
  _showerNplanes.clear();
  _showerXtent.clear();
  _vPR_seed.clear();
  stc_tools::destroy(_meas);
  stc_tools::destroy(_trajs);
  // _hadUnit.clear();
}

//*************************************************************
bool fitter::FitTrajectory(const State& seedState0, const int trajno) {
  //*************************************************************

  _m.message("+++ FitTrajectory function ++++",bhep::VERBOSE);

  // Trajectory traj1 = _traj;
  //_traj2 = _traj;

  bool ok = true;
  bool ok0, ok_quality;
  bool ok1 = false;

  // Experimenting with removing hits in a different magnetic field region.
  //_traj.nodes().erase(_traj.nodes().begin(),_traj.nodes().begin()+4);

  //cout<<"Removing nodes"<<endl;
  //cout<<"_traj.size(): "<<_traj.size()<<endl;
  /*
    for(int i = _traj.size()-1; i>=0; i--)
    {
    if(_geom.getBField(_traj.node(i).measurement().position())[0] > 0)
    {
    cout<<"Removing node: "<<i<<endl;
    _traj.nodes().erase(_traj.nodes().begin()+i) ;
    }
    }
  */
  //cout<<"Removing nodes done"<<endl;
  
  //if(_traj.size()>=16){
  //  _traj.nodes().erase(_traj.nodes().end()-1,_traj.nodes().end());
  //}
  
  //_traj.nodes().erase(_traj.nodes().end()-3,_traj.nodes().end());
  //_traj.nodes().erase(_traj.nodes().begin(),_traj.nodes().begin()+1);
  
  //cout<<"_traj.size(): "<<_traj.size()<<endl;
  //cout<<"Removing nodes done"<<endl;
  
  //_traj.sort_nodes(RP::z, 1);

  _traj2 = _traj;

  //  _traj.nodes()[0]->reset();
  //_traj.nodes()[1]->reset();
  //_traj.nodes()[2]->reset();
  //_traj.nodes()[3]->reset();

  /// fit the trajectory                
  ok0 = man().fitting_svc().fit(seedState0, _traj);
  //ok0 = man().fitting_svc().fit(seedState0, _traj, "particle/helix");
  //ok0 = man().fitting_svc().fit(seedState0, _traj, "kalman");
  
  // Check the quality if the traj is fitted
  if(ok0) ok_quality = CheckQuality(_traj); 
  
  if(!ok_quality)
    {
      cout<<"bad quality"<<endl;
    }
  
  ///refit the trajectory only when the quality is not good
  if (_refit && ok0 && !ok_quality){    
    State seedState1;
    ComputeSeedRefit(_traj, seedState1);
    ok1 = man().fitting_svc().fit(seedState1,_traj);
  }
  else
    {
      ok1=true;
    }

  ///check number of fitted nodes in traj
  _fitCheck =0;
  vector<Node*>::iterator nDIt;
  for (nDIt = _traj.nodes().begin();nDIt!=_traj.nodes().end();nDIt++)
    {
      /*
      cout<<"_traj.nodes Fitted: "<<(*nDIt)->status("fitted")<<" "
	  <<(*nDIt)->measurement().position()[2]<<endl;
      */
    if ( (*nDIt)->status("fitted") )
      _fitCheck++;

    }
  _forwardFitCheck = _fitCheck;
  
  //cout<<"fitCheck="<<_fitCheck<<" for trajectory size "<<_traj.size()<<endl;
  
  ///check for reseeding 
  double low_fit_cut;
  if (_intType == 2)
    low_fit_cut = _lowFit2;
  else
    low_fit_cut = _lowFit1;
  /*
  cout<<"intType="<<_intType<<endl;
  cout<<"ok0="<<ok0<<endl;
  cout<<"ok1="<<ok1<<endl;
  cout<<"(double)_fitCheck/(double)_traj.size()="<<(double)_fitCheck/(double)_traj.size()<<endl;
  cout<<"low_fit_cut="<< low_fit_cut<<endl;
  */

  ///reseed the trajectory
  if (_intType!= 5){   // not CA
    if (_intType != 2){ // not all planes are single occ  
      if ( ( !ok0 || !ok1 || (double)_fitCheck/(double)_traj.size() < low_fit_cut ) && _fitCheck > 0)
	// _reseed_ok = ReseedTrajectory(_traj2, trajno);
	_reseed_ok = ReseedTrajectory(trajno);
    } 
    else if (((double)_fitCheck/(double)_traj.size() < low_fit_cut ) && _fitCheck > 0){      
      ///if traj contains all single occ planes
      // _reseed_ok = ReseedTrajectory(_traj2,trajno);      
      //_reseed_ok = ReseedTrajectory(trajno); 
      _reseed_ok=true;
    }
    else _pr_count++;    
  }

  //if reseed successful
  if (_reseed_ok){    
    _traj = _traj2;
  }
  else if (!ok0) ok=false;

  // std::cout<<"copied trajectory =" << _traj<<std::endl;
  // std::cout<<"Trajectory status is "<<_traj.status("fitted")<<std::endl;
  // std::cout<<"First fitted node is "<<_traj.first_fitted_node()<<std::endl;

  if(ok) _m.message("***inside FitTrajectory, ok=", bhep::VERBOSE);
  _fitCheck = 0;
  for (nDIt = _traj.nodes().begin();nDIt!=_traj.nodes().end();nDIt++)
    {/*
      cout<<"_traj.nodes Fitted reseed: "<<(*nDIt)->status("fitted")<<" "
	  <<(*nDIt)->measurement().position()[2]<<endl;
     */
    if ( (*nDIt)->status("fitted") )
      _fitCheck++;

    }

  _reseedFitCheck = _fitCheck;
  
  //cout<<"fitCheck="<<_fitCheck<<" for trajectory size in refit "<<_traj.size()<<endl;
  
  ///length of the traj
  if(_traj.status(RP::fitted) && _fitCheck > 0){
	
    ok = man().matching_svc().compute_length(_traj, _length);///
    //man().matching_svc().compute_length(_traj, _length);
    //std::cout<<"_traj length = "<<_length<<std::endl;
  }

  //cout<<"Traj nodes size= "<<_traj.nodes().size()<<endl;

  //cout<<"nodes"<<endl;
  /*
  for(int i=0;i<_traj.nodes().size();i++)
    {
      //cout<<*(_traj.nodes()[i])<<endl;
      cout<<_traj.nodes()[i]->measurement().position()[2]<<endl;
      cout<<_traj.nodes()[i]->status("fitted")<<endl;
    }
  */
  //cout<<"Done"<<endl;


  //std::cout<<"Final trajectory =" << _traj<<std::endl;
  // traj.set_status("fitted", ok);

  return ok; 
  
}
  
//*************************************************************
bool fitter::ReseedTrajectory(const int trajno){
  //*************************************************************
  _m.message("****inside ReseedTrajectory *****************traj.nmeas = ",_traj2.size(),"\n", bhep::VERBOSE);
  
  bool ok, ok1;
  _reseed_called = true;
  
  //State backSeed;  

  //State backSeed = get_classifier().get_patRec_seed();

  State backSeed = _vPR_seed[_pr_count];
  
  cout<<"Printing backSeed"<<endl;
  cout<<backSeed<<endl;

  /*
  //Want to re-seed with diagonal matrix.
  HyperVector HV1 = backSeed.hv();//.keepDiagonalMatrix();
  HV1.keepDiagonalMatrix();
  backSeed.set_hv( HV1 );
  */
  cout<<"z pos in fitter_reseed: "<<backSeed.vector()[2]<<endl; 
    
  ///sort nodes in reverse order
  _traj2.sort_nodes(RP::z, -1);
  //_traj2.sort_nodes(RP::z, 1);

  cout<<"z pos in track:"<< _traj2.nodes()[0]->measurement().vector()[2]<<endl;

  //ComputeSeed(_traj2,backSeed);

  //(void) get_classifier().get_patternRec_seed(backSeed,_traj2);


  // a copy of the track
  // Trajectory traj1 = traj;
  _traj3 = _traj2;

  ///fit the traj
  try  {
    ok = man().fitting_svc().fit(backSeed,_traj2);
  } catch (const char* msg ) {
    ok = false;
  }

  /// compute seed for refitting and refit the traj
  if (ok && _refit){
    // Check the quality if the traj is fitted
    if (!CheckQuality(_traj2)){ 
      State seedState1;
      ComputeSeedRefit(_traj2,seedState1);
      try {
	//ok1 = man().fitting_svc().fit(seedState1,_traj3);
	ok1 = man().fitting_svc().fit(backSeed,_traj3);
      } catch (const char* msg ) {
	ok1 = false;
      }
    }
  }

  if (ok1){
    _traj2 = _traj3;
    ok=true;
  }

  // sort nodes back

  if(ok1){
  _traj2.sort_nodes(RP::z, 1);
  }

  // _traj2 = traj;
  ///increase the count to set backseed from PR vector
  _pr_count++;
  
  return ok;
}

//*************************************************************
void fitter::ComputeSeedRefit(const Trajectory& traj, State& seedState) {
  //*************************************************************
  
  _m.message("Going to calculate seed for refit...",bhep::VERBOSE);
  
  //--------- refit using a new seed --------//
  /// what is the differance in this new state than the earlier seed??	
  seedState = traj.state(traj.first_fitted_node());
  
  EVector v = seedState.vector();
  EMatrix C0 = seedState.matrix();
  
  ApplyCovarianceFactor(_facRef,C0);
  
  HyperVector HV(v,C0,RP::slopes_curv_z);
  HV.keepDiagonalMatrix();
  
  /// seedstate.set_hv( HV );
  seedState.set_hv( HV );
}

//*************************************************************
void fitter::rec_had_energy(){
  //*************************************************************
 
  // loop over all of the tracks. Ignore potential muon tracks.
  double dxdz=0, dydz=0;
  double hadEdep=0;
  
  if(_trajs.size() > 1){
    // for(int j=0; j<2; j++){
    if(1) { 
      int j = 0;
      dxdz=0, dydz=0;
      hadEdep=0;
      _transEdx[j] = 0;

      for(int i=0; i<(int)_trajs.size(); i++){
	if(i == _muonindex[j]) continue;

	// first sum the directions of the non-muon tracks
	Trajectory& traj = *_trajs[i];///
	_fitCheck = 0;
	vector<Node*>::iterator nDIt;
	for (nDIt = traj.nodes().begin();nDIt!=traj.nodes().end();nDIt++)
	  if ( (*nDIt)->status("fitted") )
	    _fitCheck++;

	if(traj.quality("fitted") && _fitCheck > 0){
	  State currstate = traj.state(traj.first_fitted_node());
	  EVector v = currstate.vector();
	  dxdz += v[3];
	  dydz += v[4];
	}
      }
      _hadEng[j] = 0; // hadEdep * ( 1 + widthI/widthS/rel_dedx/rel_dens);
      double hnorm = sqrt(1 + pow(dxdz,2) + pow(dydz,2));
      _hadUnit[j][0] = dxdz/hnorm; 
      _hadUnit[j][1] = dydz/hnorm;
      _hadUnit[j][2] = 1./hnorm;

      // now repeat the loop over trajectories to sum over energy deposition
      // create a transverse energy profile
      double vertz = 0;
      if(_trajs[_muonindex[j]]->quality("vertZ"))
	vertz = _trajs[_muonindex[j]]->quality("vertZ");
      for(int i=0; i<(int)_trajs.size(); i++){
	if(i == _muonindex[j]) continue;
	// first sum the directions of the non-muon tracks
	Trajectory& traj = *_trajs[i];
	// Add together the energies of all of the hits in the track
	std::vector<Node*> hits = traj.nodes();
	EVector zunit = EVector(3,0); zunit[2] = 1;
	EVector uxz = crossprod(_hadUnit[j],zunit);
	for ( int iHit=0; iHit < (int)hits.size(); iHit++){
	  EVector hpos = EVector(3,0);
	  hpos[0] = hits[iHit]->measurement().position()[0];
	  hpos[1] = hits[iHit]->measurement().position()[1];
	  hpos[2] = hits[iHit]->measurement().position()[2] - vertz;
	  
	  hadEdep += get_classifier().correctEdep((hits[iHit]->measurement().hv("energy").vector()[0])*MeV,hpos[0],hpos[1],hpos[2]);
	  EVector pvec = EVector(3,0);
	  pvec[0] = dxdz * hpos[2];
	  pvec[1] = dydz * hpos[2];
	  pvec[2] = hpos[2];
	  double posdiff = dot(hpos - pvec,uxz)/hpos[2];
	  _transEdx[j] += get_classifier().correctEdep((hits[iHit]->measurement().hv("energy").vector()[0])*MeV,hpos[0],hpos[1],hpos[2]) * posdiff;
	}
      }
      
      // else if(_trajs.size() == 1){ // assume there is one muon track and "hidden" hadronization
      //  _traj = _trajs[0];
      
      if(_hadmeas.size() > 2){
	for(int iHHit=0;iHHit<(int)_hadmeas.size();iHHit++){
	  hadEdep += (_hadmeas[iHHit]->get_eng())*GeV;
	}
      }
      // Inflated definition assuming a continuous, uniform energy loss through detector
      _hadEng[j] = hadEdep;
      
    }
  } 
}

//*************************************************************
void fitter::rec_had_edep(int j){
  //*************************************************************
  double hadEdep = 0;
  EVector hadCentroid = EVector(3,0);
  std::vector<cluster*> hadHits;
  // create a new trajectory based on the hadron hits alone.
  
  double minZ=999999.9, maxZ=-99999.9;
  
  _nonMuonEdep.push_back(0.0);
  _nonMuonHits.push_back(0);
  _showerDir.push_back(EVector(3,0));
  _showerVertZ.push_back(0.0);
  _showerNplanes.push_back(0);
  _showerXtent.push_back(0.0);

  if((int)_trajs.size() != 0 && j != (int)_trajs.size()) {

    //sort the trajectories
    sort( _trajs.begin(), _trajs.end(), sortTrajByLength());

    //loop over trajectories
    EVector hadCentroid = EVector(3,0);
    for(int i=0; i<(int)_trajs.size(); i++){
      
      hadHits = _hadmeas;
      // cout<<i<<" _muonindex[i] ="<< _muonindex[i]<<endl; 
      // if(i == _muonindex[i]) continue;
      if(i == j) continue;
      
      Trajectory& traj = *_trajs[i];
      
      std::vector<Node*> hits = traj.nodes();
      
      //Add total non-muon hits
      _nonMuonHits[j] += hits.size();
      
      for ( int iHit=0; iHit < (int)hits.size(); iHit++){
	EVector hpos = EVector(3,0);
	hpos[0] = hits[iHit]->measurement().position()[0];
	hpos[1] = hits[iHit]->measurement().position()[1];
	hpos[2] = hits[iHit]->measurement().position()[2] ;
	
	if(hpos[2] < minZ) minZ = hpos[2];
	if(hpos[2] > maxZ) maxZ = hpos[2];
	hadEdep += get_classifier().correctEdep((hits[iHit]->measurement().hv("energy").vector()[0])*MeV,hpos[0],hpos[1],hpos[2]);
	
	//sum of position to calculate centroid
	hadCentroid[0] += hpos[0];
	hadCentroid[1] += hpos[1];
	hadCentroid[2] += hpos[2];
	
      }
      _m.message("i=",i," hits =",hits.size(),"  _nonMuonHits=",_nonMuonHits[j]," hadEdep = ",hadEdep," hadHits=",hadHits.size(),bhep::VERBOSE);
    }
  }
  else
    hadHits = _meas;
 
  // total edep for non-muon hits
  _nonMuonEdep[j] += hadEdep ;
  _nonMuonHits[j] += hadHits.size();
  _m.message(" hadHits=",hadHits.size(),bhep::VERBOSE); 
      
  if(hadHits.size() !=0){
    for(int ih=0;ih<(int)hadHits.size();ih++){
      RecObject* ro = dynamic_cast<RecObject*>(hadHits[ih]);
      _hadTrajs.add_node(Node(*ro));
      EVector hadPos = EVector(3,0);
      hadPos[0] = hadHits[ih]->position()[0];
      hadPos[1] = hadHits[ih]->position()[1];
      hadPos[2] = hadHits[ih]->position()[2];
      
      if(hadPos[2] < minZ) minZ = hadPos[2];
      if(hadPos[2] > maxZ) maxZ = hadPos[2];
      //edep
      _nonMuonEdep[j] += get_classifier().correctEdep((hadHits[ih]->get_eng())*MeV,hadPos[0],hadPos[1],hadPos[2]) ;
      
      //sum of position to calculate centroid
      hadCentroid[0] += hadPos[0];
      hadCentroid[1] += hadPos[1];
      hadCentroid[2] += hadPos[2];
    }
  }
  // else hadCentroid = EVector(3,0);
  
  //calculate the centroid
  if (_nonMuonHits[j]) hadCentroid /= _nonMuonHits[j];
  else hadCentroid = EVector(3,0);
  
  _m.message(" hadEdep= ",hadEdep,"  nonMuonEdep= ",_nonMuonEdep[j],"  nonMuonHits= ",_nonMuonHits[j],bhep::VERBOSE);
  
  //direction w.r.t reconstructed vertex
  double dxdz=0, dydz=0, dz = 0, norm = 0;
  EVector vertex = EVector(3,0);
  
  // reconstructed vertex
  if(_trajs.size() != 0 
     && _trajs[0]->quality("fitcheck") > 0 && _trajs[0]->quality("fitted")) 
    // Use the first node of the longest trajectory.
    vertex = _trajs[0]->state(_trajs[0]->first_fitted_node()).vector();
  else
    // Use the first node in the event.
    vertex = _meas[0]->vector(); 
  
  
  _m.message("vertex   =",vertex[0],"   ",vertex[1],"  ",vertex[2], bhep::VERBOSE);
  _m.message("centroid =",hadCentroid[0],"   ",hadCentroid[1],"  ",hadCentroid[2],bhep::VERBOSE);
  
  //direction
  dz   = (hadCentroid[2] -vertex[2]);
  dxdz = (hadCentroid[0] -vertex[0])/dz;
  dydz = (hadCentroid[1] -vertex[1])/dz;
  
  //unit direction
  norm = sqrt(1 + pow(dxdz,2) + pow(dydz,2));
  
    
  if(norm != 0) {
    _showerDir[j][0] = dxdz/norm; 
    _showerDir[j][1] = dydz/norm;
    _showerDir[j][2] = 1./norm;
  }
  else _showerDir[j] = EVector(3,0);
  
  _showerVertZ[j]   = minZ;
  _showerNplanes[j] = int((maxZ - minZ)/(_geom.getPieceWidth()));
  _showerXtent[j]   = maxZ - minZ;

  // Can only add the hadron trajectory once
  
  
  _hadTrajs.set_quality("failType",10);
  _hadTrajs.set_quality("intType",10);
  _hadTrajs.set_quality("nplanes",int((maxZ - minZ)/(_geom.getPieceWidth())));
  _hadTrajs.set_quality("freeplanes",0);
  _hadTrajs.set_quality("reseed",0);
  _hadTrajs.set_quality("xtent",maxZ - minZ);
  _hadTrajs.set_quality("initialqP", 0.0);
  _hadTrajs.set_quality("fitted", 0);
  _hadTrajs.set_quality("vertZ", minZ);
  _hadTrajs.set_quality("fitcheck", 0);
    
  // trajs.push_back(hadTraj);
  _m.message("Rec hadron Unit direction components:", _showerDir[j][0], _showerDir[j][1], _showerDir[j][2],bhep::VERBOSE);
  
}


//*************************************************************
bool fitter::CheckQuality(const Trajectory& traj){
  //*************************************************************
    
 _m.message("+++CheckQuality++++",bhep::VERBOSE);

  bool ok = true;
    
  if (traj.quality()>_chi2fit_max) ok=false;

  cout<<"traj.quality()="<<traj.quality()<<endl;
  cout<<"_chi2fit_max="<<_chi2fit_max<<endl;
       
  return ok;

}

//*************************************************************
bool fitter::CreateMeasurements(const bhep::particle& p) {
  //*************************************************************
 
  _m.message("+++ CreateMeasurements function ++++",bhep::VERBOSE);
  
  Reset();
  bool ok = true;

  //string detect = _store.fetch_sstore("detect");
  const vector<bhep::hit*> hits = p.hits( _detect ); 
  //std::cout<<"Hits size creat meas= "<<hits.size()<<std::endl;
  //Cluster or directly make measurements.
  if ( _doClust && hits.size() != 0 ){
    // Make clusters
    _clusters->execute( hits, _meas ); 
  }
  else {
    // Create a cluster of each hit (without clustering)
    for(size_t j=0; j< hits.size(); j++){
      //---------- create measurement ---------------//
      cluster* mnt = GetMeasurement(*hits[j]);
      
      _meas.push_back(mnt); 
      
      _m.message("Measurement added:",*mnt,bhep::VVERBOSE);
    }//end of loop over hits
  }
  return ok;
}

//*************************************************************
bool fitter::CreateSingleTrajectory(Trajectory& traj) {
  //*************************************************************
 
  _m.message("+++ CreateSingleTrajectory function ++++",bhep::VERBOSE);
  //--------- add measurements to trajectory --------//
  ///create the trajectory  
  std::vector<cluster*>::iterator it1;
  for (it1 = _meas.begin();it1 != _meas.end();it1++){
    RecObject* ro = new RecObject();
    ro = dynamic_cast<RecObject*>(*it1);
    Node temp;
    temp.set_measurement(*ro);
    traj.add_node(temp);
    // measurement( *(*it1) );
  }
  _m.message("Trajectory created:",traj,bhep::VVERBOSE);

  return true;
}


//*************************************************************
bool fitter::CheckValidTraj(const Trajectory& traj) {
  //*************************************************************
  //cout<<" +++++inside fitter:: CheckValidTraj func "<<endl;  
  bool returnBool = true;

  //--------- Reject too many hits --------//
  
  if ((int)traj.size() < _lowPass) { 
    _failType = 1;
    returnBool = false;
  }
  return returnBool;
}

//*****************************************************************************
double fitf(Double_t *x,Double_t *par) { 
  //*****************************************************************************

  double z = x[0]; 
  double fitval = par[0]+par[1]*z+par[2]*z*z;

  return fitval ;
}

//*************************************************************
int fitter::GetQ(const Trajectory& traj){
  //*************************************************************
  double q = 0;
  
  if (_model.compare("particle/helix")==0)
    { 
      q = traj.state(traj.last_fitted_node()).vector()[dim-1];
      
      if (q<0) q=-1; else q=1;     
    }
  return (int) q;
}

//*************************************************************
cluster*  fitter::GetMeasurement(bhep::hit& hit){
  //*************************************************************
    
  _m.message("+++ getMeasurement function ++++",bhep::VERBOSE);
    
  //---- generate a virtual plane to hold the hit ----//
    
  bhep::Point3D bhit_pos = hit.x(); 

  string meastype = _geom.getMeasType();
  EMatrix cov = _geom.getCov();
  //pnumber++;

  //----- generate repack hit from bhep one ----//
    
  EVector hit_pos(2,0);
  hit_pos[0]=bhit_pos[0];
  hit_pos[1]=bhit_pos[1];

  EVector meas_pos(3,0);
  meas_pos[0] = hit_pos[0];
  meas_pos[1] = hit_pos[1];
  meas_pos[2] = bhit_pos[2];
    
  cluster* me = new cluster();
  me->set_name(meastype);
  me->set_hv(HyperVector(hit_pos,cov,RP::xyz));
  me->set_name("volume", "mother");
  me->set_position( meas_pos );
  //Add the hit energy deposit as a key to the Measurement.
  const dict::Key Edep = "E_dep";
  const dict::Key EdepVal = bhep::to_string( hit.ddata("TotalEng") );
  me->set_name(Edep, EdepVal);
  if (_patternRec){
    const dict::Key motherP = "MotherParticle";
    //const dict::Key mothName = hit.mother_particle().name();
    const dict::Key mothName = hit.sdata( "true_moth" );
    me->set_name(motherP, mothName);
  }
  
  return me; 
}

//*************************************************************
void fitter::Finalize() {
  //*************************************************************
   
  get_classifier().Finalize();
  Reset();

  if ( _doClust ) delete _clusters;
  
}

//*************************************************************
void fitter::ComputeSeed(const Trajectory& traj, State& seedState, int firsthit) {
  //*************************************************************

  _m.message("+++ computeSeed function ++++",bhep::VERBOSE);

  //use position slightly offset from first meas as seed 
  ///_lastIso is the total no of candidate muon hits inside Traj in free section 
  
  /*
  if ( (double)(traj.quality("lastIso"))/(double)traj.size() > _min_iso_prop )
    firsthit = (int)traj.size() - (int)(traj.quality("lastIso"));
  */

  EVector v(6,0), v2(1,0);
  EMatrix C(6,6,0), C2(1,1,0);
    
  v[0] = traj.nodes()[firsthit]->measurement().position()[0];
  v[1] = traj.nodes()[firsthit]->measurement().position()[1];
  v[2] = traj.nodes()[firsthit]->measurement().position()[2];  

  cout<<"z pos in fitter_seed: "<<v[2]<<endl; 

  // Estime the momentum from range
  ComputeMomFromRange( traj, (int)traj.size(), firsthit, v);

  double pSeed;
  //double wFe = _geom.get_Fe_prop();
  //Approximate p from plot of p vs. no. hits, then approx. de_dx from this.
  //if (v[5] == 0) { //pSeed = (double)(0.060*traj.nmeas())*bhep::GeV;
  //pSeed = (13300-11200*wFe) + (-128+190*wFe)*(double)traj.size();
  //v[5] = 1.0/pSeed;
  //}
  


  v[3] = 0;
  v[4] = 0;


  /*
  // Create and fill the state vector 
  EVector v(6,0); 
  v[0]= 0; // x position of the state 
  v[1]= 0; // y position 
  v[2]= 0; // z position 
  v[3]= 0; // x slope (dx/dz) 
  v[4]= 0; // y slope (dy/dz) 
  v[5]= 1; // q/p (charge over momentum) 
  
  // Create and fill the state matrix 
  EMatrix C(6,6,0); 
  C[0][0] = C[1][1] = 1; // square of the position error 
  C[2][2] = 0 ; // no error in z since this is the running coordinate 
  C[3][3] = C[4][4] = 0.1; // square of the slope error 
  C[5][5] = 0.1; // square of 1/p error 
  
  // Create and fill the State it self 
  State state; 
  // The State representation 
  state.set_name(RP::rep, RP::slopes_curv_z); 
  // The main HyperVector 
  state.set_hv(HyperVector(v,C)); 
  // The secondary sense HyperVector with sense=1 and no error 
  state.set_hv(RP::sense, HyperVector(1)); 
  */


  C[0][0] = 8.5 * 8.5 * cm * cm;
  C[1][1] = 1.5 * 1.5 *cm * cm;  // Expected pos res x,y
  C[2][2] = EGeo::zero_cov()/2;  // Expected pos res z
  //C[3][3] = C[4][4] = 1.; // Expected pos dx/dz, dy/dz
  C[3][3] = 8.5 * 8.5 * cm * cm;
  C[4][4] = 1.5 * 1.5 * cm *cm; // Expected pos dx/dz, dy/dz
  C[5][5] = pow(v[5],2); // Expected pos dx/dz, dy/dz


  v2[0] = 1;

  seedState.set_name(RP::representation, RP::slopes_curv_z); 
  // The main HyperVector 
  seedState.set_hv(HyperVector(v,C)); 
  // The secondary sense HyperVector with sense=1 and no error 
  //seedState.set_hv(RP::sense, HyperVector(1)); 

  seedState.set_hv(RP::sense,HyperVector(v2,C2));


  //std::cout<<"fitter::ComputeSeed 1/v[5]="<<1./v[5]<<std::endl;

  _m.message("++ Seed estate after setSeed() in fitter:",seedState,bhep::VERBOSE);
}

//*************************************************************
void fitter::ApplyCovarianceFactor(double factor, EMatrix& C0){
  //*************************************************************
    
  //--- a large diagonal covariance matrix ---//
  C0 *= factor;
}

//*****************************************************************************
double fitf2(Double_t *x,Double_t *par) { 
  //*****************************************************************************

  double z = x[0];
  double fitval = par[0] + par[1]*z+par[2]*z*z+par[3]*z*z*z+par[4]*z*z*z*z;

  return fitval;
}

//*****************************************************************************
void fitter::ComputeMomFromRange(const Trajectory& traj, int nplanes, int firsthit, EVector& V){
  //*****************************************************************************

  //Some catchers for pointless returns.
  int fitcatch;
  //
  /// int nfit;
  /// int fitRange[3];
  const int fitpoints = nplanes - firsthit;
  ///double meanchange = 0;
  double xpos[fitpoints], ypos[fitpoints], zpos[fitpoints];
  double upos[fitpoints];/// wpos[fitpoints];
  std::vector<EVector> dr;
  std::vector<EVector> B;
  bool isContained = true, cuspfound = false;

  double Xmax = _geom.getPlaneX() - 1*cm;
  double Ymax = _geom.getPlaneY() - 1*cm;
  /// double Zmax = _geom.getPlaneZ() - 1*cm;
  //double dx[fitpoints-1], dy[fitpoints-1], dz[fitpoints-1];
  // double ax[fitpoints-2], ay[fitpoints-2], az[fitpoints-2];
  // double bx[fitpoints-2], by[fitpoints-2], bz[fitpoints-2];
  
  ///double ds0=0, ds1=0;
  double Bmean=0;
  double pathlength=0;
  int Npts=0;
  ///double initR = 0;
  double sumDR = 0;
  int minindex = nplanes - firsthit;
  ///double minR = 999999.9999;
  double pdR = 0.0;
  
  EVector Z = EVector(3,0); Z[2] = 1;
  for (int ipoint=firsthit;ipoint < nplanes;ipoint++){
    
    xpos[ipoint-firsthit] = traj.node(ipoint).measurement().position()[0];
    ypos[ipoint-firsthit] = traj.node(ipoint).measurement().position()[1];
    zpos[ipoint-firsthit] = traj.node(ipoint).measurement().position()[2]
      - traj.node(firsthit).measurement().position()[2];
    if(fabs(xpos[ipoint-firsthit]) > Xmax || fabs(ypos[ipoint-firsthit]) > Ymax)
      isContained = false;
    else if(fabs(ypos[ipoint-firsthit]) > 
	    (1 + tan(atan(1.)/2.)) * Xmax - fabs(xpos[ipoint-firsthit])) 
      isContained = false;
    EVector pos0 = EVector(3,0);
    pos0[0] = xpos[ipoint-firsthit];
    pos0[1] = ypos[ipoint-firsthit];
    pos0[2] = zpos[ipoint-firsthit];
    EVector B0 = _geom.getBField(pos0);
    //std::cout<<"B0 in fitter::ComputeMomFromRange: "<<B0[0]<<" "<<B0[1]<<" "<<B0[2]<<" at z: "<<pos0[2]<<std::endl;
    B.push_back(B0);
    Bmean += B0.norm();
    upos[ipoint-firsthit] = // sqrt(pos0[0]*pos0[0] + pos0[1]*pos0[1]);
      dot(pos0,crossprod(Z, B0))/crossprod(Z, B0).norm();
    //if(!cuspfound)
    //  if(ipoint == firsthit) initR = upos[ipoint-firsthit];
    //  else {
    //sumDR += initR - upos[ipoint-firsthit];
    //initR = upos[ipoint - firsthit];
    //  }
    Npts++;
    if ( ipoint > firsthit){
      EVector drtemp = EVector(3,0);
      drtemp[0] = xpos[ipoint-firsthit] - xpos[ipoint-firsthit-1];
      drtemp[1] = ypos[ipoint-firsthit] - ypos[ipoint-firsthit-1];
      drtemp[2] = zpos[ipoint-firsthit] - zpos[ipoint-firsthit-1];
      dr.push_back(drtemp);      
      pathlength +=  drtemp.norm();
      if ( ipoint > firsthit + 1 ) {
	int k = ipoint-firsthit-1;
	EVector dr0 = dr[k-1];
	EVector dr1 = dr[k];
	EVector ddr = dr1 + dr0;
	EVector Ddr = dr1 - dr0;
	EVector pos = EVector(3,0);
	pos[0] = xpos[k-1]; pos[1] = ypos[k-1]; pos[2] = zpos[k-1]; 
	//	EVector B = _geom.getBField(pos);
	//std::cout<<"B in fitter::ComputeMomFromRange: "<<B[0]<<" "<<B[1]<<" "<<B[2]<<std::endl;
	double dR = dot(ddr, crossprod(Z, B0))/ (crossprod(Z,B0).norm());
	double DR = dot(Ddr, crossprod(Z, B0))/ (crossprod(Z,B0).norm());
	if(pdR != 0.0){
	  if(!cuspfound && DR/fabs(DR) == pdR/fabs(pdR)){
	    // sumDR += fabs(dR) > 0.0 ? dR/fabs(dR):0.0;
	    sumDR += dR;
	    // pdR = dR; 
	    pdR = dR;
	  }
	  else if(dR/fabs(dR) != pdR/fabs(pdR)){
	    // cuspfound = true;
	    minindex = ipoint - firsthit - 1;
	    pdR = dR;
	    // std::cout<<"At cusp, sumDR = "<<sumDR<<std::endl;
	  }
	}
	else if(!cuspfound && fabs(dR) > 0){
	  // sumDR += fabs(DR) > 0.0 ? DR/fabs(DR) : 0.0;
	  sumDR += dR;
	  pdR = dR;
	}

      }
    }
  }
  Bmean /=Npts;

  // Need to get the start Z-pos of the last Scintilator module.
  double zMax = 1600; //mm
  
  zMax = _geom.getZMax()-_geom.get_Fe_prop();
  //cout<<"zMax: "<<zMax<<endl;  
  //cout<< _geom.get_Fe_prop()<<endl;
  //cout<<traj.size()<<endl;


  zMax = 1000;


  //std::cout<<"pathLength: "<<pathlength<<std::endl;
  double final_Zpos=traj.node(nplanes -1).measurement().position()[2];
  //double final_Zpos=traj.nodes()[0]->measurement().position()[2];
  //std::cout<<"Final_Zpos"<<final_Zpos<<endl;

  double p;

  double meansign = 1;

  //meansign = CalculateCharge(traj);
  //p = RangeMomentum(pathlength,traj.node(firsthit).measurement().position()[2]);
  //p=fabs(MomentumFromCurvature(traj,0,p));//-p);


  _helix.clear();
  _quad.clear();
  _lever1.clear();
  _angle1.clear();
  _lever2.clear();
  _angle2.clear();
  vector<double> debug;

  p=MomentumFromCurvature(traj,0,p,debug);//-p);

  //cout<<"debug.size()="<<debug.size()<<endl;
  if(debug.size()==6)
    {
      _helix.push_back(debug[0]);
      _quad.push_back(debug[1]);
      _lever1.push_back(debug[2]);
      _angle1.push_back(debug[3]);
      _lever2.push_back(debug[4]);
      _angle2.push_back(debug[5]);
    }

  meansign = p/fabs(p);

  p = fabs(p);

  cout<<"Momentum used in ComputeMomFromRange in fitter="<<meansign * p<<endl;

  
  int pi3=0, pi4=0;
  while(dr.at(pi3)[2] == 0.0) pi3++;
  while(dr.at(pi4)[2] == 0.0) pi4++;
  V[3] = dr.at(pi3)[0]/dr.at(pi3)[2];
  V[4] = dr.at(pi4)[1]/dr.at(pi4)[2];

  V[5] = meansign/fabs(p);

  int sign = 1;
  if(meansign==meansign){
    if(meansign!=0)
      sign = int(meansign/fabs(meansign));
    else
      sign = 0;
  }
  else
    sign = 1;

  // std::cout<<"Pathlength is "<<pathlength // <<" or "<<pathlength0
  //	   <<" with charge "<<meansign<<std::endl;
  
  _initialqP = V[5];
  // _m.message("_initialqP ="<<_initialqP,bhep::VERBOSE);
  
}

//*****************************************************************************
void fitter::ReadParam(){
  //*****************************************************************************
    
  _m.message("+++ ReadParam function of fitter ++++",bhep::VERBOSE);
        
  _model = _store.fetch_sstore("model");//"particle/helix"; 
  dim=6; // ??????
    
  if ( _store.find_istore("refit") )
    _refit=_store.fetch_istore("refit");
  else _refit=false;

  if ( _store.find_istore("patRec") )
    _patternRec=_store.fetch_istore("patRec");
  else _patternRec=false;

  if ( _store.find_istore("do_clust") )
    _doClust = _store.fetch_istore("do_clust");
  else _doClust = false;

  _facRef = _store.fetch_dstore("facRef");

  _min_seed_hits = _store.fetch_istore("min_seed_hits");
  _min_iso_prop = _store.fetch_dstore("min_iso_prop");

  _chi2fit_max = _store.fetch_dstore("chi2fit_max");
    
  _X0 = _store.fetch_dstore("x0Fe") * mm;
  //_tolerance = _store.fetch_dstore("pos_res") * cm;
  _highPass = _store.fetch_istore("high_Pass_hits");
  _lowPass = _store.fetch_istore("low_Pass_hits");
  _lowFit1 = _store.fetch_dstore("low_fit_cut0");
  _lowFit2 = _store.fetch_dstore("low_fit_cut2");
      
}
