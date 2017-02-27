
/* -*- mode: c++ -*- */
#ifndef _fitter___
#define _fitter___

#include <recpack/RecPackManager.h>
#include <recpack/RayTool.h>
#include <recpack/KalmanFitter.h>
#include <recpack/HelixEquation.h>
#include <recpack/LsqFitter.h>
#include <recpack/ParticleState.h>

#include <mind/MINDsetup.h>
#include <mind/Utilities.h>
#include <mind/MINDfitman.h>
#include <mind/cluster.h>
#include <mind/event_classif.h>
#include <mind/hit_clusterer.h>

#include <mind/line.h>

//#include <mind/super_fit.h>

#include <bhep/event.h>
#include <bhep/gstore.h>

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>

using namespace Recpack;

class fitter : public super_fit{
  
public:
  
  fitter(const bhep::gstore& pstore,
	    bhep::prlevel vlevel=bhep::NORMAL);
    
  virtual ~fitter();
    
  //------------------ main functions ---------------// 
  //bool initialize(const bhep::sstore&) ;
  void Initialize();
  bool Execute(bhep::particle& part,int evNo);
  void Finalize() ;
  //-------------------------------------------------// 
  //Getters.    
  
  std::vector<Trajectory*>&  get_trajs(){return _trajs; }///
  Trajectory& get_hadTrajs() {return _hadTrajs;}

  std::vector<EVector>& get_had_unit(){ return _hadUnit; }
  
  double* get_had_eng(){ return _hadEng; }
  double* get_had_profile() { return _transEdx; }
  std::vector<double>& get_nonMuon_edep() { return _nonMuonEdep; }
  std::vector<int>& get_had_hits() { return _nonMuonHits; }
  int get_fail_event(){ return _failEvent; }///
  bool CheckReseed(){ return _reseed_ok; }
  std::vector<EVector>& get_had_RecDir(){ return _showerDir; }
  double showerVertZ(int i) { return _showerVertZ[i]; }
  double showerNplanes(int i) { return _showerNplanes[i]; }
  double showerXtent(int i) { return _showerXtent[i]; }

  cluster* GetMeas(int num){return _meas[num];}
  cluster* GetMeasTASD(int num){return _measTASD[num];}

  std::vector<cluster*> &GetMeasVec(){ return _meas; }
  int GetNMeas(){ return (int)_meas.size();}
  int GetNMeasTASD(){ return (int)_measTASD.size();}
  
  MINDsetup GetGeom() { return _geom; }
  //
  cluster* GetMeasurement(bhep::hit& hit);
  
  int GetQ(const Trajectory& traj);

  double GetInitialqP(){ return _initialqP; }

  int* getMuonIndex(){ return _muonindex; }
  //Tempory for likelihoods.
  void set_int_type(const string name){
    get_classifier().set_int_type( name ); }

  //recpack manager
  
  RecPackManager& man(){
    return MINDfitman::instance().manager();}

  event_classif& get_classifier(){ return _classify; }

  double get_momentum_guess(int pos){return _momentum_guess_vec[pos];}
  double get_momentum_guess_length(){return _momentum_guess_vec.size();}

  vector<double> _momentum_guess_vec;

  std::vector<double> _xDir;
  std::vector<double> _yDir;
  std::vector<double> _x0;
  std::vector<double> _y0;
  std::vector<double> _xchi;
  std::vector<double> _ychi;
  std::vector<double> _hitsPerPlanes;
  std::vector<double> _avrHitsPerUsedPlanes;

  std::vector<double> _helix;
  std::vector<double> _quad;
  std::vector<double> _lever1;
  std::vector<double> _angle1;
  std::vector<double> _lever2;
  std::vector<double> _angle2;
  
  std::vector<double> GetHelix() {return _helix;}
  std::vector<double> GetQuad() {return _quad;}
  std::vector<double> GetLever1() {return _lever1;}
  std::vector<double> GetAngle1() {return _angle1;}
  std::vector<double> GetLever2() {return _lever2;}
  std::vector<double> GetAngle2() {return _angle2;}


  std::vector<double> GetXDir() {return _xDir;}
  std::vector<double> GetYDir() {return _yDir;}
  std::vector<double> GetX0() {return _x0;}
  std::vector<double> GetY0() {return _y0;}
  std::vector<double> GetXChi() {return _xchi;}
  std::vector<double> GetYChi() {return _ychi;}
  std::vector<double> GetHPP() {return _hitsPerPlanes;}
  std::vector<double> GetAHPP() {return _avrHitsPerUsedPlanes;}


protected:
  
  //void resetVirtualPlanes(); 

  //read parameters from store
  void ReadParam();

  void TASDtracker();

  void TASDtracker2();

  double CreatePlanesWithHits(const std::vector<cluster*>&  meas, std::vector<plane_info*>& planes);

  void AddHitsToTrack(const std::vector<plane_info*>& planes, Trajectory* traj);

  void AddHitsToTrack2(const std::vector<plane_info*>& planes, Trajectory* traj);


  //void CreateTrackFromHits(const std::vector<plane_info*>& planes, Trajectory* traj);
    
  //seed for fit
  void ComputeSeed(const Trajectory& traj,State& seed, int firsthit=0);
  //void ComputeMomFromParabola(const Trajectory& traj, int nplanes, int firsthit, EVector& V);
  void ComputeMomFromRange(const Trajectory& traj, int nplanes, int firsthit, EVector& V);

  //double RangeMomentum(double length,double nodeZ);
  //double MomentumFromCurvature(const Trajectory& traj, int firsthit = 0);
  double MomentumFromDeflection(const Trajectory& traj, int firsthit = 0);//not used! For low pt.

  //seed error
  void ApplyCovarianceFactor(double factor, EMatrix& C0);
 
  //fit trajectory
  bool FitTrajectory(const State& seed, const int trajno);
  //bool ReseedTrajectory(Trajectory& traj,const int trajno);
  bool ReseedTrajectory(const int trajno);
  //bool fitHadrons();
  //double eng_scale(double visEng);

  // hadron
  void rec_had_energy();
  void rec_had_edep(int j);

  //-------- get traj from event----------//
  bool readTrajectory(const bhep::particle& part);

  // Create a single trajectory with all measurements when the PR is off
  bool CreateSingleTrajectory(Trajectory& traj); 

  // Create measurements from hits
  bool CreateMeasurements(const bhep::particle& part); 


  // Check traj passes cuts for fitting.
  bool CheckValidTraj(const Trajectory& traj);


  ///calcute seed for refit
  void ComputeSeedRefit(const Trajectory& traj, State& seedState);
  //string getPlaneName(bhep::hit);
  //--------------------------------------//
  
  bool CheckQuality(const Trajectory& traj);
    
  void Reset();

protected:

  bhep::gstore _store;
  
  bhep::prlevel _level;
    
  bhep::messenger _m;
    
  MINDsetup _geom;
  
  //counter for virtual planes
  //size_t pnumber;
  
  //Parameters to define fitting method.
  bool _refit; //Do second fit.
  bool _patternRec; //Pattern recognition algorithm required?
  
  int _fitCheck;

  int _forwardFitCheck;
  int _reseedFitCheck;

  int _min_seed_hits; //Minimum isolated hits required for Prec seed.
  double _min_iso_prop;

  //bit to tell if reseed perfromed
  bool _reseed_ok;
  bool _reseed_called ;
  bool _fitted;

  


  //------------------ Physics -----------------//
    
  double _X0;
  double _rel_dedx;
  double _rel_dens;
  double _widthI;
  double _widthS;

  int dim; //dimension of model state
  //int meas_dim; //dimension of measurement
  //double _tolerance; //pos. resolution/plane tolerance
  
  State _seedstate;   
  //EVector qoverp;
  double _initialqP;
  
  string _model;  // fit model
  //string kfitter; // kind of fit
  
  //double chi2node_max;
  //int max_outliers;
  double _chi2fit_max;
  double _facRef;

  //Detector name for hit getter. made members 26/11, weird error day.
  string _detect;
  int _testBeam;
  int _highPass;
  int _lowPass;
  double _lowFit1, _lowFit2;
  int _muonindex[2];
  std::vector<double> _nonMuonEdep ;
  std::vector<int> _nonMuonHits ;
  std::vector<EVector> _showerDir;
  std::vector<double> _showerVertZ;
  std::vector<double> _showerNplanes;
  std::vector<double> _showerXtent;

  Trajectory _traj;
  Trajectory _traj2;
  Trajectory _traj3;
  std::vector<Trajectory*> _trajs;  //
  Trajectory _hadTrajs;
  

  std::vector<cluster*> _measTASD;
  std::vector<cluster*> _meas;
  std::vector<cluster*> _hadmeas;

 

  //value set to identify where a traj failed:
  //0=too few hits. 1=too many hits. 2=outside fiducial. 3=no convergence with kink.
  //4=Couldn't find seed for patrec. 5=Failed in pat rec filtering
  int _failType;
  int _failEvent;
  int _intType ;///
  int _pr_count;///
  // int reseed_count;

  //size_t nnodes;

  //Temporary fix (??). vector to store hadron unit dir vec.
  std::vector<EVector> _hadUnit;
  double _hadEng[2];
  double _transEdx[2];
  double _length;
 
  // Stuff relevant for event classification, will be uncommented when needed.
  event_classif _classify;

  // hit clustering.
  hit_clusterer *_clusters;

  bool _doClust;
  //
  std::vector<State> _vPR_seed;
  //-------------- verbosity levels ------------//

  //int vfit,vnav,vmod,vmat,vsim;
  
};

class reverseSorter{
public:
  bool operator()(const cluster* p1, const cluster* p2){
    if (p2->position()[2] < p1->position()[2]) return true;
    return false;
  }

};

class forwardSorter{
public:
  bool operator()(const cluster* p1, const cluster* p2){
    if (p2->position()[2] > p1->position()[2]) return true;
    return false;
  }

};
  
#endif

