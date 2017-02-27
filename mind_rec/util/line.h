#ifndef _line___
#define _line___

#include <recpack/RecPackManager.h>
//#include <mind/SetupSk.h>
//#include <mind/MINDfieldMapReader.h>
//#include <mind/DeDxMap.h>

#include <TH1F.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TMath.h>
#include <mind/cluster.h>

#include <math.h>

#include <vector>

using namespace Recpack;

class Line{

public:
    
  Line(cluster* downStream, cluster* upStream);
  Line(Node* downStream, Node* upStream);
  virtual ~Line();

  double CalculateR(cluster* hit);
  //double CalculateR(Node* hit);

  //void AddHits(cluster* hit){hits.push_back(hit);};
  void AddHits(Node hit){hits.push_back(hit);};

  //std::vector<cluster*> GetHits(){return hits;};
  std::vector<Node> GetHits(){return hits;};

  void Equation();

 private:

  double _kx;
  double _ky;
  double _kz;
  double _mx;
  double _my;
  double _mz;
  
  //std::vector<cluster*> hits;
  std::vector<Node> hits;
};

#endif 
 
