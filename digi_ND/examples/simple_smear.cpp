/******************************************************************
 * Program to create a digital particle for fitting from the hits *
 * in the true particles in the bhep output of GEANT4 MIND.       *
 *                                                                *
 * Temporary measure until full digitization is available.        *
 *                                                                *
 * Execute as:                                                    *
 *       ./simple_smear <param_file_name>                         *
 *                                                                *
 * Author: Andrew Laing.                                          *
 ******************************************************************/

#include <digi/root2dst.h>
#include <digi/WriteUtil.h>

//#include <bhep/EventManager2.h>
#include <bhep/gstore.h>
#include <bhep/sreader.h>
#include <bhep/reader_root.h>
#include <bhep/bhep_svc.h>

#include <stdio.h>

#include <fstream>
#include <sstream>
#include <string>
#include <bitset> 
#include <iostream>
#include <stdlib.h> 

#include <map>

#include <digi/xml_parser.h>
#include <sqlite3.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"

#include <exception>

#include "digi/MDfragmentBM.h"
#include "digi/MDpartEventBM.h"
#include "digi/MDargumentHandler.h"
#include "digi/MDdataFile.h"

using namespace std;

// Clean up file and put unpacking data in separate file!

vector<string> dataBaseVector;

int myatoi(string text)
{
 int result;          //number which will contain the result
 istringstream convert(text); // stringstream used for the conversion constructed with the contents of 'Text' 
                             // ie: the stream will start containing the characters of 'Text'
 if ( !(convert >> result) ) result = 0;  //give the value to 'Result' using the characters in the stream
    
 return result;
}

int myatoi(char text)
{
 int result = (int)text -48;//All ascii numberes start at 48;
 return result;
}


class forwardSortInZ{
 public:
  bool operator()(const bhep::hit* p1, const bhep::hit* p2){
    if (p2->x()[2] > p1->x()[2]) return true;
    return false;
  }

};

class TempHit
{
public:
  TempHit() {};
  ~TempHit() {};
  TempHit(unsigned int amplitudeIn, unsigned long int timeingIn,
  unsigned int channelIdIn,unsigned int boardIdIn){
    amplitude=amplitudeIn;
    timeing=timeingIn;
    channelId=channelIdIn;
    boardId=boardIdIn;
  };
  unsigned int amplitude;
  unsigned long int timeing; //Time relative to spill.
  unsigned int channelId;
  unsigned int boardId;
};

double ConvertToEdep(unsigned int boardId, unsigned int channelId,unsigned int amplitude)
{
  return 1.0;
}

static int Callback(void *NotUsed, int argc, char **argv, char **azColName){
  int i;
  for(i=0; i<argc; i++){
    //printf("%s = %s\n", azColName[i], argv[i] ? argv[i] : "NULL");
    dataBaseVector.push_back(argv[i]);
   }
  //printf("\n");
  return 0;
}

void ConvertToPosition(bhep::hit* inHit, unsigned int boardId, unsigned int channelId)
{
  // In convert to position2017
  double xBarStart = -970.55;
  double yBarStart = -1139.7;
  double xBarPitch = 41.3/2;
  double yBarPitch =(438.3-87.6)/2/2; //Div by two to make algorithm work.
  double detectorZ = (1635 +3 +4)/2;
  double FeDepth = 30;
  double SciDepth = 4*7.5 + 2*5;
  double Gap = 460;
  double GapFeSci = 8;
  double GapFeFe = 30;

  double usingTASD = 1;
  double addingSize = 3850;

  double mod1ZPos = -2823 + addingSize*usingTASD; //SFFFS0
  double mod2ZPos = -2585 + addingSize*usingTASD; //SFFFS0
  double mod3ZPos = -2126 + addingSize*usingTASD; //SFFFS1
  double mod4ZPos = -1895 + addingSize*usingTASD; //SFFFS1
  double mod5ZPos = -1435 + addingSize*usingTASD; //SFModule
  double mod6ZPos = -1322 + addingSize*usingTASD; //SFModule
  double mod7ZPos = -1208 + addingSize*usingTASD; //SFModule
  double mod8ZPos = -109 + addingSize*usingTASD; //SFFModule
  double mod9ZPos = 67 + addingSize*usingTASD; // SFFModule
  double mod10ZPos = 242 + addingSize*usingTASD; //SFFModule
  double mod11ZPos = 423 + addingSize*usingTASD; //SFFModule
  double mod12ZPos = 1061 + addingSize*usingTASD; //SFFFF 
  double mod13ZPos = 1340 + addingSize*usingTASD; //SFFFFS0
  double mod14ZPos = 1624 + addingSize*usingTASD; //SFFFFS0
  double mod15ZPos = 2089 + addingSize*usingTASD; //SFFFFS1
  double mod16ZPos = 2394 + addingSize*usingTASD; //SFFFFS1
  double mod17ZPos = 2855 + addingSize*usingTASD; //SFFFFS2
  double mod18ZPos = 3158 + addingSize*usingTASD; //SFFFFS2

  double moduleList[18] = {mod1ZPos, mod2ZPos, mod3ZPos, mod4ZPos, mod5ZPos, mod6ZPos, mod7ZPos, mod8ZPos, mod9ZPos, mod10ZPos, mod11ZPos, mod12ZPos, mod13ZPos, mod14ZPos, mod15ZPos, mod16ZPos, mod17ZPos, mod18ZPos};

  /*
  double mod1ZPos = -2823; //-2700.5+SciDepth/2;//-detectorZ + SciDepth/2;
  double mod2ZPos = mod1ZPos + SciDepth/2 + 2*FeDepth + 2*GapFeFe +2*GapFeSci + SciDepth/2;
  double mod3ZPos = mod2ZPos + SciDepth/2 + Gap + SciDepth/2;
  double mod4ZPos = mod3ZPos + SciDepth/2 + 2*FeDepth + 2*GapFeFe +2*GapFeSci + SciDepth/2;
  double mod5ZPos = mod4ZPos + SciDepth/2 + Gap + SciDepth/2;
  double mod6ZPos = mod5ZPos + SciDepth/2 + FeDepth + 2*GapFeSci + SciDepth/2;
  double mod7ZPos = mod6ZPos + SciDepth/2 + FeDepth + 2*GapFeSci + SciDepth/2; // detectorZ-FeDepth- GapFeSci - SciDepth/2;
  double moduleList[7] = {mod1ZPos, mod2ZPos, mod3ZPos, mod4ZPos, mod5ZPos, mod6ZPos, mod7ZPos};
  */
  double xBarZsize = 0;//7.5 *2; //xYYx
  double yBarZsize = 0;//7.5;
  
  //bool filled = false;

  for(int i=0;i<dataBaseVector.size();i+=6)
    {
      int Global_channel = myatoi(dataBaseVector[i]);
      //if(Global_channel!=channelId)

      //if(i+6>dataBaseVector.size())
      //cout<<"broken"<<Global_channel<<"\t"<<channelId<<endl;
      if(Global_channel==channelId)
	{
	  //filled = true;
	  //cerr<<"through if check"<<endl;
	  string Bar_type = dataBaseVector[i+1];
	  int Module = myatoi(dataBaseVector[i+2]);
	  string Position_label = dataBaseVector[i+3];
	  string Bar_position = dataBaseVector[i+4];
	  string Side = dataBaseVector[i+5];

	  //cerr<<boardId<<"\t"<<channelId<<"\t"<<Bar_type<<endl;
	  
	  double barPosZ =0.0;
	  int YMeasuring=0;
	  int TASD = 0;
	  int barNumber=0;
	  double barPosT=0.0;
	  //int moduleNum=1;
	  
	  barPosZ = mod1ZPos;
	  /*
	  cout<<"Bar_type="<<Bar_type<<"\t"
	      <<"Bar_position="<<Bar_position<<"\t"
	      <<"moduleNum="<<moduleNum<<"\t"
	      <<"Position_label="<<Position_label[1]<<endl;
	  */

	  double position;
	  if(Position_label.length() == 2)
	    position = myatoi(Position_label[1]);
	  else
	    position = myatoi(Position_label[1])*10+ myatoi(Position_label[2]);
	    
	  if(position == 0 || position == 32)
	    {
	      TASD = 1;
	      if(position == 0)
		{
		  barPosZ = -2663;
		  if(Bar_type=="Vertical")
		    {
		      int temp = Global_channel % 96;
		      barPosT = -((temp-32)*10 -31*10);
		      barNumber = temp-32;
		      YMeasuring=0;
		    }
		  else
		    {
		      int temp = Global_channel % 96;
		      if(temp == 32)
			{
			  barPosT = -5;
			  barNumber = 63-33;
			}
		      else
			{
			  barPosT = (temp-33)*10 -30*10 +117; //117 since aida higher.
			  barNumber = temp-33;
			}
		      YMeasuring=1;
		    }
		}
	      else //(position == 32)
		{
		  barPosZ = -2663 - 373;
		  if(Bar_type=="Vertical")
		    {
		      int temp = Global_channel % 96;
		      barPosT = -((temp-0)*10 -7*10);
		      barNumber = temp;
		      YMeasuring=0;
		    }
		  else
		    {
		      int temp = Global_channel % 96;
	
		      barPosT = (temp-8)*10 -7*10 +117; //117 since aida higher.
		      barNumber = temp-8;
		      YMeasuring=1;
		    }
		}
	    }
	  else
	    {
	      TASD =0;
	      for (int moduleNum = 1;moduleNum<=18;moduleNum++)
		{
		  if(Bar_type=="Vertical" && position == moduleNum && Bar_position =="Front") barPosZ = moduleList[moduleNum-1]-xBarZsize;
		  else if(Bar_type=="Vertical" && position == moduleNum && Bar_position =="Back") barPosZ = moduleList[moduleNum-1]+xBarZsize;
		  else if(Bar_type=="Horizontal" && position == moduleNum && Bar_position =="Front") barPosZ = moduleList[moduleNum-1]-yBarZsize;
		  else if(Bar_type=="Horizontal" && position == moduleNum && Bar_position =="Back") barPosZ = moduleList[moduleNum-1]+yBarZsize;
		}
	    
	      //cout<<"barPosZ="<<barPosZ<<endl;
	      /*
		cerr<<"Position_label="<<Position_label<<"\t"
		<<"position="<<position<<"\t"
		<<"Global_channel="<<Global_channel<<"\t"
		<<"barPosZ="<<barPosZ<<endl;
	      */
	      if(Bar_type=="Vertical")
		{
		  //int temp = Global_channel % 96;
		  
		  if(Side =="Left")
		    {
		      int temp = Global_channel % 96;
		      barPosT = (temp % 32 )*yBarPitch +yBarStart;
		      barNumber = temp % 32;
		      //barPosT = (Global_channel -1568 +32 * myatoi(Position_label[1]))*yBarPitch +yBarStart;
		      //barNumber = (Global_channel -1568 +32 - 32 * myatoi(Position_label[1]));
		    }	
		  else
		    {
		      int temp = (Global_channel-1) % 96;
		      barPosT = (temp % 32 )*yBarPitch +yBarStart;
		      barNumber = temp % 32;
		      
		      //barPosT = (Global_channel -1568-1 +32 - 32 * myatoi(Position_label[1]))*yBarPitch +yBarStart;
		      //barNumber = (Global_channel -1568-1 +32 - 32 * myatoi(Position_label[1]));
		    }
		  YMeasuring=0;
		}
	      else
		{
		  int temp = Global_channel % 96;
		  barPosT = temp*xBarPitch +xBarStart;
		  barNumber = temp;
		  /*
		    if(Side =="Left")
		    { 
		    barPosT = (Global_channel - 96 * myatoi(Position_label[1]))*xBarPitch +xBarStart;
		    barNumber = (Global_channel - 96 * myatoi(Position_label[1]));
		    }
		    else
		    {
		    barPosT = (Global_channel - 96 * myatoi(Position_label[1]))*xBarPitch +xBarStart;
		    barNumber =  (Global_channel - 96 * myatoi(Position_label[1]));
		    }
		  */
		  YMeasuring=1;
		}
	    }
	  //cout<<"barPosZ="<<barPosZ<<"\t"<<"barPosT="<<barPosT<<"\t"<<"moduleNum="<<position<<"\t"<<"YMeasuring="<<YMeasuring<<endl;
	  
	  //cerr<<"Adding properties barPosZ"<<endl;
	  inHit->add_property("barPosZ",barPosZ);
	  //inHit->add_property("barPosZ",mod1ZPos);
	  inHit->add_property("IsYBar",YMeasuring);
	  inHit->add_property("barNumber",barNumber);
	  inHit->add_property("barPosT",barPosT);
	  inHit->add_property("IsTASD",TASD);
	  inHit->add_property("moduleNum",position);
	  //inHit->add_property("moduleNum",myatoi(Position_label[1]));
	  // can we add ->x()[2] how is it done? bhep add?
	  if(YMeasuring) inHit->set_point(*(new bhep::Point3D(0,barPosT,barPosZ)));
	  else inHit->set_point(*(new bhep::Point3D(barPosT,0,barPosZ)));
	  break;
	}

      //if(!filled && i+6>dataBaseVector.size())
      //cout<<"broken"<<Global_channel<<"\t"<<channelId<<endl;

      /*
      else
	{
	  //cerr<<"else="<<Global_channel<<"\t"<<channelId<<endl;
	   // Only here to work around a bug.
	  inHit->add_property("barPosZ",0);
	  inHit->add_property("IsYBar",0);
	  inHit->add_property("IsTASD",1);
	  inHit->add_property("barNumber",0);
	  inHit->add_property("barPosT",0);
	  inHit->add_property("moduleNum",0);
	  inHit->set_point(*(new bhep::Point3D(0,0,0)));
      
	}
      */
      // Leaving convert to position2017
    }
  if(!inHit->find_idata_name("IsTASD"))
    {
      cerr<<"Broken"<<endl;
      cerr<<"else="<<boardId<<"\t"<<channelId<<"\t"<<channelId-96*boardId<<endl;
      inHit->add_property("barPosZ",0);
      inHit->add_property("IsYBar",0);
      inHit->add_property("IsTASD",1);
      inHit->add_property("barNumber",0);
      inHit->add_property("barPosT",0);
      inHit->add_property("moduleNum",0);
      inHit->set_point(*(new bhep::Point3D(0,0,0)));
    }

  //cout<<inHit->idata("IsTASD")<<endl;
  //cout<<inHit->ddata("moduleNum")<<endl;
}

vector<bhep::hit*> ConvertToHit(vector<TempHit*> hitsVector, string dataBaseName)
{
 /*
    Using information from the sqlite3 database we convert from the vector of TempHit to
    a vector of vectors of bhep hit*. The vectors are split depending on the different gtrigs.
 */
  vector<bhep::hit*>returnVector;
  sqlite3 *database;
  char *sql;
  char *zErrMsg = 0;
  int rc = sqlite3_open(dataBaseName.c_str(), &database);
  if( rc ){
    fprintf(stderr, "Can't open database: %s\n", sqlite3_errmsg(database));
    exit(0);
  }
  else{
    fprintf(stdout, "Opened database successfully\n");
  }
  
  sql = "SELECT Global_channel, Bar_type, Module, Position_label, Bar_position, Side from Scintillator_module";
 
  rc = sqlite3_exec(database, sql, Callback, 0, &zErrMsg);
  if( rc != SQLITE_OK )
    {
      fprintf(stderr, "SQL error: %s\n", zErrMsg);
      sqlite3_free(zErrMsg);
    }
  sqlite3_close(database);
  // Will fill dataBaseVector

 for(unsigned int cnt = 0; cnt<hitsVector.size(); cnt++)
  {
    //cerr<<"boardId="<<hitsVector[cnt]->boardId<<"\t"<<"channelId="<<hitsVector[cnt]->channelId<<endl;
    //cout<<"boardId="<<hitsVector[cnt]->boardId<<"\t"<<"channelId="<<hitsVector[cnt]->channelId<<endl;
    //if(hitsVector[cnt]->boardId == 0 || hitsVector[cnt]->boardId > 60) continue;
    //if(hitsVector[cnt]->boardId > 60) continue;

    returnVector.push_back(new bhep::hit());
    double tempEdep = ConvertToEdep(hitsVector[cnt]->boardId,hitsVector[cnt]->channelId,hitsVector[cnt]->amplitude);
    returnVector.back()->add_property("EnergyDep",tempEdep);
    returnVector.back()->add_property("time",(double)hitsVector[cnt]->timeing);
    
    // Right side hits, why?
    //if(hitsVector[cnt]->boardId==3 || (hitsVector[cnt]->boardId == 5 && hitsVector[cnt]->channelId % 2))
    //{
    //	returnVector.back()->add_property("rHit",1);
    //}
    //else
    //{
    returnVector.back()->add_property("rHit",0);
	//}
    //cerr<<"ConvertToPosition2017BeamTest2 start"<<endl;
    ConvertToPosition(returnVector.back(),hitsVector[cnt]->boardId,hitsVector[cnt]->channelId);
    //cerr<<"ConvertToPosition2017BeamTest2 end"<<endl;
  }
 hitsVector.clear();
 dataBaseVector.clear();

 return returnVector;
}

bool NullCheck(vector<char*> vector)
{
  bool ret = false;
  
  for(unsigned int i=0;i<vector.size();i++)
    {
      if(vector[i] == NULL)
	{
	  ret = true;
	  break;
	}
    }
  return ret;
}

void convert_to_internal(vector<char*> eventBuffers,vector<vector<TempHit*> >* vectorFinalHits)
{ 
  vector<pair<double,vector<TempHit*> > > tempHitsPerSpill[48]; 
  vector<double> TriggerNumWithHits;
  try {
 
    for(unsigned int eventNum=0;eventNum<eventBuffers.size();eventNum++)
      {
	vector<TempHit*> tempVec;
	MDfragmentBM   spill;
	spill.SetDataPtr(eventBuffers[eventNum]);
	int nTr = spill.GetNumOfTriggers();
	
	//if(spill.GetBoardId()== 3) cerr<<"before, read and 3"<<endl;
	
	cerr<<"In board="<<spill.GetBoardId()<<endl;
	//int nTr = 1000;
	//for (unsigned int triggNum=0; triggNum<TriggerNumWithHits.size(); triggNum++) {
	for (unsigned int triggNum=0; triggNum<nTr; triggNum++) {
	//for (unsigned int triggNum=2004; triggNum<2005; triggNum++) {
	  //MDpartEventBM * event = spill.GetTriggerEventPtr(TriggerNumWithHits[triggNum]);
	  MDpartEventBM * event = spill.GetTriggerEventPtr(triggNum);
	  //if(event->getNumDataWords() < 12) continue;

	  int numberOfHits=0;
	  for (int ich=0; ich<BM_FEB_NCHANNELS; ++ich)
	    {
	      int nHits = event->GetNLeadingEdgeHits(ich);
	      numberOfHits+=nHits;
	    }
	  
	  if(numberOfHits<1) continue;
	  
	  /*
	  cout<<"In board="<<spill.GetBoardId()<<"\t"
	      <<"triggNum="<<triggNum<<"\t"
	      <<"triggerTime="<<event->GetTriggerTime()<<"\t"
	      <<"triggerTag="<<event->GetTriggerTag()<<"\t"
	      <<"numberOfHits="<<numberOfHits<<endl;
	  */
	  //if(event->getNumDataWords() < 12) continue;
	  for (int ich=0; ich<BM_FEB_NCHANNELS; ++ich) {
	    
	    //Applitude filter
	    //if(event->GetHitAmplitude(ich, 'h')<1000) continue;
	    
	    int nHits = event->GetNLeadingEdgeHits(ich);
	    int otherSide = false;
	    //otherSide= event->GetNLeadingEdgeHits(ich+BM_FEB_NCHANNELS);
	    //if(nHits && otherSide)
	    if(nHits)
	      {
		//if(spill.GetBoardId()== 3) cerr<<"inside and 3"<<endl;
		for(int ih=0; ih<nHits; ih++)
		  {
		    if(event->GetHitTime(ih,ich, 'l')>0 && event->GetHitAmplitude(ich, 'h') >0)
		      //if(true)
		      {
			//if(spill.GetBoardId()== 3) cerr<<"inside, time and 3"<<endl;
			
			int GlobalChannel=ich;
			if(spill.GetBoardId()== 0) GlobalChannel+=0;
			else if(spill.GetBoardId()== 1) GlobalChannel+=96;
			else if(spill.GetBoardId()== 2) GlobalChannel+=2*96;
			else if(spill.GetBoardId()== 3) GlobalChannel+=3*96;
			else if(spill.GetBoardId()== 4) GlobalChannel+=4*96;
			else if(spill.GetBoardId()== 5) GlobalChannel+=5*96;
			else if(spill.GetBoardId()== 8) GlobalChannel+=8*96;
			else if(spill.GetBoardId()== 9) GlobalChannel+=9*96;
			else if(spill.GetBoardId()== 10) GlobalChannel+=10*96;
			else if(spill.GetBoardId()== 11) GlobalChannel+=11*96;
			else if(spill.GetBoardId()== 12) GlobalChannel+=12*96;
			else if(spill.GetBoardId()== 13) GlobalChannel+=13*96;
			else if(spill.GetBoardId()== 16) GlobalChannel+=16*96;
			else if(spill.GetBoardId()== 17) GlobalChannel+=17*96;
			else if(spill.GetBoardId()== 18) GlobalChannel+=18*96;
			else if(spill.GetBoardId()== 19) GlobalChannel+=19*96;
			else if(spill.GetBoardId()== 20) GlobalChannel+=20*96;
			else if(spill.GetBoardId()== 24) GlobalChannel+=24*96;
			else if(spill.GetBoardId()== 25) GlobalChannel+=25*96;
			else if(spill.GetBoardId()== 26) GlobalChannel+=26*96;
			else if(spill.GetBoardId()== 27) GlobalChannel+=27*96;
			else if(spill.GetBoardId()== 28) GlobalChannel+=28*96;
			else if(spill.GetBoardId()== 32) GlobalChannel+=32*96;
			else if(spill.GetBoardId()== 33) GlobalChannel+=33*96;
			else if(spill.GetBoardId()== 34) GlobalChannel+=34*96;
			else if(spill.GetBoardId()== 35) GlobalChannel+=35*96;
			else if(spill.GetBoardId()== 36) GlobalChannel+=36*96;
			else if(spill.GetBoardId()== 37) GlobalChannel+=37*96;
			else if(spill.GetBoardId()== 40) GlobalChannel+=40*96;
			else if(spill.GetBoardId()== 41) GlobalChannel+=41*96;
			else if(spill.GetBoardId()== 42) GlobalChannel+=42*96;
			else if(spill.GetBoardId()== 43) GlobalChannel+=43*96;
			else if(spill.GetBoardId()== 44) GlobalChannel+=44*96;
			else if(spill.GetBoardId()== 48) GlobalChannel+=48*96;
			else if(spill.GetBoardId()== 49) GlobalChannel+=49*96;
			else if(spill.GetBoardId()== 50) GlobalChannel+=50*96;
			else if(spill.GetBoardId()== 51) GlobalChannel+=51*96;
			else if(spill.GetBoardId()== 52) GlobalChannel+=52*96;
			else if(spill.GetBoardId()== 53) GlobalChannel+=53*96;
			else if(spill.GetBoardId()== 56) GlobalChannel+=56*96;
			else if(spill.GetBoardId()== 57) GlobalChannel+=57*96;
			else if(spill.GetBoardId()== 58) GlobalChannel+=58*96;
			else if(spill.GetBoardId()== 59) GlobalChannel+=59*96;
			else if(spill.GetBoardId()== 60) GlobalChannel+=60*96;
			else //(spill.GetBoardId()== 61) 
			  GlobalChannel+=61*96;

			/*
			//if(spill.GetBoardId()== 28)
			  //{
			    //    cerr<<"28"<<endl;
			    cout<<"spill.GetBoardId()="<<spill.GetBoardId()<<"\t"
				<<"ich="<<ich<<"\t"<<"GlobalChannel="<<GlobalChannel<<endl;
			    cout<<event->GetHitAmplitude(ich, 'h')<<"\t"
				<<event->GetHitTime(ih,ich, 'l')<<"\t"
				<<event->GetTriggerTime()<<endl;
			    //}
			    */
			if(GlobalChannel==4991) continue; // strange bug
			else if(GlobalChannel>783 && GlobalChannel<800) continue; // board 8
			else if(GlobalChannel>23 && GlobalChannel<32) continue; // board 0

			// Bug? Horizontal should not have channel 95
			if((spill.GetBoardId()== 1 ||
			    spill.GetBoardId()== 2 ||
			    spill.GetBoardId()== 3 ||
			    spill.GetBoardId()== 4 ||
			    spill.GetBoardId()== 5 ||
			    spill.GetBoardId()== 9 ||
			    spill.GetBoardId()== 10 ||
			    spill.GetBoardId()== 11 ||
			    spill.GetBoardId()== 12 ||
			    spill.GetBoardId()== 13 ||
			    spill.GetBoardId()== 16 ||
			    spill.GetBoardId()== 17 ||
			    spill.GetBoardId()== 18 ||
			    spill.GetBoardId()== 24 ||
			    spill.GetBoardId()== 25 ||
			    spill.GetBoardId()== 26 ||
			    spill.GetBoardId()== 32 ||
			    spill.GetBoardId()== 33 ||
			    spill.GetBoardId()== 34 ||
			    spill.GetBoardId()== 35 ||
			    spill.GetBoardId()== 40 ||
			    spill.GetBoardId()== 41 ||
			    spill.GetBoardId()== 42 ||
			    spill.GetBoardId()== 43 ||
			    spill.GetBoardId()== 48 ||
			    spill.GetBoardId()== 49 ||
			    spill.GetBoardId()== 50 ||
			    spill.GetBoardId()== 51 ||
			    spill.GetBoardId()== 52 ||
			    spill.GetBoardId()== 53 ||
			    spill.GetBoardId()== 56 ||
			    spill.GetBoardId()== 57 ||
			    spill.GetBoardId()== 58 ||
			    spill.GetBoardId()== 59 ||
			    spill.GetBoardId()== 60 ||
			    spill.GetBoardId()== 61) && (GlobalChannel % 96 == 95)) continue;


			tempVec.push_back(new TempHit(event->GetHitAmplitude(ich, 'h'),
						      event->GetHitTime(ih,ich, 'l'),
						      GlobalChannel,spill.GetBoardId()));
		      }
		  }
	      }
	  }//END for channels
	  if(tempVec.size())
	    //tempHitsPerSpill[eventNum].push_back(make_pair(event->GetTriggerTime(),tempVec));
	    tempHitsPerSpill[eventNum].push_back(make_pair(event->GetTriggerTag(),tempVec));
	  tempVec.clear();
	  }//END of trigger
	  spill.Clean();
	}
	//END of eventNum;
      } catch (MDexception & lExc)  {
      std::cerr <<  lExc.GetDescription() << endl
		<< "Unpacking exception\n"
		<< "Spill skipped!\n\n";
      for(unsigned int i = 0; i<48; i++)// for tempHitsPerSpill.size()
	{ 
      tempHitsPerSpill[i].clear();
    }
  } catch(std::exception & lExc) {
    std::cerr << lExc.what() << std::endl
	      << "Standard exception\n"
	      << "Spill skipped!\n\n";
    for(unsigned int i = 0; i<48; i++)// for tempHitsPerSpill.size()
    { 
      tempHitsPerSpill[i].clear();
    }
  } catch(...) {
    std::cerr << "Unknown exception occurred...\n"
	      << "Spill skipped!\n\n";
    for(unsigned int i = 0; i<48; i++)// for tempHitsPerSpill.size()
    { 
      tempHitsPerSpill[i].clear();
    }
  }


  // Use old to fill per trigger. Check if trigger has more than 12 hits.(4 left, 4right and 4 vertical)

  // Fill vector with first FEB hits (in some sense a trigger)
  /*
for(unsigned int i = 0; i<48; i++)// for tempHitsPerSpill.size()
    {
      cerr<<"tempHitsPerSpill["<<i<<"].size()="<<tempHitsPerSpill[i].size()<<endl;
    }
  */

  vector<pair<double,vector<TempHit*> > >temp1; 
  vector<pair<double,vector<TempHit*> > >temp2; 
  
  for(unsigned int i = 0; i< tempHitsPerSpill[0].size(); i++)
    {
      vector<TempHit*> tempVector;
      tempVector.assign(tempHitsPerSpill[0][i].second.begin(),tempHitsPerSpill[0][i].second.end());
      //tempVector.insert(tempVector.end(),tempHitsPerSpill[0][i].second.begin(),tempHitsPerSpill[0][i].second.end());
      temp1.push_back(make_pair(tempHitsPerSpill[0][i].first,tempVector));
    }
  //cerr<<"temp1.size()="<<temp1.size()<<endl;

  for(unsigned int i = 0; i< temp1.size(); i++)
    {
      vector<TempHit*> tempVector;
      for(unsigned int j = 1; j< 48; j++)
	{
	  for(unsigned int k = 0; k< tempHitsPerSpill[j].size(); k++)
	    {
	      //cout<<j<<"\tGtrigg time\t"<<temp1[i].first<<"\t"<<tempHitsPerSpill[j][k].first<<endl;
	      if(temp1[i].first == tempHitsPerSpill[j][k].first)
		//cout<<"match"<<endl;
		//if(fabs(temp1[j].first-tempHitsPerSpill[i][k].first)<3)
		{
		  tempVector.insert(tempVector.end(),tempHitsPerSpill[j][k].second.begin(),tempHitsPerSpill[j][k].second.end());
		  tempVector.insert(tempVector.end(),temp1[i].second.begin(),temp1[i].second.end());
		  //cout<<i<<"\t"<<j<<"\t"<<k<<endl;
		  //cout<<temp1[i].second.size()<<"\t"
		  //   <<tempHitsPerSpill[j][k].second.size()<<"\t"
		  //  <<tempVector.size()<<endl;
		}
	    }
	}
      temp2.push_back(make_pair(temp1[i].first,tempVector));
    } 
  temp1.clear();

  //cerr<<tempVector.size()<<endl;
  for(unsigned int i = 0; i< temp2.size(); i++)
    {
      //cout<<"temp2 i="<<i<<"\t"<<temp2[i].first<<"\t"<<temp2[i].second.size()<<endl;
      vector<TempHit*> tempVector;
      tempVector.clear();
      tempVector.assign(temp2[i].second.begin(),temp2[i].second.end());
      //tempVector.insert(tempVector.end(),temp1[i].second.begin(),temp1[i].second.end());
      
      // Check board 1,2,3 and 4 have hits
      //left and right
      bool checker = false;
      if(tempVector.size()>12)
	{
	  bool checker0 = false;
	  bool checker1 = false;
	  bool checker2 = false;
	  bool checker3 = false;
	  bool checker4 = false;
	  bool checker8 = false;
	  bool checker9 = false;
	  bool checker10 = false;
	  bool checker11 = false;
	  bool checker12 = false;
	  bool checker27 = false;
	  for(unsigned int j = 0; j<tempVector.size(); j++)
	    {
	      if(tempVector[j]->boardId==0) checker0=true; //mod0 left
	      else if(tempVector[j]->boardId==1) checker1=true; //mod1 left
	      else if(tempVector[j]->boardId==2) checker2=true; //mod2 left
	      else if(tempVector[j]->boardId==3) checker3=true; //mod3 left
	      else if(tempVector[j]->boardId==4) checker4=true; //mod4 left
	      else if(tempVector[j]->boardId==8) checker8=true; //mod8 left
	      else if(tempVector[j]->boardId==9) checker9=true; //mod1 right
	      else if(tempVector[j]->boardId==10) checker10=true; //mod2 right
	      else if(tempVector[j]->boardId==11) checker11=true; //mod3 right
	      else if(tempVector[j]->boardId==12) checker12=true; //mod4 right
	      else if(tempVector[j]->boardId==27) checker27=true; //vertical 123
	    }
	  checker = checker0 && checker1 && checker2
	    && checker3 && checker4 
	    && checker8 && checker9 && checker10
	    && checker11 && checker12 && checker27;
	}

      //cerr<<tempVector.size()<<endl;
      if(tempVector.size()>12 && checker)
	{
	  //cerr<<"pushing back"<<endl;
	  vectorFinalHits->push_back(tempVector);
	  /*
	  cout<<"temp2[i].first="<<temp2[i].first<<"\t"
	      <<"temp2[i].second.size()="<<temp2[i].second.size()<<endl;
	  cerr<<"temp2[i].first="<<temp2[i].first<<"\t"
	      <<"temp2[i].second.size()="<<temp2[i].second.size()<<endl; 
	  */ 
	  /*for(unsigned int j = 0; j<tempVector.size(); j++)
	    {
	      cerr<<tempVector[j]->boardId<<endl;
	      cout<<tempVector[j]->boardId<<endl;
	    }
	  */

	}
	    
      tempVector.clear();
    }
  /*
  cerr<<"vectorFinalHits->size()="<<vectorFinalHits->size()<<endl;
  for(unsigned int i = 0; i<vectorFinalHits->size();i++)
    {
      cerr<<"vectorFinalHits["<<i<<"].size()="<<vectorFinalHits[i].size()<<endl;
    }
  */
}

void HandleData(vector<string> filenames, string filepath,root2dst* cvt)
{
  vector<MDdateFile*> dfiles;
  vector<char*> eventBuffers;
  int evt_read = 0;
  
  for(unsigned int fileInt=0; fileInt<filenames.size();fileInt++)
    {
      dfiles.push_back(new MDdateFile(filenames[fileInt],filepath));
      eventBuffers.push_back(new char());
    }
  bool validFiles = true;
  for(unsigned int fileInt=0;fileInt<dfiles.size();fileInt++)
    {
      validFiles &= dfiles[fileInt]->open();
    }
  if(validFiles){
    for(unsigned int fileInt=0;fileInt<dfiles.size();fileInt++)
      {
	dfiles[fileInt]->init();
      }
    int counter =0;
    do { // Loop over all spills
      vector<vector<TempHit*> > tempVectorHits;
      for(unsigned int i=0;i<filenames.size();i++)
	{
	  eventBuffers[i] = dfiles[i]->GetNextEvent();
	}
      
      //if(counter == 2) break; //beam test 2017
      //if(counter == 30) break; //beam test 2017
      /*
	if(counter>20 && counter<30)
	{
	counter++;
	continue;
	}
      */
      if(counter == 221) break; //beam test 2017
      //if(counter == 78) break; //beam test 2017
      //counter++;
      //cerr<<"convert_to_internal2017BeamTest3 start"<<endl;
      convert_to_internal(eventBuffers,&tempVectorHits);
      //cerr<<"convert_to_internal2017BeamTest3 end"<<endl;
      //int iEvent = 0;
      //continue;
      if(tempVectorHits.size())
	{ 
	  counter++;
	  cerr<<"counter="<<counter<<endl;
	  
	  for(unsigned int event=0;event<tempVectorHits.size();event++)
	    {
	      //cerr<<"counter="<<counter<<endl;
	      //counter++;
	      if(tempVectorHits[event].size()<12) continue;
	      
	      cout<<"eventNum="<<event<<endl;
	      string filepath2="/data/neutrino05/phallsjo/SaRoMan";
	      dataBaseVector.clear();
	      //cerr<<"ConvertToHit2017BeamTest3 start"<<endl;
	      vector<bhep::hit*> unSortedhits = ConvertToHit(tempVectorHits[event],filepath+"/BabyMIND.db");
	      //cerr<<"ConvertToHit2017BeamTest3 end"<<endl;	  
	      std::vector<bhep::hit*> hits = unSortedhits;
	      //cerr<<"sorting "<<hits.size()<<endl;
	      sort( hits.begin(), hits.end(), forwardSortInZ() );
	      //cerr<<"sorting done"<<endl;
	      ptype pT = DIGI;
	      string detect = "tracking";
	      
	      if(hits.size())
		{
		  
	      vector<particle*> hitsParticle;
	      hitsParticle.push_back(new particle(pT,detect));
	      
	      for(unsigned int i = 0; i<hits.size();i++)
		{
		  hitsParticle.back()->add_hit(detect,hits[i]); //new
		}
	      //cerr<<"particles size="<<hitsParticle.back()->hits(detect).size()<<endl;
	      //cout<<"particles size="<<hitsParticle.back()->hits(detect).size()<<endl;
	      bhep::event e(evt_read);
	      e.add_property( "IntType", "CCQE" );
	      e.add_property("G4EventID",evt_read); 
	      
	      //cerr<<"gdml hit constructor"<<endl;
	      particle* digi_part = cvt->create_digital_representation( hitsParticle );
	      
	      if(digi_part->hit_map().size())
		{
		  e.add_digi_particle( digi_part );
		  //cerr<<"gdml hit constructor end"<<endl;
		  
		  
		  //cerr<<"Writing root"<<endl;
		  bhep::bhep_svc::instance()->get_writer_root().write( e, evt_read );//iEvent );
		  //cerr<<"Writing root end"<<endl;
		  evt_read++;
		  
		  //cout<<e<<endl;
		}
	      e.clear();
	      hitsParticle.clear(); 
	    }	
	    }
	}
      tempVectorHits.clear();
    }while(!NullCheck(eventBuffers));
  }
  for(unsigned int fileInt=0;fileInt<dfiles.size();fileInt++)
    {
      dfiles[fileInt]->close();
    } 
}

int main(int argc, char* argv[]) {
  
  //MDargumentHandler argh("Example of unpacking application.");
 
 
  string param_file;
  long rndmSeed;

  if (argc == 2) param_file = argv[1];
  else {

    std::cout << "Execute as: ./simple_smear <paramater file>" << std::endl;

    return -1;
  }

  //Stores for information from parameter file.
  bhep::gstore run_store, data_store;
  bhep::gstore* con_store = new bhep::gstore();
  
  //Read run parameters.
  bhep::sreader reader1(run_store);
  reader1.file(param_file);
  reader1.group("RUN");
  reader1.read();
  //
  //Read files to be processed.
  bhep::sreader reader2(data_store);
  reader2.file(param_file);
  reader2.group("DATA");
  reader2.read();
  //hit_constructor stuff
  bhep::sreader reader3(*con_store);
  reader3.file(param_file);
  reader3.group("CON");
  reader3.read();
  //
  
  int nEvents;
  double smearRes[2];

  if ( !run_store.find_istore("nEvents") ) {
    std::cout << "Parameter file must contain number of events "
	      << "to be processed in group RUN" << std::endl;
    return -1;
  }
  else nEvents = run_store.fetch_istore("nEvents");

  if ( !run_store.find_dstore("Gaus_Sigma") ) {
    std::cout << "Parameter file must contain smear sigma in cm "
	      << "as double Gaus_Sigma in group RUN" << std::endl;
    return -1;
  }
  else smearRes[0] = run_store.fetch_dstore("Gaus_Sigma");

  if ( !run_store.find_dstore("Eng_Res") ) {
    std::cout << "Parameter file must contain Energy resolution"
	      << "as double Eng_Res in group RUN" << std::endl;
    return -1;
  }
  else smearRes[1] = run_store.fetch_dstore("Eng_Res");

  if ( !con_store->find_dstore("Gen_seed") ) {
    std::cout << "Parameter file must contain generator seed "
	      << "as double Gen_seed in group RUN" << std::endl;
    return -1;
  }
  else rndmSeed = (long)con_store->fetch_dstore("Gen_seed");
  
  //EventManager2* eman = new EventManager2(data_store, bhep::NORMAL);
  //Define and open input.
  reader_root inDst;
  //inDst.open( data_store.fetch_sstore("idst_file") );
   
  //Define output.
  WriteUtil::OpenOutputDst( data_store.fetch_sstore("odst_file") );
  
  root2dst* cvt = new root2dst(bhep::NORMAL, con_store);
  
  //eman->initialize();
  
  cvt->initialize( smearRes, rndmSeed);

  //Counters.
  int iEvent;
  int evt_read = 0;

  vector<string> input_data = data_store.fetch_svstore("idst_files");


  int testBeam = run_store.fetch_istore("testBeam");

  //vector<TempHit*> tempHits; // Vector for spill

  //vector<vector<TempHit*> > tempVectorHits; // Vector for spills per trigger.

  //vector<vector<vector<TempHit*> > > tempVectorVectorHits;


  if(testBeam)
    {
      // Want to read n number files at a time
      // From, for each spill
      // Add hits togeather using the propper time. Trigg+hit.
      // Or should it be per trigger? Since simulation per particle.

      // How do the event builder?
      //char *dataBuff;
      //uint32_t* dataPtr;

      string filepath;
      string filename;
      string filename2;
      string filename3;
      string filename4;
      
      TFile rfile("histos.root", "recreate");
      TH1I  h1("h1_hit_ch", "hit channels", 100000, 0, 100000);
      TH1I  h2("h2_hit_ch", "hit channels", 100000, 0, 100000);

      TH1I  h3("h3_hit_ch", "hit channels", 100, 0, 100);
      TH1I  h4("h4_hit_ch", "hit channels", 100, 0, 100);

      TH1I  h5("h5_hit_ch", "hit channels", 16, -120, 120);
      TH1I  h6("h6_hit_ch", "hit channels", 16, -120, 120);

      TH1I  h7("h7_hit_ch", "hit channels", 100, 0, 100);
      TH1I  h8("h8_hit_ch", "hit channels", 100, 0, 100);

      TH2D  beamHistoHG1("beamHistoHG1", "beamHistoHG1", 16, 0, 16,16,0,16);
      TH2D  beamHistoHG2("beamHistoHG2", "beamHistoHG2", 160, 0, 160,160,0,160);


      //TH1I  h2("h2_ampl", "hit ampl.", 200, 0, 5000);
      

      //filepath="/data/neutrino05/phallsjo/SaRoMan";
      //filename="feb1_std_10gevmuons_extscint16_nogarbage_refdata4.daq";
      //filename="FEB1_safe_mode_3_HG40_LG55_test1.daq";
      //filename="A2-all-10KHz-EXTPULSE_10-000us-ANALOG-1.daq";
      //filename="test.daq";

      // Do files in parallell instead of in series. Check hits in all FEBs at the same time && both channels at same time.
      //ForFile(filename,filepath,&tempHits,&tempVectorHits);

      //filepath="/data/neutrino05/phallsjo/copy/SaRoMan/daq";
      filepath="/data/neutrino06/phallsjo/daq";

      vector<string> tester;
      string sFileName;
      //ifstream fList("/data/neutrino05/phallsjo/copy/SaRoMan/module.list");
      ifstream fList("/data/neutrino06/phallsjo/daq/module.list");
      while (!fList.eof()) {
        fList >> sFileName;
        cerr << sFileName << endl;
        tester.push_back(sFileName);
      }

      //vector<string> tester;
      //tester.push_back(filename);
      //tester.push_back(filename2);
      //tester.push_back(filename3);
      //tester.push_back(filename4);

      //char *dataBuff;
      //uint32_t* dataPtr;

      //cerr<<"start of handledata"<<endl;
      //cout<<"start of handledata"<<endl;
      //HandleData2(tester,filepath,cvt,&h1);
      HandleData(tester,filepath,cvt);
      //cerr<<"End of handledata"<<endl;
      //cout<<"End of handledata"<<endl;
      
      //h1.Write();
      //h2.Write();
      rfile.Write();
      rfile.Close();
      
    }
  else
    {
      for (unsigned int ifile = 0;ifile < input_data.size();ifile++){
	
	inDst.open( input_data[ifile] );
	iEvent = 0;
	
	//  for (int iEvent = 0;iEvent < nEvents;iEvent++) {
	while( !inDst.eof( iEvent ) && evt_read < nEvents ){
	  // std::cout << "Event: " << evt_read << std::endl;
	  
	  bhep::event& e = inDst.read_event( iEvent );//eman->read();
	  std::cout << "Event: " << evt_read <<", "<<e.event_number()<< std::endl;
	  //Get Vector of all true particles from event.
	  vector<bhep::particle*> particles = e.true_particles();
	  //
	  
	  particle* digi_part = cvt->create_digital_representation( particles );
	  
	  e.add_digi_particle( digi_part );

	  //cout<<e<<endl;
	  
	  bhep::bhep_svc::instance()->get_writer_root().write( e, evt_read );//iEvent );
	  e.clear();
	  
	  iEvent++;
	  evt_read++;
	}
	
	inDst.close();
      }
      cvt->print(); // Used to print the debug histograms to a root file.
    }
 

  WriteUtil::CloseOutputDst();
  //inDst.close();

  return 0;
}
