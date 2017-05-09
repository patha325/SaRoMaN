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


void ConvertToPosition(bhep::hit* inHit, unsigned int boardId, unsigned int channelId, Xml_parser* parse, xmlNode* root_element)
{
  //<channelMap boardID="0"  channelID="9"  barNumber="128"  barPosZ="0"  isYBar="0"  moduleNum="1"  barPosT="0"/>

  //cout<<"HERE"<<endl;

  /*
    Using channelId and boardId, from the parsedXML file return the position data.
  */
  double barPosZ =0.0;
  int isYBar=0;
  int barNumber=0;
  double barPosT=0.0;
  int moduleNum=0;

  xmlNode *cur_node = parse->GetNode(root_element,"boardID",boardId,"channelID",channelId);
  //cout<<root_element->name()<<endl;
  //cout<<root_element->value()<<endl;
  //cout<<boardId<<endl;
  //cout<<channelId<<endl;

  if(cur_node != NULL)
    {
      xmlAttr* attribute = cur_node->properties;
      
      while(attribute)
	{
	  //cout<<"Attribute name="<<(char*)attribute->name<<endl;
	  if(strncmp ((char*)attribute->name,"barPosZ",7) == 0)
	    {
	      barPosZ = atof((char*)xmlNodeGetContent(attribute->children));
	    }
	  else if(strncmp ((char*)attribute->name,"isYBar",6) == 0)
	    {
	      isYBar = atoi((char*)xmlNodeGetContent(attribute->children));
	      //cout<<"boardID="<<boardId<<endl;
	      //cout<<"isYBar="<<isYBar<<endl;
	    }
	  else if(strncmp ((char*)attribute->name,"barNumber",9) == 0)
	    {
	      barNumber = atoi((char*)xmlNodeGetContent(attribute->children));
	    }
	  else if(strncmp ((char*)attribute->name,"barPosT",7) == 0)
	    {
	      barPosT = atof((char*)xmlNodeGetContent(attribute->children));
	    }
	  else if(strncmp ((char*)attribute->name,"moduleNum",9) == 0)
	    {
	      moduleNum = atoi((char*)xmlNodeGetContent(attribute->children));
	    }
	  attribute = attribute->next;
	}
      attribute =0;
    }
  else
    {
      cout<<"Node is NULL!"<<endl;
    }


  if(isYBar)
    {
      barNumber = (barNumber -118)*2;
    }
  else
    {
      barNumber = (barNumber - 34)*2 -1;
    }

  inHit->add_property("barPosZ",barPosZ);
  inHit->add_property("IsYBar",isYBar);
  inHit->add_property("barNumber",barNumber);
  inHit->add_property("barPosT",barPosT);
  inHit->add_property("moduleNum",moduleNum);
  // can we add ->x()[2] how is it done? bhep add?
  if(isYBar) inHit->set_point(*(new bhep::Point3D(0,barPosT,barPosZ)));
  else inHit->set_point(*(new bhep::Point3D(barPosT,0,barPosZ)));

  cur_node= 0;

  //void set_point(const bhep::Point3D&  x)

  //point(double x, double y, double z, string view="XYZ"):

  //void set_point(const bhep::Point3D&  x)
}

double ConvertToEdep(unsigned int boardId, unsigned int channelId,unsigned int amplitude)
{
  return 100.0;
}

vector<bhep::hit*>ConvertToHit(vector<TempHit*> hitsVector, string xmlName)
{
  /*
    Parses the xml file and uses the information to convert from the vector of TempHit to
    a vector of vectors of bhep hit*. The vectors are split depending on the different gtrigs.
  */
  Xml_parser parse(xmlName);
  xmlNode *root_element = parse.ParseXML();
  vector<bhep::hit*>returnVector;

  int counterFEB18 = 0;
  int counterFEB38 = 0;

  int counterFEB10 = 0;
  int counterFEB30 = 0;
  int counterFEB115 = 0;
  int counterFEB315 = 0;


  for(unsigned int cnt = 0; cnt<hitsVector.size(); cnt++)
  {
    returnVector.push_back(new bhep::hit());
    double tempEdep = ConvertToEdep(hitsVector[cnt]->boardId,hitsVector[cnt]->channelId,hitsVector[cnt]->amplitude);
    returnVector.back()->add_property("EnergyDep",tempEdep);
    returnVector.back()->add_property("time",(double)hitsVector[cnt]->timeing);

    if((15 < hitsVector[cnt]->channelId && hitsVector[cnt]->channelId < 32) ||
       (47 < hitsVector[cnt]->channelId && hitsVector[cnt]->channelId < 64) ||
       (79 < hitsVector[cnt]->channelId && hitsVector[cnt]->channelId < 96))
      {
	returnVector.back()->add_property("rHit",1);
      }
    else
      {
	returnVector.back()->add_property("rHit",0);
      }

    /*
    if(hitsVector[cnt]->boardId == 1 && hitsVector[cnt]->channelId == 8) counterFEB18++;
    if(hitsVector[cnt]->boardId == 3 && hitsVector[cnt]->channelId == 8) counterFEB38++;

    if(hitsVector[cnt]->boardId == 1 && hitsVector[cnt]->channelId == 0) counterFEB10++;
    if(hitsVector[cnt]->boardId == 3 && hitsVector[cnt]->channelId == 0) counterFEB30++;

    if(hitsVector[cnt]->boardId == 1 && hitsVector[cnt]->channelId == 15) counterFEB115++;
    if(hitsVector[cnt]->boardId == 3 && hitsVector[cnt]->channelId == 15) counterFEB315++;
    */
    ConvertToPosition(returnVector.back(),hitsVector[cnt]->boardId,hitsVector[cnt]->channelId,&parse,root_element);
  }
  /*
  cout<<"counterFEB18="<<counterFEB18<<endl;
  cout<<"counterFEB38="<<counterFEB38<<endl;

  cout<<"counterFEB10="<<counterFEB10<<endl;
  cout<<"counterFEB30="<<counterFEB30<<endl;

  cout<<"counterFEB115="<<counterFEB115<<endl;
  cout<<"counterFEB315="<<counterFEB315<<endl;
  */

  root_element = 0;

  return returnVector;
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
void ConvertToPosition2017BeamTest(bhep::hit* inHit, unsigned int boardId, unsigned int channelId)
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
  double mod1ZPos = -detectorZ + SciDepth/2;
  double mod2ZPos = mod1ZPos + SciDepth/2 + 2*FeDepth + 2*GapFeFe +2*GapFeSci + SciDepth/2;
  double mod3ZPos = mod2ZPos + SciDepth/2 + Gap + SciDepth/2;
  double mod4ZPos = mod3ZPos + SciDepth/2 + 2*FeDepth + 2*GapFeFe +2*GapFeSci + SciDepth/2;
  double mod5ZPos = mod4ZPos + SciDepth/2 + Gap + SciDepth/2;
  double mod6ZPos = mod5ZPos + SciDepth/2 + FeDepth + 2*GapFeSci + SciDepth/2;
  double mod7ZPos = mod6ZPos + SciDepth/2 + FeDepth + 2*GapFeSci + SciDepth/2; // detectorZ-FeDepth- GapFeSci - SciDepth/2;
  double moduleList[7] = {mod1ZPos, mod2ZPos, mod3ZPos, mod4ZPos, mod5ZPos, mod6ZPos, mod7ZPos};
  double xBarZsize = 7.5 *2; //xYYx
  double yBarZsize = 7.5;
  
  for(int i=0;i<dataBaseVector.size();i+=6)
    {
      int Global_channel = myatoi(dataBaseVector[i]);
      //cout<<"if_check"<<Global_channel<<"\t"<<channelId<<endl;
      if(Global_channel==channelId)
	{
   	string Bar_type = dataBaseVector[i+1];
   	int Module = myatoi(dataBaseVector[i+2]);
   	string Position_label = dataBaseVector[i+3];
   	string Bar_position = dataBaseVector[i+4];
   	string Side = dataBaseVector[i+5];

	double barPosZ =0.0;
	int YMeasuring=0;
	int barNumber=0;
	double barPosT=0.0;
	//int moduleNum=1;
	  
	//barPosZ = mod1ZPos;

 	if((myatoi(Position_label[1])==2 || myatoi(Position_label[1])==4) && Side =="right") Side = "left";
	else if((myatoi(Position_label[1])==2 || myatoi(Position_label[1])==4) && Side =="left") Side = "right";	
	if(Bar_type=="vertical")
	  {
	    if(Side =="left")
	      {
		barPosT = (Global_channel -1568 +32 - 32 * myatoi(Position_label[1]))*yBarPitch +yBarStart;
		barNumber = (Global_channel -1568 +32 - 32 * myatoi(Position_label[1]));
	      }	
	    else
	      {
		barPosT = (Global_channel -1568-1 +32 - 32 * myatoi(Position_label[1]))*yBarPitch +yBarStart;
		barNumber = (Global_channel -1568-1 +32 - 32 * myatoi(Position_label[1]));
	      }
	    YMeasuring=0;
	  }
	else
	  {
	    if(Side =="left")
	      { 
		barPosT = (Global_channel - 96 * myatoi(Position_label[1]))*xBarPitch +xBarStart;
		barNumber = (Global_channel - 96 * myatoi(Position_label[1]));
	      }
	    else
	      {
		barPosT = (Global_channel - 672 +96 - 96 * myatoi(Position_label[1]))*xBarPitch +xBarStart;
		barNumber =  (Global_channel - 672 +96 - 96 * myatoi(Position_label[1]));
	      }
	    YMeasuring=1;
	  }
	/*
	for (int moduleNum = 1;moduleNum<=7;moduleNum++)
	  {
	    if(Bar_type=="vertical" && myatoi(Position_label[1]) == moduleNum && Bar_position =="front") barPosZ = moduleList[moduleNum-1]-xBarZsize;
	    else if(Bar_type=="vertical" && myatoi(Position_label[1]) == moduleNum && Bar_position =="back") barPosZ = moduleList[moduleNum-1]+xBarZsize;
	    else if(Bar_type=="horizontal" && myatoi(Position_label[1]) == moduleNum && Bar_position =="front") barPosZ = moduleList[moduleNum-1]-yBarZsize;
	    else if(Bar_type=="horizontal" && myatoi(Position_label[1]) == moduleNum && Bar_position =="back") barPosZ = moduleList[moduleNum-1]+yBarZsize;
	  }
	*/

	//inHit->add_property("barPosZ",barPosZ);
	inHit->add_property("barPosZ",mod1ZPos);
	  inHit->add_property("IsYBar",YMeasuring);
	  inHit->add_property("barNumber",barNumber);
	  inHit->add_property("barPosT",barPosT);
	  inHit->add_property("IsTASD",0);
	  inHit->add_property("moduleNum",myatoi(Position_label[1]));
	  // can we add ->x()[2] how is it done? bhep add?
	  if(YMeasuring) inHit->set_point(*(new bhep::Point3D(0,barPosT,barPosZ)));
	  else inHit->set_point(*(new bhep::Point3D(barPosT,0,barPosZ)));
	}
      else
	{
	  /*
	  inHit->add_property("barPosZ",0);
	  inHit->add_property("IsYBar",0);
	  inHit->add_property("barNumber",0);
	  inHit->add_property("barPosT",0);
	  inHit->add_property("moduleNum",0);
	  inHit->set_point(*(new bhep::Point3D(0,0,0)));
	  */
	}
    }
 // Leaving convert to position2017
}

vector<bhep::hit*> ConvertToHit2017BeamTest(vector<TempHit*> hitsVector, string dataBaseName)
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
    returnVector.push_back(new bhep::hit());
    double tempEdep = ConvertToEdep(hitsVector[cnt]->boardId,hitsVector[cnt]->channelId,hitsVector[cnt]->amplitude);
    returnVector.back()->add_property("EnergyDep",tempEdep);
    returnVector.back()->add_property("time",(double)hitsVector[cnt]->timeing);

    if(hitsVector[cnt]->boardId==3 || (hitsVector[cnt]->boardId == 5 && hitsVector[cnt]->channelId % 2))
      {
	returnVector.back()->add_property("rHit",1);
      }
    else
      {
	returnVector.back()->add_property("rHit",0);
      }

    ConvertToPosition2017BeamTest(returnVector.back(),hitsVector[cnt]->boardId,hitsVector[cnt]->channelId);
 
  }

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

void convert_to_internal2017BeamTest(vector<char*> eventBuffers,vector<vector<TempHit*> >* tempVectorHits,TH1I* h1)
{
  vector<pair<double,vector<TempHit*> > > tempHitsPerSpill[3]; 
  try {
    for(unsigned int eventNum=0;eventNum<eventBuffers.size();eventNum++)
      {
	vector<TempHit*> tempVec;
	
	MDfragmentBM   spill;
	spill.SetDataPtr(eventBuffers[eventNum]);
	
	//MDpartEventBM *event;
	int nTr = spill.GetNumOfTriggers();

	for (unsigned int triggNum=0; triggNum<nTr; ++triggNum) {
	  MDpartEventBM * event = spill.GetTriggerEventPtr(triggNum);
	  if(event->getNumDataWords() < 12) continue;

	  
	  for (int ich=0; ich<BM_FEB_NCHANNELS; ++ich) {
	    int nHits = event->GetNLeadingEdgeHits(ich);
	    
	    int otherSide = false;

	    otherSide= event->GetNLeadingEdgeHits(ich+BM_FEB_NCHANNELS);

	    if (nHits)
	      {
		if(event->GetHitTime(0,ich, 'l')>0)
		  {
		    h1->Fill(triggNum);
		  }
	      }
	    
	    if(nHits && otherSide)
	      {
		for(int ih=0; ih<nHits; ih++)
		  {
		    if(event->GetHitTime(ih,ich, 'l')>0)
		      {
			//h1->Fill(triggNum);
			//counter++;
			
			//2017BeamTest
			int GlobalChannel=ich;
			if(spill.GetBoardId()== 4) GlobalChannel+=96;
			if(spill.GetBoardId()== 3) GlobalChannel+=672;
			if(spill.GetBoardId()== 5) GlobalChannel+=1568;

			tempVec.push_back(new TempHit(event->GetHitAmplitude(ich, 'h'),
						      event->GetHitTime(ih,ich, 'l'),
						      GlobalChannel,spill.GetBoardId()));
			//tempPairVec.push_back(make_pair(events[i]->GetTriggerTime(),
			//			  (new TempHit(events[i]->GetHitAmplitude(ich, 'h'),
			//				       events[i]->GetHitTime(ih,ich, 'l'),
			//				       ich,spills[i]->GetBoardId()))));
		      }
		  }
	      }
	  }//END for channels
	  //tempHitsPerSpill[i].push_back(make_pair(triggNum,tempVec));
	  tempHitsPerSpill[eventNum].push_back(make_pair(event->GetTriggerTime(),tempVec));
	  tempVec.clear();
	}//END of trigger
	spill.Clean();
      } //END of eventNum;
  } catch (MDexception & lExc)  {
    std::cerr <<  lExc.GetDescription() << endl
	      << "Unpacking exception\n"
	      << "Spill skipped!\n\n";
  } catch(std::exception & lExc) {
    std::cerr << lExc.what() << std::endl
	      << "Standard exception\n"
	      << "Spill skipped!\n\n";
  } catch(...) {
    std::cerr << "Unknown exception occurred...\n"
	      << "Spill skipped!\n\n";
  }
  
  cerr<<"tempHitsPerSpill[0].size()="<<tempHitsPerSpill[0].size()<<endl;
  cerr<<"tempHitsPerSpill[1].size()="<<tempHitsPerSpill[1].size()<<endl;
  cerr<<"tempHitsPerSpill[2].size()="<<tempHitsPerSpill[2].size()<<endl;
  //cerr<<"tempHitsPerSpill[3].size()="<<tempHitsPerSpill[3].size()<<endl;
  
  vector<pair<double,vector<TempHit*> > >temp1; 
  vector<pair<double,vector<TempHit*> > >temp2; 

  for(unsigned int i = 0; i< tempHitsPerSpill[0].size(); i++)
    {
      for(unsigned int j =0; j<tempHitsPerSpill[1].size(); j++)
	{
	  if(tempHitsPerSpill[0][i].first == tempHitsPerSpill[1][j].first)
	    {
	      vector<TempHit*> tempVector;
	      tempVector.insert(tempVector.end(),tempHitsPerSpill[0][i].second.begin(),tempHitsPerSpill[0][i].second.end());
	      tempVector.insert(tempVector.end(),tempHitsPerSpill[1][j].second.begin(),tempHitsPerSpill[1][j].second.end());
	      temp1.push_back(make_pair(tempHitsPerSpill[0][i].first,tempVector));
	      //cerr<<"match="<<tempHitsPerSpill[0][i].first<<endl;
	    }
	}
    }

  for(unsigned int i = 0; i< tempHitsPerSpill[2].size(); i+=2)
    {
      //if(tempHitsPerSpill[2][i].first == tempHitsPerSpill[2][i+1].first)
      //{
	  vector<TempHit*> tempVector;
	  tempVector.insert(tempVector.end(),tempHitsPerSpill[2][i].second.begin(),tempHitsPerSpill[2][i].second.end());
	  //tempVector.insert(tempVector.end(),tempHitsPerSpill[2][i+1].second.begin(),tempHitsPerSpill[2][i+1].second.end());
	  temp2.push_back(make_pair(tempHitsPerSpill[2][i].first,tempVector));
	  //cerr<<"match="<<tempHitsPerSpill[2][i].first<<endl;
	  //}
    }
  cerr<<"temp1.size()="<<temp1.size()<<endl;
  cerr<<"temp2.size()="<<temp2.size()<<endl;
  
  for(unsigned int i = 0; i< temp1.size(); i++)
    {
      for(unsigned int j =0; j<temp2.size(); j++)
	{
	  
	  if(temp1[i].first == temp2[j].first)
	    {
	      vector<TempHit*> tempVector;
	      tempVector.insert(tempVector.end(),temp1[i].second.begin(),temp1[i].second.end());
	      tempVector.insert(tempVector.end(),temp2[j].second.begin(),temp2[j].second.end());
	      tempVectorHits->push_back(tempVector);
	      //cerr<<"match="<<tempHitsPerSpill[2][i].first<<endl;
	    }
	  
	}
    }
  
  cerr<<tempVectorHits->size()<<endl;
}


void convert_to_internal(vector<char*> eventBuffers,
	  vector<vector<TempHit*> >* tempVectorHits,
			 TH1I* h1, double* noCut, double* firstCut,double* secondCut,double* thirdCut,double* forthCut,double* fifthCut)
{
  vector<pair<double,vector<TempHit*> > > tempHitsPerSpill[4]; 
  try {
    for(unsigned int eventNum=0;eventNum<eventBuffers.size();eventNum++)
      {
	vector<TempHit*> tempVec;
	
	MDfragmentBM   spill;
	spill.SetDataPtr(eventBuffers[eventNum]);
	
	//MDpartEventBM *event;
	int nTr = spill.GetNumOfTriggers();
	
	int asic1=0;
	int asic2=0;
	int asic3=0;
	
	for (unsigned int triggNum=0; triggNum<nTr; ++triggNum) {
	  MDpartEventBM * event = spill.GetTriggerEventPtr(triggNum);
	  
	  //if(event->getNumDataWords() < 2+4*6) continue;
	  (*noCut)++;
	  if(event->getNumDataWords() < 12) continue;
	  (*firstCut)++;

	  
	  for (int ich=0; ich<BM_FEB_NCHANNELS; ++ich) {
	    int nHits = event->GetNLeadingEdgeHits(ich);
	    
	    int otherSide = false;

	    if(ich<16) otherSide= event->GetNLeadingEdgeHits(ich+16);

	    if(15<ich&&ich<32) otherSide= event->GetNLeadingEdgeHits(ich-16);
	    
	    if(31<ich&&ich<48) otherSide= event->GetNLeadingEdgeHits(ich+16);

	    if(47<ich&&ich<64) otherSide= event->GetNLeadingEdgeHits(ich-16);

	    if(63<ich&&ich<80) otherSide= event->GetNLeadingEdgeHits(ich+16);

	    if(79<ich&&ich<96) otherSide= event->GetNLeadingEdgeHits(ich-16);


	    if (nHits)
	      {
		if(event->GetHitTime(0,ich, 'l')>0)
		  {
		    h1->Fill(triggNum);
		  }
		(*secondCut)++;
	      }
	    
	    if(nHits && otherSide)
	      {
		(*thirdCut)++;
		for(int ih=0; ih<nHits; ih++)
		  {
		    if(event->GetHitTime(ih,ich, 'l')>0)
		      {
			//h1->Fill(triggNum);
			if(ich<32) asic1++;
			
			if(31<ich && ich<64) asic2++;
			
			if(63<ich && ich<96) asic3++;
			
			//counter++;
			tempVec.push_back(new TempHit(event->GetHitAmplitude(ich, 'h'),
						      event->GetHitTime(ih,ich, 'l'),
						      ich,spill.GetBoardId()));
			//tempPairVec.push_back(make_pair(events[i]->GetTriggerTime(),
			//			  (new TempHit(events[i]->GetHitAmplitude(ich, 'h'),
			//				       events[i]->GetHitTime(ih,ich, 'l'),
			//				       ich,spills[i]->GetBoardId()))));
		      }
		  }
	      }
	  }//END for channels
	  
	  (*forthCut)++;
	  if(asic1 > 1 && asic2 > 1 && asic3 > 1) 
	    (*fifthCut)++;
	    //tempHitsPerSpill[i].push_back(make_pair(triggNum,tempVec));
	    tempHitsPerSpill[eventNum].push_back(make_pair(event->GetTriggerTime(),tempVec));
	  tempVec.clear();
	}//END of trigger
spill.Clean();
      } //END of eventNum;
  } catch (MDexception & lExc)  {
    std::cerr <<  lExc.GetDescription() << endl
	      << "Unpacking exception\n"
	      << "Spill skipped!\n\n";
  } catch(std::exception & lExc) {
    std::cerr << lExc.what() << std::endl
	      << "Standard exception\n"
	      << "Spill skipped!\n\n";
  } catch(...) {
    std::cerr << "Unknown exception occurred...\n"
	      << "Spill skipped!\n\n";
  }
  
  cerr<<"tempHitsPerSpill[0].size()="<<tempHitsPerSpill[0].size()<<endl;
  cerr<<"tempHitsPerSpill[1].size()="<<tempHitsPerSpill[1].size()<<endl;
  cerr<<"tempHitsPerSpill[2].size()="<<tempHitsPerSpill[2].size()<<endl;
  cerr<<"tempHitsPerSpill[3].size()="<<tempHitsPerSpill[3].size()<<endl;
  
  vector<pair<double,vector<TempHit*> > >temp1; 
  vector<pair<double,vector<TempHit*> > >temp2; 

  for(unsigned int i = 0; i< tempHitsPerSpill[0].size(); i++)
    {
      for(unsigned int j =0; j<tempHitsPerSpill[1].size(); j++)
	{
	  if(tempHitsPerSpill[0][i].first == tempHitsPerSpill[1][j].first)
	    {
	      vector<TempHit*> tempVector;
	      tempVector.insert(tempVector.end(),tempHitsPerSpill[0][i].second.begin(),tempHitsPerSpill[0][i].second.end());
	      tempVector.insert(tempVector.end(),tempHitsPerSpill[1][j].second.begin(),tempHitsPerSpill[1][j].second.end());
	      temp1.push_back(make_pair(tempHitsPerSpill[0][i].first,tempVector));
	      //cerr<<"match="<<tempHitsPerSpill[0][i].first<<endl;
	    }
	}
    }

  for(unsigned int i = 0; i< tempHitsPerSpill[2].size(); i++)
    {
      for(unsigned int j =0; j<tempHitsPerSpill[3].size(); j++)
	{
	  if(tempHitsPerSpill[2][i].first == tempHitsPerSpill[3][j].first)
	    {
	      vector<TempHit*> tempVector;
	      tempVector.insert(tempVector.end(),tempHitsPerSpill[2][i].second.begin(),tempHitsPerSpill[2][i].second.end());
	      tempVector.insert(tempVector.end(),tempHitsPerSpill[3][j].second.begin(),tempHitsPerSpill[3][j].second.end());
	      temp2.push_back(make_pair(tempHitsPerSpill[2][i].first,tempVector));
	      //cerr<<"match="<<tempHitsPerSpill[2][i].first<<endl;
	    }
	}
    }
  cerr<<"temp1.size()="<<temp1.size()<<endl;
  cerr<<"temp2.size()="<<temp2.size()<<endl;
  
  for(unsigned int i = 0; i< temp1.size(); i++)
    {
      for(unsigned int j =0; j<temp2.size(); j++)
	{
	  
	  if(temp1[i].first == temp2[j].first)
	    {
	      vector<TempHit*> tempVector;
	      tempVector.insert(tempVector.end(),temp1[i].second.begin(),temp1[i].second.end());
	      tempVector.insert(tempVector.end(),temp2[j].second.begin(),temp2[j].second.end());
	      tempVectorHits->push_back(tempVector);
	      //cerr<<"match="<<tempHitsPerSpill[2][i].first<<endl;
	    }
	  
	}
    }
  
  cerr<<tempVectorHits->size()<<endl;
}


void HandleData2(vector<string> filenames, string filepath,root2dst* cvt,TH1I* h1)
{
  double noCut = 0;
  double firstCut = 0;
  double secondCut = 0;
  double thirdCut = 0;
  double forthCut = 0;
  double fifthCut = 0;
  double sixthCut = 0;
  double seventhCut = 0;
  double eightCut = 0;
  double ninthCut = 0;


  vector<MDdateFile*> dfiles;
  vector<char*> eventBuffers;
  int evt_read = 0;

  TFile rfile("new.root", "recreate");
  TTree* statTree = new TTree("tree", "Tree");
  int _evNo;

  std::vector<double> _evNoVec;

  std::vector<vector <double> > _R;
  std::vector<vector <double> > _Y;
  std::vector< vector<double> > _XPos;
  std::vector< vector<double> > _YPos;
  std::vector< vector<double> > _ZPos;
  std::vector< vector<double> > _HTimeL;
  std::vector< vector<double> > _HTimeR;
  std::vector< vector<double> > _HTime;


  std::vector<double> _RTemp;
  std::vector<double> _YTemp;

  std::vector<double> _XPosTemp;
  std::vector<double> _YPosTemp;
  std::vector<double> _ZPosTemp;
  std::vector<double> _HTimeLTemp;
  std::vector<double> _HTimeRTemp;
  std::vector<double> _HTimeTemp;

  statTree->Branch("MC_Evt", &_evNo, "MC_EventNo/I");
  statTree->Branch("EvNoVec", &_evNoVec, 32000,0);
  statTree->Branch("R", &_R,32000,0);
  statTree->Branch("Y", &_Y,32000,0);
  statTree->Branch("XPos", &_XPos,32000,0);
  statTree->Branch("YPos", &_YPos, 32000,0);
  statTree->Branch("ZPos", &_ZPos,32000,0);
  statTree->Branch("TimeL", &_HTimeL,32000,0);
  statTree->Branch("TimeR", &_HTimeR,32000,0);
  statTree->Branch("Time", &_HTime,32000,0);


  for(unsigned int fileInt=0; fileInt<filenames.size();fileInt++)
    {
      dfiles.push_back(new MDdateFile(filenames[fileInt],filepath));
      eventBuffers.push_back(new char());
    }
  //dfiles.push_back(new MDdateFile(filenames[0],filepath));

  // Open the file and loop over event->.
  //char *eventBuffer;


  bool validFiles = true;
  for(unsigned int fileInt=0;fileInt<dfiles.size();fileInt++)
    {
      validFiles &= dfiles[fileInt]->open();
    }
  if(validFiles){
  //if ( dfiles[0]->open() ) { // There is a valid files to unpack
    //dfiles[0]->init();
    for(unsigned int fileInt=0;fileInt<dfiles.size();fileInt++)
      {
	dfiles[fileInt]->init();
      }
    int counter =0;
    int xEv(0);
    do { // Loop over all spills
      vector<vector<TempHit*> > tempVectorHits;
      for(unsigned int i=0;i<filenames.size();i++)
	{
	  eventBuffers[i] = dfiles[i]->GetNextEvent();
	}


      //eventBuffer =  dfiles[0]->GetNextEvent();
      //if(counter == 2) break; //beam test 2017
      counter++;
      cout<<"dfile.NSpills()="<<dfiles[0]->NSpills()<<endl;
      //if(eventBuffer== NULL) break;
    
      //convert_to_internal(eventBuffers,&tempVectorHits,h1,&noCut,&firstCut,&secondCut,&thirdCut,&forthCut,&fifthCut);
      
      convert_to_internal2017BeamTest(eventBuffers,&tempVectorHits,h1);

    ++xEv;
    //       } while (xEv < 5);
    //} while ( eventBuffers[0] );

    //reader_root inDst;
    int iEvent;
    //int evt_read = 0;
    
    //inDst.open( input_data[0] );
    iEvent = 0;

    _evNoVec.clear();
    _R.clear();
    _Y.clear();
    _XPos.clear();
    _YPos.clear();
    _ZPos.clear();
    _HTimeL.clear();
    _HTimeR.clear();
    _HTime.clear();

    for(unsigned int event=0;event<tempVectorHits.size();event++)
      {
	_evNo = event;
	sixthCut++;
	if(tempVectorHits[event].size()<12) continue;
	seventhCut++;
	cout<<"eventNum="<<event<<endl;
	string filepath2="/data/neutrino05/phallsjo/SaRoMan";
	//vector<bhep::hit*> hits = ConvertToHit(tempVectorHits[event],filepath2+"/test2.xml");
	//vector<bhep::hit*> unSortedhits = ConvertToHit(tempVectorHits[event],filepath2+"/test2.xml");
	dataBaseVector.clear();
	vector<bhep::hit*> unSortedhits = ConvertToHit2017BeamTest(tempVectorHits[event],filepath+"/BabyMIND.db");

	std::vector<bhep::hit*> hits = unSortedhits;
	sort( hits.begin(), hits.end(), forwardSortInZ() );
	
	ptype pT = DIGI;
	string detect = "tracking";
	vector<particle*> hitsParticle;
	hitsParticle.push_back(new particle(pT,detect));

	_evNoVec.push_back(event);

	_RTemp.clear();
	_YTemp.clear();
	_XPosTemp.clear();
	_YPosTemp.clear();
	_ZPosTemp.clear();
	_HTimeLTemp.clear();
	_HTimeRTemp.clear();
	_HTimeTemp.clear();

	double initTime = hits[0]->ddata("time");

	double initTimeL = 0;
	double initTimeR = 0;

	if( hits[0]->idata("rHit"))
	  {
	    initTimeR = hits[0]->ddata("time");
	    initTimeL = hits[1]->ddata("time");
	  }
	else
	  {
	    initTimeL = hits[0]->ddata("time");
	    initTimeR = hits[1]->ddata("time");
	  }
	  

	for(unsigned int i = 0; i<hits.size();i++)
	  {
	    //hitsParticle.back()->add_hit(detect,hits[i]);

	    eightCut++;

	    if(fabs(hits[i]->ddata("time") - initTime) > 10) continue;//filter
	      
	    ninthCut++;


	    _ZPosTemp.push_back(hits[i]->ddata("barPosZ"));
	    if(hits[i]->idata("IsYBar"))
	      {
		_YPosTemp.push_back(hits[i]->ddata("barPosT"));
	      }
	    else
	      {
		_XPosTemp.push_back(hits[i]->ddata("barPosT"));
	      }

	    _YTemp.push_back(hits[i]->idata("IsYBar"));
	    _RTemp.push_back(hits[i]->idata("rHit"));
	      /*
	    if(hits[i]->idata("rHit"))
	      {
		_HTimeRTemp.push_back(hits[i]->ddata("time")-initTimeR);
	      }
	    else
	      {
		_HTimeLTemp.push_back(hits[i]->ddata("time")-initTimeL);
	      }
	      */
	      /*
	    if(fabs(hits[i]->ddata("time") - initTime) > 3000)
	      {
		cerr<<hits[i]->ddata("time")<<endl;
		cerr<<initTime<<endl;
	      }
	      */
	    _HTimeTemp.push_back(hits[i]->ddata("time") - initTime);


	    hitsParticle.back()->add_hit(detect,hits[i]); //new

	  }
	_R.push_back(_RTemp);
	_Y.push_back(_YTemp);
	
	_XPos.push_back(_XPosTemp);
	_YPos.push_back(_YPosTemp);
	_ZPos.push_back(_ZPosTemp);
	_HTimeL.push_back(_HTimeLTemp);
	_HTimeR.push_back(_HTimeRTemp);
	_HTime.push_back(_HTimeTemp);
	
	/*
	  int _evNo;
	  std::vector< vector<double> > _XPos
	  std::vector< vector<double> > _YPos;
	  std::vector< vector<double> > _ZPos;
	  std::vector< vector<double> > _HTimeLHTimeL;
	  std::vector< vector<double> > _HTimeR;
	*/
	
	cout<<"particles size="<<hitsParticle.back()->hits(detect).size()<<endl;
	
	bhep::event e(evt_read);
	e.add_property( "IntType", "CCQE" );
	e.add_property("G4EventID",evt_read); 
	
	particle* digi_part = cvt->create_digital_representation( hitsParticle );
	e.add_digi_particle( digi_part );
	
	bhep::bhep_svc::instance()->get_writer_root().write( e, evt_read );//iEvent );
	evt_read++;
	
	//cout<<e<<endl;
	e.clear();
	
	hitsParticle.clear();
	
      }
    statTree->Fill();
    


    }while(!NullCheck(eventBuffers));
  }
  //dfiles[0]->close();
  for(unsigned int fileInt=0;fileInt<dfiles.size();fileInt++)
    {
      dfiles[fileInt]->close();
    }

		   

  //MDdateFile dfile(filenames[0], filepath);

  rfile.Write();
  //rfile.Close();

      cerr<<"CUTS:"<<endl;
      cerr<<noCut<<endl;
      cerr<<firstCut<<endl;
      cerr<<secondCut<<endl;
      cerr<<thirdCut<<endl;
      cerr<<forthCut<<endl;
      cerr<<fifthCut<<endl;
      cerr<<sixthCut<<endl;
      cerr<<seventhCut<<endl;
      cerr<<eightCut<<endl;
      cerr<<ninthCut<<endl;
      cerr<<"CUTS done"<<endl;

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
      //cout<<"tempHits="<<tempHits.size()<<endl;

      //filename2="feb3_std_10gevmuons_extscint16_nogarbage_refdata4.daq";
      //filename2="FEB3_safe_mode_3_HG40_LG55_test1.daq";
      //ForFile(filename,filepath,&tempHits,&tempVectorHits);

      //void ForFilesV(vector<string> filenames, string filepath, vector<TempHit*>* tempHits,vector<vector<TempHit*> >* tempVectorHits)
      /*
      filepath="/data/neutrino05/phallsjo/SaRoMan";
      filename="FEB1_safe_mode_3_HG40_LG55_test1.daq";
      filename2="FEB3_safe_mode_3_HG40_LG55_test1.daq";
      filename3="FEB4_safe_mode_3_HG40_LG55_test1.daq";
      filename4="FEB5_safe_mode_3_HG40_LG55_test1.daq";

      filepath="/data/neutrino05/phallsjo/SaRoMan/1GeV/Hadrons/1GeV_hadrons_OR96t";
      filename="FEB1_OR96t_fpga_allch_adcen_only_holdHG65_holdLG65_hadrons_1GeV_t2.daq";
      filename2="FEB3_OR96t_fpga_allch_adcen_only_holdHG65_holdLG65_hadrons_1GeV_t2.daq";
      filename3="FEB4_OR96t_fpga_allch_adcen_only_holdHG65_holdLG65_hadrons_1GeV_t2.daq";
      filename4="FEB5_OR96t_fpga_allch_adcen_only_holdHG65_holdLG65_hadrons_1GeV_t2.daq";
      */


      //filepath="/data/neutrino05/phallsjo/SaRoMan";
      //filename="FEB1_AllCh_minusCh95_safe_mode_10_muons.daq";
      //filename2="FEB3_AllCh_minusCh95_safe_mode_10_muons.daq";
      //filename3="FEB5_AllCh_minusCh29317895_safe_mode_10_muons.daq";
      //filename4="FEB4_AllCh_safe_mode_10_muons.daq";


      //filename3="feb5_std_10GeVmuons_extscint16_nogarbage_refdata4.daq";
      //filename4="feb4_std_10GeVmuons_extscint16_nogarbage_refdata4.daq";

      //filepath="/data/neutrino05/phallsjo/SaRoMan/Safe_mode_10";

      //filename="FEB1_AllCh_minusCh95_safe_mode_10_muons.daq";
      //filename2="FEB3_AllCh_minusCh95_safe_mode_10_muons.daq";
      //filename3="FEB4_AllCh_safe_mode_10_muons.daq";
      //filename4="FEB5_AllCh_minusCh29317895_safe_mode_10_muons.daq";

      filepath="/data/neutrino05/phallsjo/copy/SaRoMan/daq";
      //filename="5GeVc_muons/ModulePosition7/FEB3/FEB3_modpos7_magnet_ON_beam_ON_muons_5GeVc_t1.daq";
      //filename2="5GeVc_muons/ModulePosition7/FEB4/FEB4_modpos7_magnet_ON_beam_ON_muons_5GeVc_t1.daq";
      //filename3="5GeVc_muons/ModulePosition7/FEB5/FEB5_modpos7_magnet_ON_beam_ON_muons_5GeVc_t1.daq";

      //filename="2GeVc_muons/FEB3/FEB3_mod1_beam_ON_muons_2GeV_t1.daq";
      //filename2="2GeVc_muons/FEB4/FEB4_mod1_beam_ON_muons_2GeV_t1.daq";
      //filename3="2GeVc_muons/FEB5/FEB5_mod1_beam_ON_muons_2GeV_t1.daq";

      //filename="2GeVc_hadrons/FEB3/FEB3_mod1_beam_ON_hadrons_2GeV_t1.daq";
      //filename2="2GeVc_hadrons/FEB4/FEB4_mod1_beam_ON_hadrons_2GeV_t1.daq";
      //filename3="2GeVc_hadrons/FEB5/FEB5_mod1_beam_ON_hadrons_2GeV_t1.daq";

      vector<string> tester;
      string sFileName;
      ifstream fList("/data/neutrino05/phallsjo/copy/SaRoMan/module.list");
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

 cerr<<"start of handledata"<<endl;
 HandleData2(tester,filepath,cvt,&h1);
 cerr<<"End of handledata"<<endl;

 //h1.Write();
 //h2.Write();
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
