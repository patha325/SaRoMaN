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

#include "TFile.h"
#include "TH1D.h"

#include <exception>

#include "digi/MDfragmentBM.h"
#include "digi/MDpartEventBM.h"
#include "digi/MDargumentHandler.h"
#include "digi/MDdataFile.h"

using namespace std;

class TempHit
{
public:
  TempHit() {};
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
    }
  else
    {
      cout<<"Node is NULL!"<<endl;
    }

  inHit->add_property("barPosZ",barPosZ);
  inHit->add_property("IsYBar",isYBar);
  inHit->add_property("barNumber",barNumber);
  inHit->add_property("barPosT",barPosT);
  inHit->add_property("moduleNum",moduleNum);
  // can we add ->x()[2] how is it done? bhep add?
  if(isYBar) inHit->set_point(*(new bhep::Point3D(0,barPosT,barPosZ)));
  else inHit->set_point(*(new bhep::Point3D(barPosT,0,barPosZ)));


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

  for(unsigned int cnt = 0; cnt<hitsVector.size(); cnt++)
  {
    returnVector.push_back(new bhep::hit());
    double tempEdep = ConvertToEdep(hitsVector[cnt]->boardId,hitsVector[cnt]->channelId,hitsVector[cnt]->amplitude);
    returnVector.back()->add_property("EnergyDep",tempEdep);
    returnVector.back()->add_property("time",(double)hitsVector[cnt]->timeing);
    ConvertToPosition(returnVector.back(),hitsVector[cnt]->boardId,hitsVector[cnt]->channelId,&parse,root_element);
  }

  return returnVector;
}

void ForFile(string filename, string filepath, vector<TempHit*>* tempHits,vector<vector<TempHit*> >* tempVectorHits)
{
    MDdateFile dfile(filename,filepath);

      char *eventBuffer;
      if ( dfile.open() ) { // There is a valid files to unpack
	dfile.init();
	
	int xEv(0);

	int spillCnt = 1;

	do { // Loop over all spills
	  eventBuffer =  dfile.GetNextEvent();
	  try {
	    MDfragmentBM   spill;
	    spill.SetDataPtr(eventBuffer);
	    
	    if(spillCnt ==2) break;
	    
	    spillCnt++;
	    // For each spill and all files, fill vector of hits.	 
	    
	    MDpartEventBM event;
	    int nTr = spill.GetNumOfTriggers();
	    for (int i=0; i<6000;++i){//nTr; ++i) {
	      vector<TempHit*> tempVec;
	      event = spill.GetTriggerEventPtr(i);
	      //           event->Dump();
	      for (int ich=0; ich<BM_FEB_NCHANNELS; ++ich) {
		int nHits = event.GetNLeadingEdgeHits(ich);
		//if (nHits)
		//{
		    //h1.Fill(ich, nHits);
		    //int q = event.GetHitAmplitude(ich, 'h');
		    //h2.Fill(q);
		//}
		
		for(int ih=0; ih<nHits; ih++)
		  {
		    if(event.GetHitTime(ih,ich, 'l') >0 && event.GetHitAmplitude(ich, 'h') > 1000)
		      {
			tempHits->push_back(new TempHit(event.GetHitAmplitude(ich, 'h'),
							event.GetHitTime(ih,ich, 'l'),ich,spill.GetBoardId()));
			//event.GetHitTime(ih,ich, 'l')+4000*i,ich,spill.GetBoardId()));
			//cout<<"spill.GetBoardId="<<spill.GetBoardId()<<endl;

			tempVec.push_back(new TempHit(event.GetHitAmplitude(ich, 'h'),
							event.GetHitTime(ih,ich, 'l'),ich,spill.GetBoardId()));
		      }
		  }
		
	      }
	      tempVectorHits->push_back(tempVec);
	    }
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
	  
	  ++xEv;
	  //       } while (xEv < 5);
	} while ( eventBuffer );
      }
      
      dfile.close();
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

  vector<TempHit*> tempHits; // Vector for spill

  vector<vector<TempHit*> > tempVectorHits; // Vector for spills per trigger.


  if(testBeam)
    {

      // Want to read n number files at a time
      // From, for each spill
      // Add hits togeather using the propper time. Trigg+hit.
      // Or should it be per trigger? Since simulation per particle.

      // How do the event builder?
      char *dataBuff;
      uint32_t* dataPtr;

      string filepath;
      string filename;
      filepath="/data/neutrino05/phallsjo/SaRoMan";
      filename="feb1_std_10gevmuons_extscint16_nogarbage_refdata4.daq";
      //filename="A2-all-10KHz-EXTPULSE_10-000us-ANALOG-1.daq";
      //filename="test.daq";
      
      //TFile rfile("histos.root", "recreate");
      //TH1I  h1("h1_hit_ch", "hit channels", 100, 0, 100);
      //TH1I  h2("h2_ampl", "hit ampl.", 200, 0, 5000);
      
      ForFile(filename,filepath,&tempHits,&tempVectorHits);
      cout<<"tempHits="<<tempHits.size()<<endl;

      filename="feb3_std_10gevmuons_extscint16_nogarbage_refdata4.daq";
      ForFile(filename,filepath,&tempHits,&tempVectorHits);
      
      //h1.Write();
      //h2.Write();
      //rfile.Close();
      //dfile.close();
      //delete dataBuff;

      cout<<"tempHits="<<tempHits.size()<<endl;
      cout<<tempVectorHits.size()<<endl;
      for(unsigned int i=0;i<tempVectorHits.size();i++)
	{
	  if(tempVectorHits[i].size()>0)
	    {
	      cout<<tempVectorHits[i].size()<<endl;
	    }
	}

//vector<bhep::hit*>ConvertToHit(vector<TempHit> hitsVector, string xmlName)
      vector<bhep::hit*> hits = ConvertToHit(tempHits,filepath+"/test2.xml");


      ptype pT = DIGI;
      string detect = "tracking";
      vector<particle*> hitsParticle;
      hitsParticle.push_back(new particle(pT,detect));
      for(unsigned int i = 0; i<hits.size();i++)
	{
	  hitsParticle.back()->add_hit(detect,hits[i]);
	}

      inDst.open( input_data[0] );
      iEvent = 0;

      bhep::event& e = inDst.read_event( iEvent );


      //bhep::event e;
      particle* digi_part = cvt->create_digital_representation( hitsParticle );
      e.add_digi_particle( digi_part );
      bhep::bhep_svc::instance()->get_writer_root().write( e, evt_read );//iEvent );

      // Bhep hits into vector of particles. 
      //particle(ptype type, string name);
      //particle()
      //void add_hit(string detector, hit* ht){
      //ptype pT = DIGI;
      //string detect = "tracking";

      cout<<e<<endl;

	
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

	  cout<<e<<endl;
	  
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
