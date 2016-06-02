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

using namespace std;

// Union used to decode the FEB communication data fields
typedef union
{
    unsigned int value;       // Packed field value

    struct              // Common field type ID
    {
        unsigned int    : 28;   // padding
        unsigned int ID   : 4;      // Field type ID
    };

    struct              // Spill header
    {
        unsigned int Tag  : 16;   // Spill Tag
        unsigned int    : 5;    // padding
        unsigned int BID  : 7;    // Board ID
    }SpH;
    
    union             // Spill trailer
    {
    struct              // Spill trailer 1 Tag
    {
       unsigned int Tag : 16;   // Spill Tag
    };
    struct              // Spill trailer 1 temperature
    {
       unsigned int Hum : 8;    // Humidity
       unsigned int Temp  : 10;   // Temperature
    };
    struct              // Spill trailer 1 common
    {
       unsigned int     : 20;   // padding
       unsigned int SbT : 1;    // Spill trailer 1 subtype
       unsigned int BID : 7;    // Board ID
    };
    struct              // Spill trailer 2
    {
       unsigned int SpTime: 28;   // Spill time
    };
    }SpT;
    
    struct              // Event header
    {
      unsigned int TrTag  : 28; // Trigger Tag
    }EvH;
    
    union             // Event trailer
    {
    struct              // Event trailer 1
    {
       unsigned int TrTag : 28;   // Trigger Tag
    };
    struct              // Event trailer 2
    {  
       unsigned int TrTime: 20;   // Trigger time
       unsigned int     : 3;    // padding
       unsigned int HitCnt: 5;    // Hit counts within event
    };
    }EvT;
    
    union             // Event data
    {
    struct              // Time
    {
       unsigned int HitTime: 16;    // Hit time
    };
    struct              // Amplitude
    {  
       unsigned int Ampl  : 12;   // Amplitude
       unsigned int AmplID: 4;    // Amplitude ID
    };
    struct              // Common
    {
       unsigned int     : 16;   // padding
       unsigned int HitID : 5;    // Hit ID
       unsigned int ChID  : 7;    // Channel ID
    };
    }EvD;
} CommunicationField_t;       /* 4 bytes size */

vector<bhep::hit*> HandleTestBeamData()//(char *sInputFileName)
{
  vector<bhep::hit*> returnVector;

  char sInputFileName[42] = "A2-1+2-HV135-HG55-THD280-TC1-HOLD6-1.daq";

  FILE *pDAQDataFile = NULL;
  CommunicationField_t FieldTest;

  pDAQDataFile = fopen(sInputFileName, "rb");

  std::vector<CommunicationField_t> fieldVec;
  
  while(fread(&FieldTest.value, sizeof(unsigned int), 1, pDAQDataFile) == 1)
  {
    // Process
     /*
    fprintf(pTXTDataFile,
      "\n%d\t%d\t%d\t%d\t%d\t%d\t",
        FieldTest.ID,
        FieldTest.EvD.ChID,
        FieldTest.EvD.HitID,
        FieldTest.EvD.AmplID,
        FieldTest.EvD.Ampl,
        FieldTest.EvD.HitTime);
        */

      //cout<<FieldTest.EvD.HitTime<<endl;
      //fieldVec.push_back(FieldTest);

    // filtering, simply emulating the c code.
    if(FieldTest.EvD.Ampl < 900) continue;

    if(FieldTest.EvD.AmplID != 3) continue;

    if(FieldTest.EvD.ChID != 66) continue;
        
    returnVector.push_back(new bhep::hit());
    
    //returnVector.back()->add_data("time",FieldTest.EvD.HitTime);
    returnVector.back()->add_property("time",(double)FieldTest.EvD.HitTime);

    //cout<<"FieldTest.EvD.HitID: "<<FieldTest.EvD.HitID<<endl;

    //returnVector.back().add_data(name,value);
    //ddata( "barPosZ" )
    //idata( "IsYBar" )
    //ddata( "EnergyDep" )
    //ddata( "time" );
    //idata( "barNumber" )
    //ddata( "barPosT" )

    // Destinguish between real and simulated, add a data structure TrueData?


    }

   fclose(pDAQDataFile);

  return returnVector;
}

int main(int argc, char* argv[]) {
  
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
  /*
  //If file ending is .root do nothing extra
  if (input_data[0].find(".root")!=npos)
    {
      HandleTestBeamData();
    }
  */

  //vector<bhep::hit*> testVec = HandleTestBeamData();
  /*
  for(int cnt=0;cnt<testVec.size();cnt++)
    {
      //cout<<testVec[cnt]->ddata("time")<<endl;
    }
  */
  //cout<<"size:"<<testVec.size()<<endl;

  //If file ending is something else, parse it into a root file.

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
      
      bhep::bhep_svc::instance()->get_writer_root().write( e, evt_read );//iEvent );
      e.clear();
      
      iEvent++;
      evt_read++;
    }

    inDst.close();
  }

  cvt->print(); // Used to print the debug histograms to a root file.

  WriteUtil::CloseOutputDst();
  //inDst.close();

  return 0;
}
