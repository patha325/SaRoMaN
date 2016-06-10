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

#include <locale>

#include <libxml/parser.h>
#include <libxml/tree.h>

using namespace std;

class SpillHeader
{
public:
  SpillHeader() {};
  SpillHeader(unsigned long int in)
  {
    id = (in >> 28)& 0xF;
    boardId = (in >> 20) & 0x7F;
    spillTag = in & 0xFFFF;
      
  };
  unsigned int id;
  unsigned int boardId;
  unsigned int spillTag;
};

class EventHeader
{
public:
  EventHeader() {};
  EventHeader(unsigned long int in)
  {
    id = (in >> 28)& 0xF;
    triggerTag = in & 0xFFFFFFF; //27-0

  };
  unsigned int id;
  unsigned long int triggerTag;
};

class HitAmplitude
{
public:
  HitAmplitude() {};
  HitAmplitude(unsigned long int in)
  {
    hitAmplitudeId = (in >> 28)& 0xF;
    channelId = (in >> 21)& 0x7F;
    hitId = (in >> 16)& 0x1F;
    amplitudeId = (in >> 12)& 0xF;
    amplitudeMeasurement = in & 0xFFF;
  };
  unsigned int hitAmplitudeId;
  unsigned int channelId;
  unsigned int hitId;
  unsigned int amplitudeId;
  unsigned int amplitudeMeasurement;
};

class HitTime
{
public:
  HitTime() {};
  HitTime(vector<unsigned long int> inV)
  {
    hitTimeId = (inV[0] >> 28)& 0xF;
    channelId = (inV[0] >> 21)& 0x7F;
    hitId = (inV[0] >> 16)& 0x1F;
    hitTime  = inV[0]& 0xFFFF;

    if(inV.size() <2)
      {
	cout<<"No hits!!!"<<endl;
      }
    else
      {
	vector<unsigned long int> ampV(inV.begin()+1,inV.end());
	for(unsigned int cnt=0;cnt<ampV.size();cnt++)
	  {
	    hitAmplitudeVec.push_back(HitAmplitude(ampV[cnt]));
	  }
      }
  };
  unsigned int hitTimeId;
  unsigned int channelId;
  unsigned int hitId;
  unsigned int hitTime;

  vector<HitAmplitude> hitAmplitudeVec;

};

class EventTrailer
{
public:
  EventTrailer() {};
  EventTrailer(vector<unsigned long int> inV)
  {
    //unsigned int temp[] = {(inV[0] << 28)& 0xF,(inV[1] << 28)& 0xF};
    //eventTrailerId = temp;
    //memcpy(&eventTrailerId, &temp, sizeof eventTrailerId);

    //eventTrailerId = {(inV[0] << 28)& 0xF,(inV[1] << 28)& 0xF};
    triggerTag = inV[0] & 0xFFFFFFF;
    hitInEvt = (inV[1] >> 23)& 0x1F;
    id = (inV[1] >> 19)& 0x3;
    triggerTime = inV[1] & 0x3FFFF;

  };
  //unsigned int eventTrailerId[2];
  unsigned long int triggerTag;
  unsigned int hitInEvt;
  unsigned int id;
  unsigned int triggerTime;
};

class SpillTrailer
{
public:
  SpillTrailer() {};
  SpillTrailer(vector<unsigned long int> inV)
  {
    //unsigned int temp[] = {(inV[0] << 28)& 0xF,(inV[1] << 28)& 0xF,(inV[2] << 28)& 0xF};
    //unsigned int temp2[] = {(inV[0] << 21)& 0x7F,(inV[1] << 21)& 0x7F};

    //spillTrailerId = temp;
    //boardId =temp2;

    //memcpy(&spillTrailerId, &temp, sizeof spillTrailerId);
    //memcpy(&boardId, &temp2, sizeof boardId);

  //spillTrailerId = {(inV[0] << 28)& 0xF,(inV[1] << 28)& 0xF,(inV[2] << 28)& 0xF};
  //boardId = {(inV[0] << 21)& 0x7F,(inV[1] << 21)& 0x7F};

    spillTag = inV[0] & 0xFFFF;
    temperature = (inV[1] >> 8)& 0x3FF;
    humidity = inV[1] & 0xFF;
    spillTime = inV[2] & 0xFFFFFFF;

  };
  //unsigned int spillTrailerId[3];
  //unsigned int boardId[2];
  unsigned int spillTag;
  unsigned int temperature;
  unsigned int humidity;
  unsigned long int spillTime;

};

class HitStructure
{
public:
  HitStructure() {};
  HitStructure(vector<unsigned long int> longV)
  {
    if(longV.size() < 8)
      {
      cout<<"No hit info"<<endl;
      }
    else
      {

    cout<<"Starting constructor"<<endl;
    header = SpillHeader(longV[0]);
    cout<<"Created SpillHeader"<<endl;
    eventHeader = EventHeader(longV[1]);
    cout<<"Created EventHeader"<<endl;
    triggerTag = longV[1] & 0xFFFFFFF; //27-0

    vector<unsigned long int> hitsV(longV.begin()+2,longV.end()-5);
    /*
    std::cout << "HITID = " << longV[2]<< std::endl;
    std::cout << "HITID = " << std::bitset<32>((longV[2])) << std::endl;
    std::cout << "HITID = " << std::bitset<5>((longV[2] >> 16) & 0x1F) << std::endl;

    std::cout << "HITID = " << std::bitset<32>((longV[3])) << std::endl;
    std::cout << "HITID = " << std::bitset<5>((longV[3] >> 16) & 0x1F) << std::endl;
    */
    // Create hitTimeVec.
    // When hit id changes, new hit.

    unsigned int baseHitId = (hitsV[0] >> 16) & 0x1F;
    vector<unsigned long int> temp;

    for(unsigned int cnt = 0; cnt<hitsV.size(); cnt++)
      {
	//std::cout << "LONG = " << std::bitset<32>(hitsV[cnt]) << std::endl;
	//std::cout << "HITID = " << std::bitset<5>((hitsV[cnt] >> 16) & 0x1F) << std::endl;
	if((hitsV[cnt] >> 16) & 0x1F == baseHitId)
	  {
	    temp.push_back(hitsV[cnt]);
	  }
	else
	  {
	    hitTimeVec.push_back(HitTime(temp));
	    temp.clear();
	    // Next hit
	    baseHitId = (hitsV[cnt] >> 16) & 0x1F;
	    temp.push_back(hitsV[cnt]);
	  }

      }

    cout<<"Created HitTimeVec"<<endl;

    vector<unsigned long int> eventV(longV.end()-5,longV.end()-3);

    eventTrailer = EventTrailer(eventV);

    cout<<"Created EventTrailer"<<endl;

    vector<unsigned long int> spillV(longV.end()-3,longV.end());

    spillTrailer = SpillTrailer(spillV);

    cout<<"Created SpillTrailer"<<endl;
    
      }

  };
  // want to take in a vector to produce constructor?
  SpillHeader header;
  EventHeader eventHeader;
  vector<HitTime> hitTimeVec;
  EventTrailer eventTrailer;
  SpillTrailer spillTrailer;

  unsigned long int triggerTag;


};







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


vector<bhep::hit*> ConvertToHit(HitStructure inStructure)
{
  vector<bhep::hit*> returnVector;


// Read in 2 txt files that contains info needed to convert.
// one for edep per channel per board
// one for x/y z ybar, num 


  for(unsigned int cnt=0;cnt<inStructure.hitTimeVec.size();cnt++)
    {
      for(unsigned int inCnt=0;inCnt<inStructure.hitTimeVec[cnt].hitAmplitudeVec.size();inCnt++)
	{
	  //ddata( "EnergyDep" )
	  //inStructure.hitTimeVec[cnt].hitAmplitudeVec[inCnt].amplitudeMeasurement
	  //inStructure.hitTimeVec[cnt].hitAmplitudeVec[inCnt].channelId
	  //inStructure.spillHeader.boardId
	  
	  //ddata( "barPosZ" )
	  //idata( "IsYBar" )
	  //idata( "barNumber" )
	  //ddata( "barPosT" )
	  //inStructure.hitTimeVec[cnt].hitAmplitudeVec[inCnt].channelId
	  //inStructure.spillHeader.boardId

	  //ddata( "time" );
	  //inStructure.hitTimeVec[cnt].hitTime
	}

    }


  

    //returnVector.back()->add_data("time",FieldTest.EvD.HitTime);
    //returnVector.back()->add_property("time",(double)FieldTest.EvD.HitTime);

    //cout<<"FieldTest.EvD.HitID: "<<FieldTest.EvD.HitID<<endl;

    //returnVector.back().add_data(name,value);
  

  return returnVector;
}

static void print_element_names(xmlNode * a_node)
{
    xmlNode *cur_node = NULL;

    for (cur_node = a_node; cur_node; cur_node = cur_node->next) {
        if (cur_node->type == XML_ELEMENT_NODE) {
            printf("node type: Element, name: %s\n", cur_node->name);

	    xmlAttr* attribute = cur_node->properties;

	    while(attribute)
	      {
		printf("node type: Element, properties: %s\n", attribute->name);
		// if the attribute name ==...
		// if that content == .. 
		printf("node type: Element, content: %s\n", xmlNodeGetContent(attribute->children));
		attribute = attribute->next;
	      }


        }

        print_element_names(cur_node->children);
    }
}


void parseXML(void)
{
    xmlDoc *doc = NULL;
    xmlNode *root_element = NULL;


    /*
     * this initialize the library and check potential ABI mismatches
     * between the version it was compiled for and the actual shared
     * library used.
     */
    LIBXML_TEST_VERSION

    /*parse the file and get the DOM */
    doc = xmlReadFile("testXML.xml", NULL, 0);

    if (doc == NULL) {
        printf("error: could not parse file %s\n", "testXML.xml");
    }

    /*Get the root element node */
    root_element = xmlDocGetRootElement(doc);
    cout<<"Printing"<<endl;
    print_element_names(root_element);
    cout<<"End of print"<<endl;

    //


    /*free the document */
    xmlFreeDoc(doc);

    /*
     *Free the global variables that may
     *have been allocated by the parser.
     */
    xmlCleanupParser();







}





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
    cout<<"value"<<FieldTest.value<<endl;

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

vector<bhep::hit*> HandleTestBeamDataNew()//(char *sInputFileName)
{
  parseXML();

  vector<bhep::hit*> returnVector;
  
  char buffer[20];

  //unsigned char *data=new unsigned char [20];


  //string inFile = "A2-1+2-HV135-HG55-THD280-TC1-HOLD6-1.daq";
  
  
  //std::ifstream inFile("A2-1+2-HV135-HG55-THD280-TC1-HOLD6-1.daq",ios::binary | ios::in);
  std::ifstream inFile("A2-1+2-HV135-HG55-THD280-TC1-HOLD6-1.daq",ios::binary |ios::in);

  //std::locale loc("");
  //cout<<std::locale("").name().c_str()<<endl;

  //inFile.imbue(loc);

  //inFile.read(buffer, sizeof(buffer));

  //inFile.read(buffer, sizeof(buffer));

  //inFile.read( (char*)data,20);
  /*
  while (!inFile.eof()) {
    unsigned char a,b,c,d;
    inFile >> a >> b >> c >> d;
    std::cout <<std::hex
	      <<static_cast<unsigned>(a)
	      <<static_cast<unsigned>(b)
	      <<static_cast<unsigned>(c)
	      <<static_cast<unsigned>(d)
	      <<endl;
  }
  */

  
  std::string line;
  while (std::getline(inFile, line))
    {
    //cout<<line<<endl;
      //cout<<line.length()<<endl;

      cout<<std::hex<<line<<endl;

      /*
      for(int i = 0; i<line.size();i+=4)
	{
	  cout<<std::bitset<4>(atoi(line.substr(i,1).c_str()))
	      <<std::bitset<4>(atoi(line.substr(i+1,1).c_str()))
	      <<std::bitset<4>(atoi(line.substr(i+2,1).c_str()))
	      <<std::bitset<4>(atoi(line.substr(i+3,1).c_str()))<<endl;
	}
      */


      //cout<<line.length()<<endl;
      
      std::vector<char> bytes(line.begin(), line.end());

      //std::vector<unsigned short> bytesInt(line.begin(), line.end());
      /*
      cout<<"Begin"<<endl;
      for(int i = 0; i<bytesInt.size();i+=4)
	{
	  cout<<std::hex<<bytes[i]<<endl;
	  cout<<std::hex<<bytes[i+1]<<endl;
	  cout<<std::hex<<bytes[i+2]<<endl;
	  cout<<std::hex<<bytes[i+3]<<endl;
	  cout<<"done"<<endl;
	  
	  cout<<std::bitset<8>(bytes[i])
	      <<std::bitset<8>(bytes[i+1])
	      <<std::bitset<8>(bytes[i+2])
	      <<std::bitset<8>(bytes[i+3])
	      <<endl;
	  
	}
      cout<<"END"<<endl;
      */


      /*
      cout<<"Begin"<<endl;
      for(int i = 0; i<bytesInt.size();i+=4)
	{
	  cout<<std::hex<<bytesInt[i]<<endl;
	  cout<<std::hex<<bytesInt[i+1]<<endl;
	  cout<<std::hex<<bytesInt[i+2]<<endl;
	  cout<<std::hex<<bytesInt[i+3]<<endl;
	  cout<<"done"<<endl;
	  
	  cout<<std::bitset<8>(bytesInt[i])
	      <<std::bitset<8>(bytesInt[i+1])
	      <<std::bitset<8>(bytesInt[i+2])
	      <<std::bitset<8>(bytesInt[i+3])
	      <<endl;
	  
	}
      cout<<"END"<<endl;
      */


      //cout<<"TESTING"<<endl;
      //cout<<bytes[0]<<bytes[1]<<bytes[2]<<bytes[3]<<endl;
      //cout<<"END TESTING"<<endl;

      //cout<<std::bitset<8>(bytes[0])<<std::bitset<8>(bytes[1])<<std::bitset<8>(bytes[2])<<std::bitset<8>(bytes[3])<<endl;
      //cout<<bytes.size()<<endl;
      //cout<<bytes[0]<<endl;
      //std::cout << "LONG = " << std::bitset<8>(bytes[0]) << std::endl;

  
      std::vector<unsigned long int> longV;

      for(int i = 0; i<bytes.size();i+=4)
	{
	  /*
	  cout<<bytes[i]<<bytes[i+1]<<bytes[i+2]<<bytes[i+3]<<endl;
	  cout<<std::bitset<8>(bytes[i])
	      <<std::bitset<8>(bytes[i+1])
	      <<std::bitset<8>(bytes[i+2])
	      <<std::bitset<8>(bytes[i+3])
	      <<endl;
	  */
	  longV.push_back(bytes[i] << 24 | bytes[i+1] << 16 |  bytes[i+2] << 8 | bytes[i+3] << 0);

	  //cout<<longV.back()<<endl;
	  //cout<<std::bitset<32>(longV.back())<<endl;
	}
  
      HitStructure test(longV);

      //cout<<longV.size()<<endl;
      //cout<<longV[0]<<endl;

      //longV.clear();
	
      //cout<<0xFF<<endl;
      //cout<<(0xFF >> 2)<<endl;
      //cout<<((0xFF >> 4) & 0x7) <<endl;

      //std::cout << "LONG = " << std::bitset<32>(longV[0]) << std::endl;
      //std::cout << "LONG = " << std::bitset<4>(longV[0]) << std::endl;
      
  }
  
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
  //(void) HandleTestBeamData();
  (void) HandleTestBeamDataNew();

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
