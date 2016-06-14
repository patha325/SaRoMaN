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


using namespace std;

class SpillHeader
{
public:
  SpillHeader() {};
  SpillHeader(unsigned long int in)
  {
    boardId = (in >> 20) & 0x7F;
    spillTag = in & 0xFFFF;
      
  };
  unsigned int boardId;
  unsigned int spillTag;
};

class GtrigHeader
{
public:
  GtrigHeader() {};
  GtrigHeader(unsigned long int in)
  {
    triggerTag = in & 0xFFFFFFF; //27-0

  };
  unsigned long int triggerTag;
};

class HitAmplitude
{
public:
  HitAmplitude() {};
  HitAmplitude(unsigned long int in)
  {
    channelId = (in >> 21)& 0x7F;
    r = (in >> 20)& 0x1;
    hitId = (in >> 16)& 0xF;
    amplitudeId = (in >> 12)& 0xF;
    amplitudeMeasurement = in & 0xFFF;
  };
  unsigned int channelId;
  unsigned int r;
  unsigned int hitId;
  unsigned int amplitudeId;
  unsigned int amplitudeMeasurement;
};

class HitTime
{
public:
  HitTime() {};
  HitTime(unsigned long int in)
  {
    channelId = (in >> 21)& 0x7F;
    edge = (in >> 20)& 0xF;
    hitId = (in >> 16)& 0xF;
    hitTime  = in & 0xFFFF;
  };
  unsigned int channelId;
  unsigned int edge;
  unsigned int hitId;
  unsigned int hitTime;

  vector<HitAmplitude*> hitAmplitudeVec;

};

class GtrigTrailer
{
public:
  GtrigTrailer() {};
  GtrigTrailer(vector<unsigned long int> inV)
  {
    triggerTag = inV[0] & 0xFFFFFFF;
    hitInEvt = (inV[1] >> 23)& 0x1F;
    id = (inV[1] >> 19)& 0x3;
    triggerTime = inV[1] & 0x3FFFF;

  };
  unsigned long int triggerTag;
  unsigned int hitInEvt;
  unsigned int id;
  unsigned int triggerTime;
};

class Gtrig
{
public:
  Gtrig() {};

  GtrigHeader* gtrigHeader;

  // channelId as key.
  map<unsigned int,HitTime*> hitTimeMap;
  //mymap.insert ( std::pair<unsigned int, HitTime>('a',100) );

  //vector<HitTime> hitTimeVec;
  GtrigTrailer* gtrigTrailer;

};

class SpillTrailer
{
public:
  SpillTrailer() {};
  SpillTrailer(vector<unsigned long int> inV)
  {
    boardId = (inV[0] << 21)& 0x7F;
    spillTag = inV[0] & 0xFFFF;
    temperature = (inV[1] >> 8)& 0x3FF;
    humidity = inV[1] & 0xFF;
    spillTime = inV[2] & 0xFFFFFFF;
  };
  unsigned int boardId;
  unsigned int spillTag;
  unsigned int temperature;
  unsigned int humidity;
  unsigned long int spillTime;
};

class HitStructure
{
public:
  HitStructure() {};
  // want to take in a vector to produce constructor?
  SpillHeader* header;
  vector<Gtrig*> gtrigVec;
  SpillTrailer* trailer;

  void clear() {gtrigVec.clear();};

};

vector<bhep::hit*> ConvertToHit(HitStructure inStructure)
{

  vector<bhep::hit*> returnVector;

  /*
  for(unsigned int cnt = 0;cnt<inStructure.gtrigVec.size();cnt++
    {
      // Filter depending on channel
      if(inStructure.gtrigVec[cnt]->hitTimeMap.begin()->first == 0)
	{

	}
      
      // Filter depending on hittime
      if(inStructure.gtrigVec[cnt]->hitTimeMap.begin()->second->hitTime == 0)
	{

	}



    }
  */


  /*

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
  
  */
  return returnVector;
}




vector<bhep::hit*> HandleTestBeamDataNew()//(char *sInputFileName)
{
  Xml_parser parse("testXML.xml");
  xmlNode *root_element = parse.ParseXML();
  //parse.print_element_names(root_element);

  //cout<<"Print done"<<endl;
  
  xmlNode *cur_node = parse.GetNode(root_element,"boardID",0,"channelID",1);
  
  if(cur_node != NULL)
    {
      xmlAttr* attribute = cur_node->properties;
      
      while(attribute)
	{
	  printf("node type: Element, properties: %s\n", attribute->name);
	  printf("node type: Element, content: %s\n", xmlNodeGetContent(attribute->children));
	  attribute = attribute->next;
	}
    }
  else
    {
      cout<<"NULL"<<endl;
    }
  
  vector<bhep::hit*> returnVector;
  
  //unsigned char *data=new unsigned char [20];

  //string inFile = "A2-1+2-HV135-HG55-THD280-TC1-HOLD6-1.daq";
  //std::ifstream inFile("A2-1+2-HV135-HG55-THD280-TC1-HOLD6-1.daq",ios::binary | ios::in);
  //std::ifstream inFile("A2-1+2-HV135-HG55-THD280-TC1-HOLD6-1.daq",ios::binary |ios::in);

  std::ifstream inFile("A2-all-10KHz-EXTPULSE_10-000us-ANALOG-1.daq",std::ios::binary);

  //std::ifstream inFile("A2-all-10KHz-EXTPULSE_10-000us-ANALOG-1.daq",ios::binary |ios::in);
 

  // New data format

  // what is a line? Perhaps better to read in char by char.
  // if first 4 bits = spill header, create object
  // keep on reading in lines or chars and filling object with them depending on id.
  // stop when first 4 bits equal spill foter.

  // Each "spill object" will contain x number of triggers with y number of hit times and z number of hit amplitudes.

  // Order in excel, xyyyyyyyzzzzzzzzzzz. When adding hit amplitude z find hit time to add it to.

  vector<HitStructure> spills;
  HitStructure temp;
  vector< unsigned long int> tempV;
  bool start = false;

  while (!inFile.eof()) {
    unsigned long int longWord;
    unsigned char a,b,c,d;
    unsigned char hex1,hex2,hex3,hex4,hex5,hex6,hex7,hex8;
    unsigned char id;
    //inFile >> a >> b >> c >> d;
    // Bug in hardware, backwards ordering.
    
    cin.flags(ios_base::showbase| ios_base::hex);
    inFile >> std::noskipws;
    inFile >> a >> b >> c >> d;
    /*
    hex1 = (a >> 4) & 0xF;
    hex2 = a & 0xF;
    hex3 = (b >> 4) & 0xF;
    hex4 = b & 0xF;
    hex5 = (c >> 4) & 0xF;
    hex6 = c & 0xF;
    hex7 = (d >> 4) & 0xF;
    hex8 = d & 0xF;

    cout<<"hex1="<<hex1<<" "<<std::bitset<4>(hex1)<<endl;
    cout<<"hex2="<<hex2<<" "<<std::bitset<4>(hex2)<<endl;
    cout<<"hex3="<<hex3<<" "<<std::bitset<4>(hex3)<<endl;
    cout<<"hex4="<<hex4<<" "<<std::bitset<4>(hex4)<<endl;
    cout<<"hex5="<<hex5<<" "<<std::bitset<4>(hex5)<<endl;
    cout<<"hex6="<<hex6<<" "<<std::bitset<4>(hex6)<<endl;
    cout<<"hex7="<<hex7<<" "<<std::bitset<4>(hex7)<<endl;
    cout<<"hex8="<<hex8<<" "<<std::bitset<4>(hex8)<<endl;
    */

    id = (d >> 4) & 0xF;

    longWord = a << 24 | b << 16 | c << 8 | d << 0;
    //cout<<hex<<longWord<<endl;

    longWord =0;

    longWord =  a << 0 | b << 8 | c << 16 | d << 24;

    //longWord = a << 24 | b << 16 | c << 8 | d << 0;
    // Bug in hardware, backwards ordering.
    //longWord = a << 0 | b << 8 | c << 16 | d << 24;

    //cout<<"a="<<a<<endl;
    //cout<<std::bitset<8>(a)<<endl;
    /*
    cout<<"a: "<<hex<<a<<endl;
    cout<<"b: "<<hex<<b<<endl;
    cout<<"c: "<<hex<<c<<endl;
    cout<<"d: "<<hex<<d<<endl;
    cout<<hex<<longWord<<endl;
    */
    //cout<<std::bitset<32>(longWord)<<endl;
    //cout<<std::bitset<4>(id)<<endl;
    //cout<<hex<<id<<endl;

    if(id!=0 && !start) continue;

    if(id == 0 && !start)
      {
	//spill header
	temp.header = new SpillHeader(longWord);
	temp.gtrigVec.push_back(new Gtrig());
	start = true;
      }
    else if(id == 1 && start)
      {
	//Gtrig header
	temp.gtrigVec.back()->gtrigHeader = new GtrigHeader(longWord);
      }
    else if(id == 2 && start)
      {
	//Hit time
	// Key using channelID
	temp.gtrigVec.back()->hitTimeMap.insert(std::pair<unsigned int, HitTime*>(
										 ((unsigned int) longWord >> 21) & 0x7F,
										 new HitTime(longWord)));
      }
    else if(id == 3 && start)
      {
	//Hit amplitude
	// Find using channelID
	std::map<unsigned int,HitTime*>::iterator it;


	it = temp.gtrigVec.back()->hitTimeMap.find((longWord >> 21) & 0x7F);
	if (it !=temp.gtrigVec.back()->hitTimeMap.end())
	  {
	    temp.gtrigVec.back()->hitTimeMap.find(
						  (longWord >> 21) & 0x7F
						  )->second->hitAmplitudeVec.push_back(new HitAmplitude(longWord));
	  }
	else
	  {
	    cout<< "Hit amplitude without hit time?!"<<endl;
	    cout<<"ID="<<(unsigned int)((longWord >> 21) & 0x7F)<<endl;
	  }
      }
    else if(id == 4 && start)
      {
	//gtrig trailer 1
	tempV.push_back(longWord);
      }
    else if(id == 5 && start)
      {
	//gtrig trailer 2
	tempV.push_back(longWord);

	if(tempV.size() != 2)
	  {
	    cout<< "gtrig 2 before 1 bad data!"<<endl;
	    break;
	  }

	temp.gtrigVec.back()->gtrigTrailer = new GtrigTrailer(tempV);
	tempV.clear();
      }
    else if(id == 6 && start)
      {
	//spill trailer 1
	tempV.push_back(longWord);
      }
    else if(id == 7 && start)
      {
	//spill trailer 2
	tempV.push_back(longWord);
	
	if(tempV.size() != 3)
	  {
	    cout<< "spill trailer 2 before 2 spill trailer 1 bad data!"<<endl;
	    break;
	  }
	temp.trailer = new SpillTrailer(tempV);
	tempV.clear();
	spills.push_back(temp);
	start = false;
	cout<<"done"<<endl;
      }
    else
      {
	cout<< "wrong id bad data!"<<endl;
	break;
      }


    /*
    std::cout <<std::hex
	      <<static_cast<unsigned>(a)
	      <<static_cast<unsigned>(b)
	      <<static_cast<unsigned>(c)
	      <<static_cast<unsigned>(d)
	      <<endl;
    */
  }
  cout<<"Created HitStructures"<<endl;
  cout<<"spills.size()="<<spills.size()<<endl;
  for(unsigned int cnt =0;cnt<spills.size();cnt++)
    {
      cout<<"spill "<<cnt<<endl;
      for(unsigned int inCnt=0; inCnt<spills[cnt].gtrigVec.size();inCnt++)
	{
	  cout<<"gtrigVec "<<inCnt<<endl;
	  for(map<unsigned int,HitTime*>::iterator it = spills[cnt].gtrigVec[inCnt]->hitTimeMap.begin();
	      it != spills[cnt].gtrigVec[inCnt]->hitTimeMap.end(); ++it)
	    {
	      cout<<"Key="<<it->first<<endl;
	      cout<<"hitTime="<<it->second->hitTime<<endl;
	      cout<<"channelId="<<it->second->channelId<<endl;
	    }
	  cout<<"gtrigVec done"<<endl;
	}
      cout<<"spill done"<<endl;
    }
  cout<<"loop done"<<endl;



  // Old data format
  /*
  std::string line;
  while (std::getline(inFile, line))
    {
    //cout<<line<<endl;
      //cout<<line.length()<<endl;

      cout<<std::hex<<line<<endl;


      //cout<<line.length()<<endl;
      
      std::vector<char> bytes(line.begin(), line.end());

      std::vector<unsigned long int> longV;

      for(int i = 0; i<bytes.size();i+=4)
	{

	  longV.push_back(bytes[i] << 24 | bytes[i+1] << 16 |  bytes[i+2] << 8 | bytes[i+3] << 0);

	  //cout<<longV.back()<<endl;
	  //cout<<std::bitset<32>(longV.back())<<endl;
	}
  
      HitStructure test(longV);

  }
  */
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
