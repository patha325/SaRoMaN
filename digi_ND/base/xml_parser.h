#ifndef _XML_PARSER__
#define _XML_PARSER__

#include <libxml/parser.h>
#include <libxml/tree.h>
#include <iostream>
#include <string>
#include <cstring>

using namespace std;

class Xml_parser
{
public:
  Xml_parser(string infile);
  ~Xml_parser();
  void print_element_names(xmlNode * a_node);
  xmlNode* ParseXML(void);

  xmlNode* GetNode(xmlNode* node,string inName1, int inValue1, string inName2, int inValue2); 

 private:
  string _infile;
  xmlDoc* _doc;
};


#endif
