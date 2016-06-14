#include <xml_parser.h>

Xml_parser::Xml_parser(string infile)
{
  _infile=infile;
  _doc = NULL;
}

Xml_parser::~Xml_parser()
{
  /*free the document */
  xmlFreeDoc(_doc);
  /*
   *Free the global variables that may
   *have been allocated by the parser.
   */
  xmlCleanupParser();
}

void Xml_parser::print_element_names(xmlNode * a_node)
{
    xmlNode *cur_node = NULL;
    char str[] = "boardID";

    for (cur_node = a_node; cur_node; cur_node = cur_node->next) {
        if (cur_node->type == XML_ELEMENT_NODE) {
            printf("node type: Element, name: %s\n", cur_node->name);

	    xmlAttr* attribute = cur_node->properties;

	    while(attribute)
	      {
		if (strncmp ((char*)attribute->name,str,7) == 0)
		  {
		    printf("node type: Element, properties: %s\n", attribute->name);
		    
		    
		    // if the attribute name ==...
		    // if that content == .. 
		    printf("node type: Element, content: %s\n", xmlNodeGetContent(attribute->children));
		    attribute = attribute->next;
		  }
		else
		  {
		    attribute = attribute->next;
		  }
	      }

        }

        print_element_names(cur_node->children);
    }
}


xmlNode* Xml_parser::ParseXML(void)
{
    xmlNode *root_element = NULL;

    /*
     * this initialize the library and check potential ABI mismatches
     * between the version it was compiled for and the actual shared
     * library used.
     */
    LIBXML_TEST_VERSION

    /*parse the file and get the DOM */
      _doc = xmlReadFile(_infile.c_str(), NULL, 0);

    if (_doc == NULL) {
        printf("error: could not parse file %s\n", "testXML.xml");
    }

    /*Get the root element node */
    root_element = xmlDocGetRootElement(_doc);

    return root_element;
}

xmlNode* Xml_parser::GetNode(xmlNode* node,string inName1, int inValue1, string inName2, int inValue2)
{
    xmlNode *cur_node = NULL;
    bool first = false;
    bool second = false;

    //cout<<"test"<<endl;

    for (cur_node = node->children; cur_node; cur_node = cur_node->next) {
      //cout<<"test2"<<endl;
      first = false;
      second = false;
        if (cur_node->type == XML_ELEMENT_NODE) {
	  //cout<<"test3"<<endl;
	  //printf("node type: Element, name: %s\n", cur_node->name);

	    xmlAttr* attribute = cur_node->properties;

	    while(attribute)
	      {
		if(strncmp ((char*)attribute->name,inName1.c_str(),strlen(inName1.c_str())) == 0)
		  {
		    if(atoi((char *)xmlNodeGetContent(attribute->children))==inValue1)
		      {
			first = true;
		      }
		  }
		if(strncmp ((char*)attribute->name,inName2.c_str(),strlen(inName2.c_str())) == 0)
		  {
		    if(atoi((char *)xmlNodeGetContent(attribute->children))==inValue2)
		      {
			second = true;
		      }
		  }
		
		if(first && second)
		  {
		    break;
		  }
		else
		  {
		    attribute = attribute->next;
		  }		  
	      }
	    if(first && second)
	      {
		break;
	      }
        }
    }
    return cur_node;
}
