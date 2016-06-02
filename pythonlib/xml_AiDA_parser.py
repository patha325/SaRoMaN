#######################################################################################################################
#Created by Patrik Hallsjo @ University of Glasgow
#Need automatic dating through GIT, 
#Modified on 2/6-2016
#Created on 2/6-2016
#######################################################################################################################
#General python import
#######################################################################################################################
from xml.dom import minidom
import sys
#######################################################################################################################
#Importing own python files
#######################################################################################################################

#######################################################################################################################
#Class generation
#######################################################################################################################
class xml_AiDA_parser:
    
    def __init__(self,parsingFile,outFile):
        self.parsingFile = parsingFile
        self.outFile = outFile
        self.xmldoc = minidom.parse(self.parsingFile)

    def Parse_file(self):
        #Variables used
        volumeNodeList = self.xmldoc.getElementsByTagName('volume')
        #print 'volumeNodeList'
        #print volumeNodeList
        #print 'done'
        
        solidRefRefList = []
        volumeRefRefList = []

        boxNodeList = self.xmldoc.getElementsByTagName('box')
        boxXList =[]
        boxYList =[]
        boxZList =[]
        boxXListEval =[]
        boxYListEval =[]
        boxZListEval =[]
        boxDict = {}
        posDict = {}

        
        constantNodeList = self.xmldoc.getElementsByTagName('constant')
        constantDict= {}

        #Get constants
        for constantNode in constantNodeList:
            try:
                constantDict[str(constantNode.attributes['name'].value)] = float(constantNode.attributes['value'].value)
            except ValueError: #If read element is not well defined. Could intead eval the full dict.
                constantDict[str(constantNode.attributes['name'].value)] = eval(str(constantNode.attributes['value'].value),constantDict)

                continue

        for volumeNode in volumeNodeList:
            if(volumeNode.attributes['name'].value == 'EMRDetector'):
                replicaVolList= volumeNode.getElementsByTagName("replicavol")
                auxiliaryList= volumeNode.getElementsByTagName("auxiliary")
                 
            if(volumeNode.attributes['name'].value == 'EMR'):
                posRefList= volumeNode.getElementsByTagName("position")
                 
            

        print replicaVolList[0].attributes['number'].value
        for auxiliary in auxiliaryList:
            if(auxiliary.attributes['auxtype'].value == 'BarLength'):
                print auxiliary.attributes['auxvalue'].value
        print posRefList[0].attributes['z'].value
                      


                
        outfile = open(self.outFile,'w+')
        outfile.write("SOLID\t%s\t%s\t%s\t%s\n"%('S',1000,1000,195*2))
        outfile.write("POS\t%s\t%s\t%s\t%s\n"%('S',0,0,0))
        outfile.close()

#######################################################################################################################
