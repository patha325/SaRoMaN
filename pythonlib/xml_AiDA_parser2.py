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
class xml_AiDA_parser2:
    
    def __init__(self,parsingFile,outFile):
        self.parsingFile = parsingFile
        self.outFile = outFile
        self.xmldoc = minidom.parse(self.parsingFile)

    def Parse_partial_file(self,parsingFileName):
        print parsingFileName

        self.save = self.xmldoc

        self.xmldoc = minidom.parse(parsingFileName)

        self.Parse_file()

        self.xmldoc = self.save

    def Parse_file(self):

        physNodeVolRefList = []
        initZPosList = []
        replicaList = []
        widthList = []
        barLengthList = []
        centerZPosList = []

        
        worldRef = self.xmldoc.getElementsByTagName('setup')[0].getElementsByTagName('world')[0].attributes['ref'].value

        #print worldRef

        volumeNodeList = self.xmldoc.getElementsByTagName('volume')

        #print volumeNodeList

        for volumeNode in volumeNodeList:
            if(volumeNode.attributes['name'].value == worldRef):
                physNodeList=volumeNode.getElementsByTagName('physvol')
                if(len(physNodeList) ==1): #If 1, old gdml format, 2 repeating structure.
                    old = 1
                else :
                    old = 0

        if(old == 0):
            #Find the initial and center z positions from the physical volumes.
            for physNode in physNodeList:
                if(len(physNode.getElementsByTagName('volumeref')) ==0):
                    self.Parse_partial_file(physNode.getElementsByTagName('file')[0].attributes['name'].value)
                    continue
                #print physNode.getElementsByTagName('volumeref')[0].attributes['ref'].value
                physNodeVolRefList.append(physNode.getElementsByTagName('volumeref')[0].attributes['ref'].value)
                #print physNode.getElementsByTagName('volumeref')[0].attributes['ref'].value[-8:]
                if(physNode.getElementsByTagName('volumeref')[0].attributes['ref'].value[-8:] != 'Detector'):
                    initZPosList.append(physNode.getElementsByTagName('position')[0].attributes['z'].value)
                else:
                    centerZPosList.append(physNode.getElementsByTagName('position')[0].attributes['z'].value)


            #Find the number of replica volumes, z width of modules and bar xy lengths from the physical volume names.
            for physNodeVolRef in physNodeVolRefList:
                for volumeNode in volumeNodeList:
                    if(volumeNode.attributes['name'].value == physNodeVolRef):
                        if volumeNode.getElementsByTagName("replicavol"):
                            replicaList.append(volumeNode.getElementsByTagName("replicavol")[0].attributes['number'].value)
                        if volumeNode.getElementsByTagName("width"):
                            widthList.append(volumeNode.getElementsByTagName("width")[0].attributes['value'].value)
                        auxiliaryList= volumeNode.getElementsByTagName("auxiliary")
                        for auxiliary in auxiliaryList:
                            if(auxiliary.attributes['auxtype'].value == 'BarLength'):
                                barLengthList.append(auxiliary.attributes['auxvalue'].value)


            if(len(initZPosList)!=len(replicaList)):
                print 'Can not parse this gdml file'
            else:
                zlength = abs(eval(initZPosList[0]))*2*10
                xlength = eval(barLengthList[0])*10
                ylength = eval(barLengthList[0])*10
                outfile = open(self.outFile,'w+')
                outfile.write("SOLID\t%s\t%s\t%s\t%s\n"%('S',xlength,ylength,zlength))
                outfile.write("POS\t%s\t%s\t%s\t%s\n"%('S',0,0,centerZPosList[0]))
                outfile.close()

        if(old == 1):
            self.Old_parse_file()
                        
    def Old_parse_file(self):
        #Variables used
        volumeNodeList = self.xmldoc.getElementsByTagName('volume')
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

        
        #Get volume ref & position name
        for volumeNode in volumeNodeList:
            if(volumeNode.attributes['name'].value == 'Detector'):
                volumeRefList= volumeNode.getElementsByTagName("volumeref")
                for volumeRef in volumeRefList:
                    volumeRefRefList.append(volumeRef.attributes['ref'].value)
                posNameRefList= volumeNode.getElementsByTagName("position")
                for posNameRef in posNameRefList:
                    name_val = str(posNameRef.attributes['name'].value)
                    x_val = eval(str(posNameRef.attributes['x'].value),constantDict)
                    y_val = eval(str(posNameRef.attributes['y'].value),constantDict)
                    z_val = eval(str(posNameRef.attributes['z'].value),constantDict)
                    posDict[name_val] = ((x_val,y_val,z_val))
                    
        #Get solid ref
        for volumeRefRef in volumeRefRefList: 
            for volumeNode in volumeNodeList:
                if(volumeNode.attributes['name'].value == volumeRefRef):
                    solidRefList= volumeNode.getElementsByTagName("solidref")
                    for solidRef in solidRefList:
                        solidRefRefList.append(solidRef.attributes['ref'].value)

        #Get solid ref xyz
        for solidRefRef in solidRefRefList: 
            for boxNode in boxNodeList:
                if(boxNode.attributes['name'].value == solidRefRef):
                    x_val = eval(str(boxNode.attributes['x'].value),constantDict)
                    y_val = eval(str(boxNode.attributes['y'].value),constantDict)
                    z_val = eval(str(boxNode.attributes['z'].value),constantDict)
                    boxDict[str(solidRefRef)] = ((x_val,y_val,z_val))
                
        outfile = open(self.outFile,'w+')
        for key,value in boxDict.iteritems():
            outfile.write("SOLID\t%s\t%s\t%s\t%s\n"%(key[:-6],value[0],value[1],value[2]))
        for key,value in posDict.iteritems():
            outfile.write("POS\t%s\t%s\t%s\t%s\n"%(key[:-4],value[0],value[1],value[2]))
        outfile.close()                   




    



#######################################################################################################################
