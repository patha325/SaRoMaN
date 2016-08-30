import os

class generate_EMR_GDML:

    def __init__(self):
        self.sensitivedetector = 'EMR'
        self.NbBarsPerPlane = 90
        self.NbPlanes = 45
        self.NbOfBars = self.NbPlanes * self.NbBarsPerPlane 
        self.unit = 'cm'
        self.BarWidth = 1.0
        self.BarHeight = 0.75
        self.BarLength = 90.0
        self.Gap = 0.05
        self.ActiveDepth = self.NbPlanes * (self.BarHeight + self.Gap)
        self.PlaneWidth  = self.NbBarsPerPlane * (self.BarWidth + self.Gap)
        self.HoleRad = 0.3
        self.FiberCladdingExtRadius = 0.12
        self.AddWLSFiber = 1

        self.outfilename = 'half_AiDA.gdml'
    
        
    def print_baseline_EMR(self):
        data = '''<?xml version="1.0" encoding="UTF-8" standalone="no" ?>
<gdml xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="http://service-spi.web.cern.ch/service-spi/app/releases/GDML/schema/gdml.xsd">
    
  <define>
    <rotation name="RotateZ" z="90" unit="deg"/>
  </define>
        

  <materials>
    <isotope N="1" Z="1" name="H10xa133600">
      <atom unit="g/mole" value="1.00782503081372"/>
    </isotope>
    <isotope N="2" Z="1" name="H20xa133d90">
      <atom unit="g/mole" value="2.01410199966617"/>
    </isotope>
    <element name="H0xa133b30">
      <fraction n="0.999885" ref="H10xa133600"/>
      <fraction n="0.000115" ref="H20xa133d90"/>
    </element>  
    <isotope N="12" Z="6" name="C120xa126a70">
      <atom unit="g/mole" value="12"/>
    </isotope>
    <isotope N="13" Z="6" name="C130xa126b20">
      <atom unit="g/mole" value="13.0034"/>
    </isotope>
    <element name="C0xa126810">
      <fraction n="0.9893" ref="C120xa126a70"/>
      <fraction n="0.0107" ref="C130xa126b20"/>
    </element>
    <material name="POLYSTYRENE" state="solid">
      <D unit="g/cm3" value="1.06"/>
      <fraction n="0.077418" ref="H0xa133b30"/>
      <fraction n="0.922582" ref="C0xa126810"/>
    </material>
    <isotope N="14" Z="7" name="N140xa126de0">
      <atom unit="g/mole" value="14.0031"/>
    </isotope>
    <isotope N="15" Z="7" name="N150xa126e80">
      <atom unit="g/mole" value="15.0001"/>
    </isotope>
    <element name="N0xa126b80">
      <fraction n="0.99632" ref="N140xa126de0"/>
            <fraction n="0.00368" ref="N150xa126e80"/>
    </element>
    <isotope N="16" Z="8" name="O160xa1270f0">
      <atom unit="g/mole" value="15.9949"/>
    </isotope>
    <isotope N="17" Z="8" name="O170xa127180">
      <atom unit="g/mole" value="16.9991"/>
    </isotope>
    <isotope N="18" Z="8" name="O180xa1271f0">
      <atom unit="g/mole" value="17.9992"/>
    </isotope>
    <element name="O0xa126ec0">
      <fraction n="0.99757" ref="O160xa1270f0"/>
      <fraction n="0.00038" ref="O170xa127180"/>
      <fraction n="0.00205" ref="O180xa1271f0"/>
    </element>
    <isotope N="36" Z="18" name="Ar360xa1270a0">
      <atom unit="g/mole" value="35.9675"/>
    </isotope>
    <isotope N="38" Z="18" name="Ar380xa1274b0">
      <atom unit="g/mole" value="37.9627"/>
    </isotope>
    <isotope N="40" Z="18" name="Ar400xa1275a0">
      <atom unit="g/mole" value="39.9624"/>
    </isotope>
    <element name="Ar0xa127230">
      <fraction n="0.003365" ref="Ar360xa1270a0"/>
      <fraction n="0.000632" ref="Ar380xa1274b0"/>
      <fraction n="0.996003" ref="Ar400xa1275a0"/>
    </element>
    <material name="AIR" state="gas">
      <D unit="g/cm3" value="0.00120479"/>
      <fraction n="0.000124000124000124" ref="C0xa126810"/>
      <fraction n="0.755267755267755" ref="N0xa126b80"/>
      <fraction n="0.231781231781232" ref="O0xa126ec0"/>
      <fraction n="0.0128270128270128" ref="Ar0xa127230"/>
    </material>       
    <material name="Glue" state="solid">
      <D unit="g/cm3" value="1.36"/>
      <composite n="261" ref="C0xa126810"/>
      <composite n="304" ref="H0xa133b30"/>
      <composite n="40" ref="O0xa126ec0"/>
    </material>
    <material name="Acrylic" state="solid">
      <D value="1.18" unit="g/cm3"/>
      <composite n="5" ref="C0xa126810"/>
      <composite n="8" ref="H0xa133b30"/>
      <composite n="2" ref="O0xa126ec0"/>
    </material>
  </materials>
  
  <solids>
    <box x="108.40" y="108.40" z="%(ActiveDepth)s" lunit="cm" name="AiDA_solid"/>
    <box x="108.40" y="108.40" z="%(ActiveDepth)s" lunit="cm" name="Calorimeter_solid"/>
        <box x="158.40" y="158.40" z="1.75" lunit="cm"
	 name="Plane_solid"/>
    <box x="100" y="100" z="1.5" lunit="cm"
	 name="PlanePair_solid"/>

    <tube rmax="%(HoleRad)s" z="%(BarLength)s" lunit="cm" deltaphi="360" aunit="degree" name="hole_solid"/>
    <tube rmax="%(FiberCladdingExtRadius)s" z="%(BarLength)s" lunit="cm" deltaphi="360" aunit="degree" name="Cladding_solid"/>
    <tube rmax="0.114" z="%(BarLength)s" lunit="cm" deltaphi="360" aunit="degree" name="FiberCore_solid"/>
    
    <box x="%(BarLength)s" y="%(BarWidth)s" z="%(BarHeight)s" lunit="cm" name="activeX_bar_solid"/>
    <box x="%(BarWidth)s" y="%(BarLength)s" z="%(BarHeight)s" lunit="cm" name="activeY_bar_solid"/>
    <box x="%(BarLength)s" y="%(BarLength)s" z="%(BarHeight)s" lunit="cm"
	 name="activeX_plane_solid"/>
    <box x="%(BarLength)s" y="%(BarLength)s" z="%(BarHeight)s" lunit="cm"
	 name="activeY_plane_solid"/>
    <box x="%(PlaneWidth)s" y="%(PlaneWidth)s" z="%(ActiveDepth)s" lunit="cm"
	 name="active_TASD_solid"/> 

  </solids>
  
  <structure>       
    <!--  -->
    <volume name="BarX">
      <materialref ref="POLYSTYRENE"/>
      <solidref ref="activeX_bar_solid"/>
      <auxiliary auxtype="SD" auxvalue="Bar"/>
      <auxiliary auxtype="Color" auxvalue="Blue"/>
      <!-- <auxiliary auxtype="Invisible" auxvalue="1"/> -->
      <auxiliary auxtype="EMField" auxvalue="NULL"/>
    </volume>
    <volume name="BarY">
      <materialref ref="POLYSTYRENE"/>
      <solidref ref="activeY_bar_solid"/>
      <auxiliary auxtype="SD" auxvalue="Bar"/>
      <auxiliary auxtype="Color" auxvalue="Blue"/>
      <!-- <auxiliary auxtype="Invisible" auxvalue="1"/>  -->
      <auxiliary auxtype="EMField" auxvalue="NULL"/>
    </volume>
        '''%vars(self)

        print data
        return data

        
    def generate_parvol(self):
        data = '''
    <volume name="TASD">
      <materialref ref="AIR"/>
      <solidref ref="active_TASD_solid"/>
      <paramvol ncopies="%(NbOfBars)s">
         <volumeref ref="BarX"/>
         <parameterised_position_size>'''%vars(self)

        # Calculate bar definition by copy number
        params = {}
        for i in range(0, self.NbOfBars):
            
            params['copy'] = i + 1
            params['zpos'] = ( int(i / self.NbBarsPerPlane) + 0.5 ) \
                             * (self.BarHeight + self.Gap) - self.ActiveDepth / 2.
            params['z']    = self.BarHeight
            # even planes measure x
            if int(i / self.NbBarsPerPlane) % 2 == 0:
                params['xpos'] = (i%self.NbBarsPerPlane + 0.5 \
                                  - self.NbBarsPerPlane / 2.) *\
                                 (self.BarWidth + self.Gap)
                params['ypos'] = 0
                params['x']   = self.BarWidth
                params['y']   = self.BarLength
            else:
                params['ypos'] = (i%self.NbBarsPerPlane + 0.5 \
                                  - self.NbBarsPerPlane / 2.) *\
                                 (self.BarWidth + self.Gap)
                params['xpos'] = 0
                params['y']   = self.BarWidth
                params['x']   = self.BarLength
                
            

            data += '''
        <parameters number="%(copy)s">
            <position name="bar%(copy)s" x="%(xpos)s" y="%(ypos)s" z="%(zpos)s" unit="cm"/>
            <box_dimensions z="%(z)s" x="%(x)s" y="%(y)s" lunit="cm"/>
        </parameters>'''%params

        # complete the definition of the EMR calorimeter
        data += '''
      </parameterised_position_size>
    </paramvol>
    <auxiliary auxtype="NbOfBars" auxvalue="%(NbOfBars)s"/> 
    <auxiliary auxtype="unit" auxvalue="cm"/>
    <auxiliary auxtype="BarWidth" auxvalue="%(BarWidth)s"/>
    <auxiliary auxtype="BarHeight" auxvalue="%(BarHeight)s"/>
    <auxiliary auxtype="BarLength" auxvalue="%(BarLength)s"/>
    <auxiliary auxtype="Gap" auxvalue="%(Gap)s"/>
    <auxiliary auxtype="HoleRad" auxvalue="%(HoleRad)s"/>
    <auxiliary auxtype="FiberCladdingExtRadius" auxvalue="%(FiberCladdingExtRadius)s"/>
    <auxiliary auxtype="AddWLSFiber" auxvalue="1"/>
  </volume>
    <volume name="AiDA">
      <materialref ref="AIR"/>
      <solidref ref="AiDA_solid"/>
      <physvol name="TASD_phys">
        <volumeref ref="TASD"/>
        <position name="TASD_pos" unit="mm"
                  x="0.0" y="0.0" z="0.0"/>
      </physvol>
    </volume>
  </structure>
  
  <setup name="Default" version="1.0">
    <world ref="AiDA"/>
  </setup>
</gdml>'''%vars(self)

        return data
    
    def print_EMR_to_File(self):

        data = self.print_baseline_EMR()
        data += self.generate_parvol()
        
        outfile = open(self.outfilename,'w+')
        outfile.write(data)
        outfile.close()

if __name__ == "__main__":
    gen = generate_EMR_GDML()
    gen.print_EMR_to_File()

    
