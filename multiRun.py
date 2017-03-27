#######################################################################################################################
#General python import
#######################################################################################################################
import sys
import random
import shutil
import time
#######################################################################################################################
#Importing own python files
#######################################################################################################################
import saroman
#######################################################################################################################
#Class generation
#######################################################################################################################

#######################################################################################################################
#File specific functions
#######################################################################################################################

if __name__ == "__main__":
    s=saroman.saroman()
    #s.home = '/afs/phas.gla.ac.uk/user/p/phallsjo/SaRoMaN'
    #s.exec_base = os.path.join(self.home, 'SaRoMaN')
    #s.out_base  = os.path.join(self.home, 'batch')
    #s.seed = 500 * random.random()
    s.third_party_support = '/data/neutrino05/phallsjo/third_party'
    #s.generate_field_map = False
    s.seed = 1000
    s.Nevts = 1000
    s.part = 'mu+'#'14'
    s.pid = 13

    def my_range(start, end, step):
        while start <= end:
            yield start
            start += step
    
    for x in my_range(5.1, 6.1, 0.1):
        s.part_eng_min = x #4.0
        s.part_eng_max = x #4.0
    
        s.Handle_commandline_input(sys.argv[1:])

        s.base ='/data/neutrino05/phallsjo/copy/SaRoMan/out/rec_out/nd_'+s.part+'CC/'

        time.sleep(2)

        shutil.move(s.base+'nd_'+s.part+'CC_'+str(s.Nevts)+'.root',
                    s.base+'MIND/'+str(s.part_eng_min)+'_'+str(s.Nevts)+'.root')

        time.sleep(2)
    

#######################################################################################################################
