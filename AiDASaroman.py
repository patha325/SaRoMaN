import saroman, os, sys

if __name__=="__main__":
    s = saroman.saroman()
    s.xml_file_path = os.path.join(s.exec_base,'AiDA_TASD.gdml')
    s.Handle_commandline_input(sys.argv[1:])
