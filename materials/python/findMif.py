#!/usr/bin/python3

import sys
import zipfile

from os import walk
from os.path import join

assert len(sys.argv) == 3
mydir = sys.argv[1]
mylist = sys.argv[2]

# read identifier list
fh = open(sys.argv[2], "r")

idList= []
for l in fh:
    idList.append( '"' +l.strip()+'"' ) 
fh.close()

print( "Identifiers: " + str( idList ))

# scan zipfiles

for (dirpath, dirnames, filenames) in walk( mydir ):
    #myfiles.extend(filenames)
    for f in filenames:
        #print( join( dirpath, f ) )

        zf = zipfile.ZipFile( join( dirpath, f ) )
        zfNames = zf.namelist()
        if len(zfNames ) > 1:
            print( "Skip (multi):" +  join( dirpath, f ))
        else:
            done = False
            lines = zf.read(zfNames[0]).decode('utf-8').split('\n')
            for l in lines:
                for id in idList:
                    if id in l:
                        print( join( dirpath, f )+"\t" + id )
                        done = True
                        break;
                if done:
                    break;
            if not done:
                print( "Skip (miss):" +  join( dirpath, f ))
