import sys
sys.path.insert(0,"/home/lukasz/AAtWork/DIP/QCB/PPI2019/scripts/python/pylib")

from os import walk
from zipfile import ZipFile

import argparse
from qcb import mif25

parser = argparse.ArgumentParser(description='MIF Reader')
parser.add_argument('--source-dir',  dest="sdir", type=str, required=True,
                    help='MIF file(s) directory. Compressed files OK.')

parser.add_argument('--acc',  dest="alist", nargs='*', type=str, required=True,
                    help='Accessions to test')
args = parser.parse_args()

mifParser = mif25.Parser()

alist = args.alist

files = []

for (dirpath, dirnames, filenames) in walk( args.sdir):
    for file in filenames:        
        fname = dirpath + '/' + file
        files.append( fname )
        
i = 0

mif = {}

for file in files:

    i = i + 1 
    aclist = {}
    source = []
    
    if file.find("unassigned") < 0:
        if file.endswith( ".zip" ):
            myzip = ZipFile( file, 'r' )

            for sl in myzip.namelist():
                if  sl.find("negative") < 0 :            
                    source.append( myzip.open( sl, 'r' ) )
        else:
            source.append( open(file, 'r' ) )
            
    if source:
        print( i, file )

        for cs in source:
            m = mifParser.parse( cs )
            ek0 = list(m["e10t"].keys())[0]
            pmid = m["e10t"][ek0]["pmid"]
            
            mif[pmid] = m

expid = {}
llist = [] 
for pmid, cmif in mif.items():

    print("processing pmid: " + pmid)
    sys.stdout.flush()

    for iid, i11n in cmif["i11n"].items():
        plist= []
        for p11t in i11n["p11t"]:
            plist.append({"erole":p11t["erole"]["label"], "label":p11t["i10r"]["label"]})
            llist.append( p11t["i10r"]["label"] )

        # test if in accesion list

        keep = False

        for facc in alist:
            for pacc in llist:
                if facc == pacc:
                    keep = True
                    break
        if keep:

            uuid = str(pmid) + "_" + i11n["rid"]

            for p1 in plist:
                key = p1["label"]

                if not key in expid:
                    expid[key]= []

                expid[key].append( uuid )

print( "key count: ",len(list(expid.items())))    
print("---- property table starts here ----")

for k, v in expid.items():
    print( k, "\t", "|".join(v), "\t", len(v) )


