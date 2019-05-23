import sys
sys.path.insert( 0, "./pylib/")  

# Note: to run from any directory use full path when pointing sys.path to
#       pylib directory above.    

###############################################################################
#  mif254 to sif/property file converter
# 
#   input: mif file (verison 254) (or directory)
#          (optional) acession list 
#  
#   output: sif or property files (iprop - interactions, nprop - proteins)
#
#   Note: Only associations and physial associations are recognized/processed
#         For simplicty mif25 library converts mif file into a generic, nested
#         list/dictionary data structure. More mature implementation might/
#         should define separate classes for participants(proteins), features, 
#         cv terms, interactions.         
#
###############################################################################

import argparse
from qcb import mif25

# default input file

file='../../ATP/mif25/29650704.mif25'

parser = argparse.ArgumentParser(description='mif2simple')
parser.add_argument( '--mif-file',  dest="mfile", type=str, required=False,
                     default="", help='Input MIF(v254) file.')

parser.add_argument( '--mif-dir',  dest="mdir", type=str, required=False,
                     default="", help='MIF file(s) directory. Compressed files OK.')

parser.add_argument( '--acc-file',  dest="acc", type=str, required=False,
                     default= "", help='Accessions to test')

parser.add_argument( '--out',  dest="out", type=str, required=False,
                     default="sif", help='Output type (default sif)')

args = parser.parse_args()

# accession fileter list

afile = args.acc

alist = []

if afile != "":
    af = open(afile,"r")
    for afl in af:
        alist.append(afl.strip())

print( alist )



# ensure only one of mif-file and mif-dir provided

if args.mfile != "" and args.mdir != "":
    sys.exit("ERROR: Only one of --mif-file and --mif-dir can be specified") 

# build input file list

files = []

if args.mdir == "":
    files.append( file )
else:
    for (dirpath, dirnames, filenames) in walk( args.mdir):
        for file in filenames:        
            fname = dirpath + '/' + file
            files.append( fname )

i = 0

mif = {}

# create mif file parser

mifParser = mif25.Parser()


# parse files

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
            source.append( open( file, 'r' ) )
            
    if source:
        print( "%d %s" % (i, file) )

        for cs in source:
            m = mifParser.parse( cs )
            ek0 = list(m["e10t"].keys())[0]
            pmid = m["e10t"][ek0]["pmid"]
            
            mif[pmid] = m
    
# go over parsed files

for pmid, cmif in mif.items():

    print( "%s %s" % (pmid, "file") )

    # process interactions

    for int_id in cmif["i11n"].keys():

        pmap = {}   # dictionary of participants (indexed by primary crossreference )

        i11n = cmif["i11n"][int_id]

        # get experimental method and interaction type

        method = i11n["e10t"][0]["intMth"]
        itype  = i11n["itype"]
    
        plist= []
        llist = []
        xlist = []
        pxlist= []

        # prepare a list of participants

        for part in cmif["i11n"][int_id]["p11t"]:        
            #print("***")
            gname = ""       
            if 'alias' in part["i10r"]:
                alist = part["i10r"]["alias"];

                for a in alist:
                    if a['type'] == 'gene name':
                        gname = a['value']
                
            plist.append( {"erole":part["erole"]["label"], 
                           "label":part["i10r"]["label"],
                           "gene": gname,
                           "type":part["i10r"]['type']['label'],
                           "pxref":part["i10r"]["pxref"]} )

            llist.append( part["i10r"]["label"] )
            pxlist.append( part["i10r"]["pxref"][0]["ac"] )

        # test if in accesion list

        keep = False
        if alist:
            for facc in alist:
                if facc in pxlist:
                    keep = True
                    break
        else:
            keep = True

        if keep:
        
            # Process associations (spoke expansion) and print out sif or property info
            # Note:  interaction type changes when spoke-expanding   

            if i11n["itype"]["label"] == "association":
            
                bait =""
                pbait = {}
                for p in plist:                
                    if p["erole"] == "bait":                
                        bait =  p["pxref"][0]["ac"]
                        pbait = p                    
                        break

                for p1 in plist:
                    if p1["pxref"][0]["ac"] != bait:

                        pmap[bait] = pbait
                        pmap[p1["pxref"][0]["ac"]] = p1

                    if p1["pxref"][0]["ac"] < bait:
                        if args.out == "sif":
                            print( " ".join((bait, ' pp ',  p1["pxref"][0]["ac"] )))
 
                        if args.out == "iprop":                       
                            ename = " ".join((bait, '(pp)',  p1["pxref"][0]["ac"] ))                       
                            print( "\t".join((ename, "physical association", method["label"])))
                        
                    else:
                        if args.out == "sif":
                            print( " ".join((p1["pxref"][0]["ac"], ' pp ', bait)) )

                        if args.out == "iprop":                       
                            ename = " ".join((p1["pxref"][0]["ac"], '(pp)', bait))                       
                            print( "\t".join((ename, "physical association", method["label"])))

            # Process physical associations (matrix expansion)
            # Note: interaction type does not change

            if i11n["itype"]["label"] == "physical association":
                for p1 in plist:
                    for p2 in plist:
                        if p1["pxref"][0]["ac"] <  p2["pxref"][0]["ac"]:

                            pmap[p1["pxref"][0]["ac"]] = p1
                            pmap[p2["pxref"][0]["ac"]] = p2


                            if args.out == "sif":
                                print( " ".join((p1["pxref"][0]["ac"], ' pp ',  p2["pxref"][0]["ac"] )))
                    
                            if args.out == "iprop":
                                ename =  " ".join((p1["pxref"][0]["ac"], '(pp)',  p2["pxref"][0]["ac"] ))
                                print( "\t".join((ename, itype["label"], method["label"])))

        if args.out == "pprop":
            for pid, p in pmap.items():
                print("\t".join((pid, p["label"],p["gene"])))
