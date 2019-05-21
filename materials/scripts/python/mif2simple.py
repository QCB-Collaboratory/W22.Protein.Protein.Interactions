import sys
sys.path.insert( 0, "/cluster1/opt/python/lib/")

###############################################################################
#  mif254 to sif/property file converter
# 
#   input: mif file (verison 254)   
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

parser = argparse.ArgumentParser(description='mif2simple')
parser.add_argument( '--file',  dest="file", type=str, required=False,
                    default="", help='Input MIF(v254) file.')

parser.add_argument('--out',  dest="out", type=str, required=False,
                    default="sif", help='Output type (default sif)')
args = parser.parse_args()

# default input file

file='ATP/mif25/29650704.mif25'

if args.file != "":
    file = args.file

# create mif file parser

mifParser = mif25.Parser()

of = open( file, 'r' )

# parse mif file

mif = mifParser.parse( of )

pmap = {}   # dictionary of participants (indexed by primary crossreference )

# process interactions

for int_id in mif["i11n"].keys():

    i11n = mif["i11n"][int_id]

    # get experimental method and interaction type

    method = i11n["e10t"][0]["intMth"]
    itype  = i11n["itype"]
    
    plist= []
    llist = []
    xlist = []

    # prepare a list of participants

    for part in mif["i11n"][int_id]["p11t"]:        

        alist = part["i10r"]["alias"];
        gname = ""
        for a in alist:
            if a['type'] == 'gene name':
                gname = a['value']
                
        plist.append( {"erole":part["erole"]["label"], 
                       "label":part["i10r"]["label"],
                       "gene": gname,
                       "type":part["i10r"]['type']['label'],
                       "pxref":part["i10r"]["pxref"]} )

        llist.append( part["i10r"]["label"] )

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

                if p1["pxref"][0]["ac"] > bait:
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
