import sys

from lxml import etree as ET
import json

#test_mif25='ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/pmid/2019/15138291.xml'
test_mif25="/home/lukasz/AAtWork/DIP/QCB/PPI2019/files/mif25/11034202.xml"

#source = sys.argv[1]
source = test_mif25

parser = ET.XMLParser( remove_blank_text=True ) 
dom = ET.parse( source, parser )
     
#<entrySet level="2" version="5" minorVersion="4" 
#          xsi:schemaLocation="http://psi.hupo.org/mi/mif http://psidev.sourceforge.net/mi/rel25/src/MIF254.xsd" 
#          xmlns="http://psi.hupo.org/mi/mif" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
#           
#         abc:xmlns="http://psi.hupo.org/mi/mif" 

ns = { 'mif': 'http://psi.hupo.org/mi/mif' }

elist = dom.xpath( '//mif:entry', namespaces = ns)

entry = elist[0] # only one entry per file in IMEx primary records data

mif = {}

# parse source
#-------------

source = entry.xpath( "./mif:source", namespaces = ns )[0]

print(ET.tostring(source))
    
mif["source"]  = {}  #  = parseSource( source )

# parse experiments
#------------------

expl = entry.xpath( ".//mif:experimentDescription",
                    namespaces = ns )

mif["e10t"] = {}

for cexp in expl:
    e10t = {} # parseExperiment( cexp )
    
    #mif["e10t"][ e10t['rid'] ] = e10t

    
# parse interactors
#------------------
    
i10rl = entry.xpath( ".//mif:interactor",
                     namespaces = ns )
mif["i10r"] = {}    

for cintr in i10rl:
    ii = {} # parseInteractor( cintr )

    # mif["i10r"][ ii['rid'] ] = ii
        
    
# parse interactions
#-------------------
        
i11nl = entry.xpath( ".//mif:interaction",
                     namespaces = ns )

mif["i11n"] = {} 
       
for cintn in i11nl:
    ii = {} # parseInteraction( mif, cintn )
    #mif["i11n"][ii['rid'] ] = ii

#print( ET.tostring(entry) )

print(json.dumps(mif, indent=2))

