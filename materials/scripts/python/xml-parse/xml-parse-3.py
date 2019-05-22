import sys
sys.path.insert(0,"/cluster1/opt/python/lib/")

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

def parseSource( msrc ):
    
    source = {}

    # source
    #-------
    
    #  <source releaseDate="2013-09-02+01:00">
    #   <names>
    #    <shortLabel>MINT</shortLabel>
    #    <fullName>MINT, Dpt of Biology, University of Rome Tor Vergata</fullName>
    #   </names>
    #   <xref>
    #    <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0471" refType="identity" refTypeAc="MI:0356"/>
    #    <secondaryRef db="intact" dbAc="MI:0469" id="EBI-1579228" refType="identity" refTypeAc="MI:0356"/>
    #   </xref>
    #   <attributeList>
    #    <attribute name="url" nameAc="MI:0614">http://mint.bio.uniroma2.it/mint</attribute>
    #    <attribute name="url">http://mint.bio.uniroma2.it/mint</attribute>
    #   </attributeList>
    #  </source>
        
    source['label'] = msrc.xpath( ".//mif:shortLabel/text()",
                                  namespaces = ns )[0]
    
    name =  msrc.xpath( ".//mif:fullName/text()",
                        namespaces = ns )
    if name:
        source['name'] = name[0]
    else:
        source['name'] = source['label']
        
    source['ns'] = msrc.xpath( ".//mif:primaryRef/@db",
                               namespaces = ns )[0]
    source['ac'] = msrc.xpath( ".//mif:primaryRef/@id",
                               namespaces = ns )[0]
    return source


def parseInteractor( dom ):

    idix = {}
    i10r = {}
    
    i10r['rid'] =  dom.xpath( "./@id",
                              namespaces = ns )[0]

    SL =  dom.xpath( "./mif:names/mif:shortLabel/text()",
                     namespaces = ns )
    if len( SL ) > 0:
        i10r[ 'label' ] = SL[0]
    else:
        i10r[ 'label' ] = ''
        
    FN =  dom.xpath( "./mif:names/mif:fullName/text()",
                     namespaces = ns )
    if len( FN ) > 0:
        i10r[ 'name' ] = FN[0]
    else:
        i10r[ 'name' ] = ''
            
    i10r[ 'alias' ] =  dom.xpath( "./mif:names/mif:alias/text()",
                                  namespaces = ns )
    
    itpDom = dom.xpath( "./mif:interactorType",
                        namespaces = ns )
    
    i10r[ 'type' ] =  parseCV( itpDom[0] )

    
    prXRefDom = dom.xpath( "./mif:xref/mif:primaryRef",
                           namespaces = ns )

    #i10r[ 'pxref' ] = parseRefList( prXRefDom )
    
    scXRefDom = dom.xpath( "./mif:xref/mif:secondaryRef",
                           namespaces = ns )

    #i10r[ 'sxref' ] = parseRefList( scXRefDom )
        
    
    
    return i10r


def parseCV( cvdom ):

    cvID = cvdom.xpath( "./mif:xref/mif:primaryRef/@id",
                        namespaces = ns )
    cvDB = cvdom.xpath( "./mif:xref/mif:primaryRef/@db",
                        namespaces = ns )
    cvDBID = cvdom.xpath( "./mif:xref/mif:primaryRef/@dbAc",
                          namespaces = ns )
    cvSL = cvdom.xpath( "./mif:names/mif:shortLabel/text()",
                        namespaces = ns )
    cvFN = cvdom.xpath( "./mif:names/mif:fullName/text()",
                        namespaces = ns )
    
    return {"ns":cvDB[0],"ac":cvID[0], "label":cvSL[0],"name":cvSL[0],"nsAc":cvDBID[0] }
    



    
# parse source
#-------------

source = entry.xpath( "./mif:source", namespaces = ns )[0]
    
mif["source"]  = parseSource( source )

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
    ii = parseInteractor( cintr )
    mif["i10r"][ ii['rid'] ] = ii
        
    
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



