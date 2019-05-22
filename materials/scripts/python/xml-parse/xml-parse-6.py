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


def parseParticipant( mif, dom ):    

    global nc

    p11t = {}

    refID =  dom.xpath( "./@id",
                        namespaces = ns )[0]
        
    # interactor
    #-----------
        
    iref = dom.xpath( "./mif:interactorRef/text()",
                      namespaces = ns )

    if len( iref ) == 1:
        iref = iref[0]
        p11t['iref'] = iref
        p11t['i10r'] = mif["i10r"][iref]
    else:
        p11t['iref'] = "N/A"
            
    if p11t['iref'] == "N/A":
            
        i10rl = entry.xpath( ".//mif:interactor",
                             namespaces = ns )
        if i10rl:
            p11t['i10r'] =  parseInteractor( i10rl[0] )

    # bio role
    #---------
        
    bioRoleDom = dom.xpath( "./mif:biologicalRole",
                            namespaces = ns )
    p11t['brole'] = parseCV( bioRoleDom[0] ) 

    
    # exp role
    #---------
        
    erole = dom.xpath( ".//mif:experimentalRole",
                       namespaces = ns )
        
    if len(erole) > 0:
        p11t['erole'] = parseCV( erole[0] )
                     
    # exp prep role
    #--------------
        
    eprep = dom.xpath( ".//mif:experimentalPreparation",
                       namespaces = ns )
    if len( eprep ) > 0:
        p11t['eprep'] = parseCV( eprep[0] )

        
    # features
    #---------
        
    ftrDom = dom.xpath( ".//mif:feature",
                        namespaces = ns )
        
    #p11t['feature'] = parseFeature( ftrDom )
        
    return p11t


def parseInteraction( mif, dom ):

    global ns

    i11n = {}
        
    refID =  dom.xpath( "./@id",
                        namespaces = ns )[0]        
    
    i11n['rid'] = refID

    imexID =  dom.xpath( "./@imexId",
                         namespaces = ns )

    if imexID:
        i11n['imex'] = imexID[0]
    else:
        i11n['imex'] = ""


    # interaction type
    #-----------------
        
    itype = parseCV( dom.xpath( ".//mif:interactionType",
                                namespaces = ns )[0] )
                                
    i11n["itype"] = itype

    
    # experiment ref list
    #--------------------

    e10tl = dom.xpath( "./mif:experimentList",
                       namespaces = ns )[0]

    e10t = []
    
    for e in e10tl:

        if e.xpath('local-name(.)') ==  'experimentRef':

            print(" Compresed MIF")
            e10t.append( mif["e10t"][ e.xpath('./text()')[0] ])

        else:
            
            print(" Expanded MIF")            
            e10t.append( parseExperiment( e ) )

    i11n["e10t"] = e10t

    
    # participant list
    #-----------------

    p11tl = dom.xpath( "./mif:participantList/mif:participant",
                       namespaces = ns )
    p11t = []
    
    for p in p11tl:
        p11t.append( parseParticipant( mif, p ) )
            
    i11n['p11t'] = p11t

    return i11n
            
def parseExperiment( dom ):

    global ns

    e10t = {}
            
    # <experimentDescription id="11175">
        
    e10t['rid'] =  dom.xpath( "./@id",
                              namespaces = ns )[0]
    
    #     <names>
    #         <shortLabel>volkman-2002-1</shortLabel>
    #         <fullName>Structure of the N-WASP EVH1 domain-WIP complex: insight into the molecular basis of Wiskott-Aldrich Syndrome.</fullName>
    #     </names>
    #     <bibref>
    #         <xref>
    #             <primaryRef db="pubmed" dbAc="MI:0446" id="12437929" refType="primary-reference" refTypeAc="MI:0358"/>
    #         </xref>
    #     </bibref>
        
    e10t['pmid'] = dom.xpath( ".//mif:bibref//mif:primaryRef/@id",
                              namespaces = ns )[0]
            
    #     <xref>
    #         <primaryRef db="mint" dbAc="MI:0471" id="MINT-722978" refType="identity" refTypeAc="MI:0356"/>
    #         <secondaryRef db="pubmed" dbAc="MI:0446" id="12437929" refType="secondary-ac" refTypeAc="MI:0360"/>
    #         <secondaryRef db="doi" dbAc="MI:0574" id="10.1016/S0092-8674(02)01076-0" refType="secondary-ac" refTypeAc="MI:0360"/>
    #         <secondaryRef db="intact" dbAc="MI:0469" id="EBI-8558582" refType="identity" refTypeAc="MI:0356"/>
    #         <secondaryRef db="mint" dbAc="MI:0471" id="MINT-5215918" refType="primary-reference" refTypeAc="MI:0358"/>
    #     </xref>
    
    # Add Xreferences if needed

    #     <hostOrganismList>
    #         <hostOrganism ncbiTaxId="-1">
    #             <names>
    #                 <shortLabel>in vitro</shortLabel>
    #                 <fullName>In vitro</fullName>
    #             </names>
    #         </hostOrganism>
    #     </hostOrganismList>

    e10t['host'] = {'ns':'taxid'}

    e10t['host']['ac'] = dom.xpath( ".//mif:hostOrganism/@ncbiTaxId",
                                    namespaces = ns )[0]
            
    e10t['host']['label'] = dom.xpath( ".//mif:hostOrganism//mif:shortLabel/text()",
                                       namespaces = ns )[0]

    name = dom.xpath( ".//mif:hostOrganism//mif:fullName/text()",
                      namespaces = ns ) 
    if name:
        e10t['host']['name'] = name[0]
    else:
        e10t['host']['name'] = e10t['host']['label']

    #     <interactionDetectionMethod>
    #         <names>
    #             <shortLabel>x-ray diffraction</shortLabel>
    #             <fullName>x-ray crystallography</fullName>
    #             <alias typeAc="MI:0303" type="go synonym">X-ray</alias>
    #             <alias typeAc="MI:1041" type="synonym">X-ray</alias>
    #         </names>
    #         <xref>
    #             <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0114" refType="identity" refTypeAc="MI:0356"/>
    #             <secondaryRef db="intact" dbAc="MI:0469" id="EBI-1272" refType="identity" refTypeAc="MI:0356"/>
    #             <secondaryRef db="pubmed" dbAc="MI:0446" id="14755292" refType="primary-reference" refTypeAc="MI:0358"/>
    #         </xref>
    #     </interactionDetectionMethod>
            
    intIdMthDom = dom.xpath( ".//mif:interactionDetectionMethod",
                             namespaces = ns )
            
    e10t['intMth'] = parseCV( intIdMthDom[0] )
                        
    #     <participantIdentificationMethod>
    #         <names>
    #             <shortLabel>predetermined</shortLabel>
    #             <fullName>predetermined participant</fullName>
    #             <alias typeAc="MI:1041" type="synonym">predetermined</alias>
    #         </names>
    #         <xref>
    #             <primaryRef db="psi-mi" dbAc="MI:0488" id="MI:0396" refType="identity" refTypeAc="MI:0356"/>
    #             <secondaryRef db="intact" dbAc="MI:0469" id="EBI-1465" refType="identity" refTypeAc="MI:0356"/>
    #             <secondaryRef db="pubmed" dbAc="MI:0446" id="14755292" refType="primary-reference" refTypeAc="MI:0358"/>
    #         </xref>
    #     </participantIdentificationMethod>
    
    prtIdMthDom = dom.xpath( ".//mif:participantIdentificationMethod",
                             namespaces = ns )
        
    e10t['prtMth'] = parseCV( prtIdMthDom[0] )            
         
    return e10t

def parseSource( msrc ):

    global ns
    
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

    global ns

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
    
    i10r[ 'pxref' ] = parseRefList( prXRefDom )
    
    scXRefDom = dom.xpath( "./mif:xref/mif:secondaryRef",
                           namespaces = ns )

    i10r[ 'sxref' ] = parseRefList( scXRefDom )

    return i10r


def parseCV( cvdom ):

    global ns
    
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
    
    
def parseRefList( ref ):

    global ns
    
    ret = []
    for r in ref:
        acID = r.xpath( "./@id",
                        namespaces = ns )
        nsID = r.xpath( ".//@db",
                        namespaces = ns )
 
        nsDB = r.xpath( ".//@dbAc",
                        namespaces = ns )
        verID  = r.xpath( "@version",
                          namespaces = ns )
        ref = {}
        if acID:
            ref['ac'] = acID[0]
        if nsID:
            ref['ns'] = nsID[0]
        if nsDB:
            ref['nsAc'] = nsID[0]
        if verID:
            ref['ver'] = verID[0]

        ret.append( ref )
            
    return ret
    
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
    e10t = parseExperiment( cexp )
    mif["e10t"][ e10t['rid'] ] = e10t

    
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
    ii = parseInteraction( mif, cintn )
    mif["i11n"][ii['rid'] ] = ii

#print( ET.tostring(entry) )

print(json.dumps(mif, indent=2))



