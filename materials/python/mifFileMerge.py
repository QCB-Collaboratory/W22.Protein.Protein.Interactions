#!/usr/bin/python3

import sys
import json

from os import walk
from os.path import join

from lxml import etree

mifNs = { 'm': 'http://psi.hupo.org/mi/mif'}

def getSource( mifData, mifTree ):

    if 'sourceList' not in mifData:
        mifData['sourceList'] = []
    
    # find record source
    mifSrcName = mifTree.xpath( "//m:source//m:shortLabel/text()",
                            namespaces=mifNs)
    # find record release date
    mifDate = mifTree.xpath( "//m:source/@releaseDate",
                             namespaces=mifNs)
    
    mifData['sourceList'].append( { "database": mifSrcName[0],
                                    "date": mifDate[0] } )
    return mifData


def getExperiments( mifData, mifTree ):

    mifData['exByRefId'] = {}

    #find experiments
    mifExList = mifTree.xpath( "//m:experimentDescription",
                               namespaces=mifNs)        
    for ex in mifExList:

        experiment = {}
    
        # experiment id attribute: unique *within individual file*
        exId = ex.xpath( "./@id",
                         namespaces=mifNs)
        mifData['exByRefId'][exId[0]] = experiment
    
        #pubmed
        exPmid = ex.xpath( "./m:bibref/m:xref/*[@db='pubmed']/@id",
                           namespaces=mifNs)
        if exPmid:
            experiment['pmid'] = exPmid[0]

        #experiment host
        exHost = ex.xpath( "./m:hostOrganismList/m:hostOrganism",
                           namespaces=mifNs)
        hostList = []
        experiment['hostList'] = hostList 
        for host in exHost:
            curHost = {}
            hostList.append( curHost )

            hostName = host.xpath( "./m:names/m:shortLabel/text()",
                                   namespaces=mifNs )        
            hostTaxid = host.xpath( "./@ncbiTaxId",
                                    namespaces=mifNs )
            
            curHost['name'] = hostName[0] 
            curHost['taxid'] = hostTaxid[0] 
            
        # detection method
        exDetMth = ex.xpath( "./m:interactionDetectionMethod",
                             namespaces=mifNs)
        detMthName = exDetMth[0].xpath( "./m:names/m:shortLabel/text()",
                                        namespaces=mifNs )
        detMthId = exDetMth[0].xpath( "./m:xref/*[@db='psi-mi']/@id",
                                      namespaces=mifNs )
            
        experiment['detMeth'] = { 'name': detMthName[0],
                                  'cv': 'psi-mi',
                                  'id': detMthId[0] }    
    return mifData

    
def getInteractors( mifData, mifTree ):
        
    mifData['irByRefId'] = {} 

    if 'irByAccId' not in mifData:
        mifData['irByAccId'] = {}

    #find interactors
    mifIrList = mifTree.xpath( "//m:interactor",
                               namespaces=mifNs)
    for ir in mifIrList:

        interactor = {}
    
        # primary identifier: *shoud* be unique across *all files*
        primaryIdDB = ir.xpath( "./m:xref/m:primaryRef/@db",
                                namespaces=mifNs)
        primaryIdAC = ir.xpath( "./m:xref/m:primaryRef/@id",
                                namespaces=mifNs)
    
        interactor['acc'] = primaryIdDB[0]+":"+primaryIdAC[0]
        mifData['irByAccId'][ interactor['acc'] ] = interactor
    
        # interactor id attribute: unique *only* within individual file 
        irId = ir.xpath( "./@id",
                         namespaces=mifNs)    
        mifData['irByRefId'][ irId[0] ] = interactor

        # short label
        label = ir.xpath( "./m:names/m:shortLabel/text()",
                          namespaces=mifNs)        
        if label:
            interactor['label'] = label[0]
        else:
            interactor['label'] =''

        # interactor (molecule) type
        irType = ir.xpath( "./m:interactorType",
                           namespaces=mifNs )
        irTypeName = irType[0].xpath( "./m:names/m:shortLabel/text()",
                                      namespaces=mifNs )
        irTypeId = irType[0].xpath( "./m:xref/*[@db='psi-mi']/@id",
                                    namespaces=mifNs )

        interactor['type'] = { 'name': irTypeName[0],
                               'cv': 'psi-mi',
                               'id': irTypeId[0] }
        
        # gene name
        gene = ir.xpath( "./m:names/m:alias[@type='gene name']/text()",
                         namespaces=mifNs)        
        if gene:
            interactor['gene'] = gene[0]
        else:
            interactor['gene'] = ''

        interactor['taxon'] = {'taxid': '', 'label': '', 'name': ''} 

        #taxon
        irHost = ir.xpath( "./m:organism",
                           namespaces=mifNs )
        if irHost:
            taxonSName = irHost[0].xpath( "./m:names/m:shortLabel/text()",
                                          namespaces=mifNs )        
            taxonFName = irHost[0].xpath( "./m:names/m:fullName/text()",
                                          namespaces=mifNs )        
            taxonId = irHost[0].xpath( "./@ncbiTaxId",
                                       namespaces=mifNs )

            if taxonSName:
                interactor['taxon']['label'] = taxonSName[0]

            if taxonFName:
                interactor['taxon']['name'] = taxonFName[0]

            if taxonId:
                interactor['taxon']['taxid'] = taxonId[0]
                        
    return mifData
            

def getEvidence( mifData, mifTree ):

    if 'evidList' not in mifData:    
        mifData['evidList'] = []
    
    #find evidence
    mifEvList = mifTree.xpath( "//m:interactionList/m:interaction",
                               namespaces=mifNs)
    for ev in mifEvList:
    
        evid = {}
        mifData['evidList'].append( evid )

        evid['source'] = mifData['sourceList'][-1]
        
        evid['expt'] = None
    
        # find experiment
        evExpRef = ev.xpath( ".//m:experimentRef/text()",
                             namespaces=mifNs )
        if evExpRef:  # compact MIF
            evid['expt'] = mifData['exByRefId'][ evExpRef[0] ]
        else:         # expanded MIF 
            evExpId = ev.xpath( ".//m:experimentDescription/@id",
                                namespaces=mifNs )
            if evExpId:
                evid['expt'] = mifData['exByRefId'][ evExpId[0] ]
         
        # find interaction type
        evITname = ev.xpath( "./m:interactionType//m:shortLabel/text()",
                             namespaces=mifNs )
        evITid = ev.xpath( "./m:interactionType/m:xref/*[@db='psi-mi']/@id",
                           namespaces=mifNs )
    
        evid['intType'] = { 'name': evITname[0],
                            'cv': 'psi-mi',
                            'id': evITid[0] }

        #find imex
        imexid = ev.xpath( "./@imexId",
                         namespaces=mifNs )
        if imexid:
            evid['imexid'] = imexid[0]
        else:
            evid['imexid'] = ''
            
        evid['partList'] = []

        # find participants
        evParts = ev.xpath( ".//m:participantList/m:participant",
                            namespaces=mifNs )

        for p in evParts:
            part = {}
            evid['partList'].append( part )

            # find interactor
            interactor = None
            intRefId = p.xpath( "./m:interactorRef/text()",
                                namespaces=mifNs )
            if intRefId:  # compact MIF
                interactor =  mifData['irByRefId'][ intRefId[0] ]
            else:         # expanded MIF 
                intId = p.xpath( ".//m:interactor/@id",
                                 namespaces=mifNs )
                if intId:
                    interactor =  mifData['irByRefId'][ intId[0] ]
            part['interactor'] = interactor

            # experimental/biological role
            # (bait, prey, ancillary, enzyme, substrate, etc.)
        
            part['roleList'] = []
            roleList = p.xpath( ".//*[local-name(.)='experimentalRole' or" +  
                                "     local-name(.)='biologicalRole'  ]",
                                namespaces=mifNs )
            for r in roleList:
                roleName = r.xpath( "./m:names/m:shortLabel/text()",
                                    namespaces=mifNs )
                roleId = r.xpath( "./m:xref/*[@db='psi-mi']/@id",
                                  namespaces=mifNs )

                if roleName[0] != 'unspecified role':
                    role = { 'name': roleName[0],
                             'cv': 'psi-mi',
                             'id': roleId[0] }                
                    part['roleList'].append( role )

    return mifData

def dropInternalRefs( mifData ):
    mifData.pop( 'exByRefId', None )
    mifData.pop( 'irByRefId', None )
        
def buildNetwork( mifData, expand=None):
    
    # expand:
    # None   - binary edges only
    # spoke  - bait/prey (spoke) expansion 
    # matrix - spoke + matrix of physical and below, ignoring ancillaries

    # network data (note: participant details are losts !!!)
    edges = {} # key: (acc1,acc2); value = {'evidList': [ev1,ev2,...], 'intType':[...]}
    
    for ev in mifData['evidList']:        
        bait = None
        preyDict = {}
        partDict = {}
        for p in ev['partList']:        
            skip = False
            for r in p['roleList']:                
                if r['name'] == 'bait':
                    bait = p['interactor']['acc']
                if r['name'] == 'prey':
                    if p['interactor']['acc'] in preyDict:
                        preyDict[ p['interactor']['acc'] ] += 1
                    else:
                        preyDict[ p['interactor']['acc'] ] = 1
                        
                if r['name'] == 'ancillary':
                    skip = True
            if not skip:
                if p['interactor']['acc'] in partDict:
                    partDict[ p['interactor']['acc'] ] += 1
                else:
                    partDict[ p['interactor']['acc'] ] = 1
        
        if len(partDict) == 2:
            (part1, part2) = ( partDict.keys() )
            if part1 > part2:
                (part2, part1) = (part1, part2) 
            if (part1, part2) in edges:
                edges[ (part1, part2) ]['evidList'].append( ev )
                edges[ (part1, part2) ]['intType'].append( ev['intType'] )
                
            elif part1 != part2:
                edges[ (part1, part2) ] = { 'evidList': [ ev ],
                                            'intType': [ev['intType']] }

        elif len(partDict) > 2:
            if expand in ['spoke','matrix']:
                if bait is not None:
                    for prey in preyDict.keys():

                        itype = {'cv': 'psi-mi', 'id': 'MI:0915',
                                 'name': 'physical association' }
                        
                        if bait > prey:
                            (part1, part2) = (prey, bait)
                        else:
                            (part1, part2) = (bait, prey)
                        if (part1, part2) in edges:
                            edges[(part1, part2)]['evidList'].append( ev )
                            edges[(part1, part2)]['intType'].append( itype )
                        elif part1 != part2:
                            edges[ (part1, part2) ] = { 'evidList': [ ev ],
                                                        'intType': [itype] }                            
                            
            if expand in ['matrix'] and ev['intType']['name'] not in ['association']:
                for part1 in partDict.keys():
                    for part2 in partDict.keys():
                        if part2 > part1:
                            if (part1, part2) in edges:
                                edges[ (part1,part2) ]['evidList'].append( ev )
                                edges[ (part1,part2) ]['intType'].append( ev['intType'] )
                            else:
                                edges[ (part1,part2) ] = {'evidList': [ ev ],
                                                          'intType': [ev['intType']] }

    return edges    


def toJson( mifData,  indent=None ):
    return json.dumps( mifData, sort_keys=True, indent=indent )


def toSif( mifData, network ):
    buffer = ''
    for edge in network:
        (p1,p2) = edge

        p1type = mifData['irByAccId'][p1]['type']['name'][0]
        p2type = mifData['irByAccId'][p2]['type']['name'][0]        
        
        buffer += "%s %s%s %s\n" % ( p1, p1type, p2type, p2 ) 
    return buffer


def toEdgePoperties( mifData, network ):
    buffer = '\t'.join( ['#acc', 'type_id', 'type_name',
                         'method_id', 'method_name', 'pmid','imexid',
                         'method_count','pmid_count','imexid_count'])+'\n'

    for edge in network:
        (p1,p2) = edge
        evList = network[edge]['evidList']
        itypeList = network[edge]['intType']

        columns =[]
        
        p1type = mifData['irByAccId'][p1]['type']['name'][0]
        p2type = mifData['irByAccId'][p2]['type']['name'][0]        
        
        columns.append("%s (%s%s) %s" % (p1, p1type, p2type, p2))

        #  type 
        typeId=''
        typeName=''
        for itype in itypeList:
            typeId += itype['id'] + '|'
            typeName += itype['name'] + '|'

        columns.extend( [typeId[:-1], typeName[:-1]] )
        
        # method
        methodId=''
        methodName=''
        methodDict={}
        
        for ev in evList:        
            methodId +=   ev['expt']['detMeth']['id'] + '|'
            methodName += ev['expt']['detMeth']['name'] + '|'
            methodDict[ ev['expt']['detMeth']['id'] ] = True
            
        columns.extend( [methodId[:-1],  methodName[:-1]] )

        # pmid
        pmid = ''
        pmidDict = {}
        for ev in evList:
            pmid += ev['expt']['pmid'] + '|'
            pmidDict[ev['expt']['pmid']] = True
        columns.append( pmid[:-1] )

        # imexid
        imexid = ''
        imexidDict = {}
        for ev in evList:
            imexid += ev['imexid'] + '|'
            imexidDict[ev['imexid']] = True

        columns.append( imexid[:-1] )

        columns.extend( [str( len(methodDict) ),
                         str( len(pmidDict) ),
                         str( len(imexidDict) )])
        
        buffer += '\t'.join(columns) + '\n'
        
    return buffer


def toNodePoperties( mifData, network ):
    buffer = '#acc\tlabel\ttype\ttaxid\ttaxon_label\ttaxon_name\tgene\n'

    nodeDict = mifData['irByAccId']
    
    for node in nodeDict:
        #acession
        line = "%s\t" % ( nodeDict[node]['acc'] )
        
        # label
        line += "%s\t" %( nodeDict[node]['label'] ) 

        # type
        line += "%s\t" %( nodeDict[node]['type']['name'] ) 

        # taxid
        line += "%s\t" %( nodeDict[node]['taxon']['taxid'] ) 

        # taxon label
        line += "%s\t" %( nodeDict[node]['taxon']['label'] ) 

        # taxon name
        line += "%s\t" %( nodeDict[node]['taxon']['name'] ) 

        # gene
        line += "%s\t" %( nodeDict[node]['gene'] ) 

        buffer += line + "\n"        
    return buffer


#------------------------------------------------------------------------------
# merge files
#------------

mifData = {} 
assert len(sys.argv) == 4
mydir = sys.argv[1]
mymth = sys.argv[2]
myout = sys.argv[3]

myfiles = []
for (dirpath, dirnames, filenames) in walk( mydir ):
    myfiles.extend(filenames)
    break

for f in myfiles:
    mifTree = etree.parse( join( mydir, f ) )

    getSource( mifData, mifTree )
    getExperiments( mifData, mifTree )
    getInteractors( mifData, mifTree )
    getEvidence( mifData, mifTree )
    dropInternalRefs( mifData )

#------------------------------------------------------------------------------
# operations 
#------------
    
if myout == 'json':
    print( toJson(mifData, indent=1) )
    
if mymth in ['none','spoke','matrix']:
    network = buildNetwork( mifData, expand=mymth )

if myout == 'sif':
    print( toSif( mifData, network) )
elif myout == 'nodeprop':
    print( toNodePoperties( mifData, network) )
elif myout == 'edgeprop':
    print( toEdgePoperties( mifData, network) )

