#!/usr/bin/python3

import sys
from lxml import etree

mifNs = { 'm': 'http://psi.hupo.org/mi/mif'}
assert len(sys.argv) == 2
mifFile = sys.argv[1]

mifTree = etree.parse( mifFile )

# find record source
mifSrcName = mifTree.xpath( "//m:source//m:shortLabel/text()",
                            namespaces=mifNs)
print( "Record Source: %s" % ( mifSrcName[0] ))

# find record release date
mifDate = mifTree.xpath( "//m:source/@releaseDate",
                         namespaces=mifNs)
print("Release Date: %s" % ( mifDate[0] ))

exByRefId = {} 

#find experiments
mifExList = mifTree.xpath( "//m:experimentDescription",
                           namespaces=mifNs)
for ex in mifExList:

    experiment = {}
    
    # experiment id attribute: unique *within individual file*
    exId = ex.xpath( "./@id",
                     namespaces=mifNs)
    exByRefId[ exId[0] ] = experiment
    
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
    
    
irByRefId = {} 
irByAccId = {}

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
    irByAccId[ interactor['acc'] ] = interactor
    
    # interactor id attribute: unique *only* within individual file 
    irId = ir.xpath( "./@id",
                     namespaces=mifNs)    
    irByRefId[ irId[0] ] = interactor

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
        
evidList = []
    
#find evidence
mifEvList = mifTree.xpath( "//m:interactionList/m:interaction",
                           namespaces=mifNs)
for ev in mifEvList:
    
    evid = {}
    evidList.append( evid )
    expt = None
    
    # find experiment
    evExpRef = ev.xpath( ".//m:experimentRef/text()",
                         namespaces=mifNs )
    if evExpRef:  # compact MIF
        expt = exByRefId[ evExpRef[0] ]
    else:         # expanded MIF 
        evExpId = ev.xpath( ".//m:experimentDescription/@id",
                            namespaces=mifNs )
        if evExpId:
            expt = exByRefId[ evExpId[0] ]

    evid['expt'] = expt

    # find interaction type
    evITname = ev.xpath( "./m:interactionType//m:shortLabel/text()",
                         namespaces=mifNs )
    evITid = ev.xpath( "./m:interactionType/m:xref/*[@db='psi-mi']/@id",
                       namespaces=mifNs )
    
    evid['IntType'] = { 'name': evITname[0],
                        'cv': 'psi-mi',
                        'id': evITid[0] }
   
    #find imex
    imexid = ev.xpath( "./@imexId",
                       namespaces=mifNs )
    if imexid:
        evid['imexid'] = imexid[0]
    else:
        evid['imexid'] = ''
 
    ptList = []
    evid['ptList'] = ptList

    # find participants
    evParts = ev.xpath( ".//m:participantList/m:participant",
                        namespaces=mifNs )

    for part in evParts:
        cpart = {}
        ptList.append(cpart)

        # find interactor
        interactor = None
        intRefId = part.xpath( "./m:interactorRef/text()",
                               namespaces=mifNs )
        if intRefId:  # compact MIF
            interactor = irByRefId[ intRefId[0] ]
        else:         # expanded MIF 
            intId = part.xpath( ".//m:interactor/@id",
                                namespaces=mifNs )
            if intId:
                interactor = irByRefId[ intId[0] ]
        cpart['interactor'] = interactor

        # experimental/biological role
        # (bait, prey, ancillary, enzyme, substrate, etc.)
        partRoleList = []
        cpart['roleList'] = partRoleList
        roleList = part.xpath( ".//*[local-name(.)='experimentalRole' or " +  
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
                partRoleList.append(role)
import json        
print( json.dumps( evidList, sort_keys=True, indent=4) )


    
    
    
