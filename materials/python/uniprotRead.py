#!/usr/bin/python3

import sys
from lxml import etree
from urllib.request import urlopen

uprotUrl = "https://www.uniprot.org/uniprot/%%ACC%%.xml"
uprotNs = { 'up': 'http://uniprot.org/uniprot'}
assert len(sys.argv) == 2
accession = sys.argv[1]

accessionUrl = uprotUrl.replace("%%ACC%%",accession)
uprotTree = etree.parse( urlopen( accessionUrl ))

# find name
uprotNameList = uprotTree.xpath("//up:protein//up:fullName/text()",
                                namespaces=uprotNs)

for name in uprotNameList:
    print( "Name: %s" % ( name))

# find accessions

uprotAccList = uprotTree.xpath("//up:accession/text()",
                               namespaces=uprotNs)
accType = "primary"
for acc in uprotAccList:
    print("Accession: %s (%s)" % ( acc, accType ))
    accType ="secondary"
    

# find gene names
    
uprotGeneNameList = uprotTree.xpath( "//up:gene/up:name",
                                     namespaces=uprotNs)
for gname in uprotGeneNameList:
    gnType = gname.xpath( "./@type" )
    gnName = gname.xpath( "./text()" )
    print("Gene Name: %s  (%s)" % (gnName[0] , gnType[0] )) 

# find sequence

uprotSequenceList = uprotTree.xpath( "//up:sequence",
                                     namespaces=uprotNs)
for seq in uprotSequenceList:
    seqMass = seq.xpath( "./@mass")
    seqVersion = seq.xpath( "./@version")
    seqStr = seq.xpath( "./text()")

    print( "Sequence (ver: %s): %s " % ( seqVersion[0], seqStr[0] ))
    print( "Mass: %6.1fkD" %(int(seqMass[0])/1000)) 
