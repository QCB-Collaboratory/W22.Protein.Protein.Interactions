#!/usr/bin/python3

import sys
from lxml import etree
from urllib.request import urlopen

uniprotUrl = "https://www.uniprot.org/uniprot/%%ACC%%.xml"

assert len(sys.argv) == 2
accession = sys.argv[1]

accessionUrl = uniprotUrl.replace("%%ACC%%",accession)

xmlDom = etree.parse( urlopen( accessionUrl ))
print( etree.tostring( xmlDom ).decode() )

