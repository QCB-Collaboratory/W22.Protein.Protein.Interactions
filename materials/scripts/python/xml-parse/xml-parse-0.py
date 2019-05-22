import sys

from lxml import etree as ET
import json

#test_mif25='ftp://ftp.ebi.ac.uk/pub/databases/intact/current/psi25/pmid/2019/15138291.xml'
test_mif25="/home/lukasz/AAtWork/DIP/QCB/PPI2019/files/mif25/11034202.xml"

source = test_mif25
#source = sys.argv[1]

parser = ET.XMLParser( remove_blank_text=True ) 
dom = ET.parse( source, parser )
     
#<entrySet level="2" version="5" minorVersion="4" 
#          xsi:schemaLocation="http://psi.hupo.org/mi/mif http://psidev.sourceforge.net/mi/rel25/src/MIF254.xsd" 
#          xmlns="http://psi.hupo.org/mi/mif" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance">
#

# equvalent to:
#<abc:entrySet abc:xmlns="http://psi.hupo.org/mi/mif" ....

ns = { 'mif': 'http://psi.hupo.org/mi/mif' }

elist = dom.xpath( '//mif:entry', namespaces = ns)

entry = elist[0] # only one entry per file in IMEx primary records data

print( ET.tostring(entry) )  # this is, technically, an array of bytes. convert to string by calling .decode("UTF-8")

