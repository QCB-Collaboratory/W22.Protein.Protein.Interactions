from lxml import etree
xml = '<protein acc="P60010"><seq format="fasta">MKYDDEW...</seq></protein>'

dom = etree.fromstring( xml )

for child in dom:
    print( child.tag.decode() )
    print( child.get("format").decode() )
    print( child.text.decode() )
    print( etree.tostring(child).decode() )
        

