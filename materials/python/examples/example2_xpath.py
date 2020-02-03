from lxml import etree
xml = '<protein acc="P60010"><seq format="fasta">MKYDDEW...</seq></protein>'

xmlDom = etree.fromstring( xml )

root = xmlDom.xpath('/protein')

for child in root[0]:
    print( child.tag.decode() )
    print( child.get("format").decode() )
    print( child.text.decode() )
    print( etree.tostring(child).decode() )
        

t1 = xmlDom.xpath('/protein/seq/text()')
print(t1[0].decode())

t2 = xmlDom.xpath('//seq/text()')
print(t2[0].decode())

e3 = xmlDom.xpath('//protein[@acc="P60010"]/seq')
print(etree.tostring(e3[0]).decode())

e4 = e3[0].xpath('./text()')
print(e4[0].decode())

