from lxml import etree

xml = '''<mif:protein xmlns:mif="http://psi.hupo.org/mi/mif" acc="P60010">
         <mif:seq>MKYDDEW...</mif:seq>
      </mif:protein>'''

xmlDom = etree.fromstring( xml )
e = xmlDom.xpath('/m:protein/m:seq',
                 namespaces={'m': 'http://psi.hupo.org/mi/mif'})

print( e[0].tag.decode() )

qname = etree.QName(e[0])
print( qname.localname.decode() )
print( qname.namespace.decode() ) 
