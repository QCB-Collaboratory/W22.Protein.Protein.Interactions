import sys
import json

from lxml import etree as ET

class Parser():

    def __init__( self ):
        self.ns = { 'mif': 'http://psi.hupo.org/mi/mif' }
        
    def parse( self, source ):

        mif = {}
        #mif["entry"] = {}
        mif["source"] = {}        # source 
        mif["i10r"] = {}          # interactors
        mif["e10t"] = {}          # experiments 
        mif["i11n"] = {}          # interactions
        
        parser = ET.XMLParser( remove_blank_text=True ) 
        dom = ET.parse( source, parser )
        
        elist = dom.xpath( '//mif:entry',
                           namespaces = self.ns)
        
        entry = elist[0]

        # parse source
        #-------------
        
        mif["source"] = self.__parseSource( entry )

        # parse experiments
        #------------------

        expl = entry.xpath( ".//mif:experimentDescription",
                            namespaces = self.ns )

        for cexp in expl:
        
            e10t = self.__parseExperiment( cexp )
            mif["e10t"][ e10t['rid'] ] = e10t              
      
        # parse interactors 
        #------------------
            
        i10rl = entry.xpath( ".//mif:interactor",
                             namespaces = self.ns )
            
        for cintr in i10rl:
            ii = self.__parseInteractor( cintr )
            mif["i10r"][ii['rid']] = ii
         
        # parse interactions
        #-------------------
                
        i11nl = entry.xpath( ".//mif:interaction",
                             namespaces = self.ns )
        
        for cintn in i11nl:
            ii = self.__parseInteraction( mif, cintn )
            mif["i11n"][ii['rid'] ] = ii

        return mif

    # private methods
    #---------------------------------------------------------------------------            

    def __parseInteraction( self, mif, dom ):

        i11n = {}
        
        refID =  dom.xpath( "./@id",
                            namespaces = self.ns )[0]        

        i11n['rid'] = refID


        imexID =  dom.xpath( "./@imexId",
                             namespaces = self.ns )

        if imexID:
            i11n['imex'] = imexID[0]
        else:
            i11n['imex'] = ""



        # participant list
        #-----------------

        p11tl = dom.xpath( "./mif:participantList/mif:participant",
                           namespaces = self.ns )
        p11t = []
        for p in p11tl:
            p11t.append( self.__parseParticipant( mif, p ) )
            
        i11n['p11t'] = p11t

        # experiment ref list
        #--------------------

        e10tl = dom.xpath( "./mif:experimentList",
                           namespaces = self.ns )[0]

        e10t = []
    
        for e in e10tl:

            if e.xpath('local-name(.)') ==  'experimentRef':
                e10t.append( mif["e10t"][ e.xpath('./text()')[0] ])
            else:
                e10t.append( self.__parseExperiment( e ) )

        i11n["e10t"] = e10t


        # interaction type
        #-----------------
        
        itype = self.__parseCV( dom.xpath( ".//mif:interactionType",
                                           namespaces = self.ns )[0] )
                                
        i11n["itype"] = itype
    
        return i11n
            

    def __parseParticipant( self, mif, dom ):    

        p11t = {}

        #print(ET.tostring(dom))

        refID =  dom.xpath( "./@id",
                            namespaces = self.ns )[0]
        
        # interactor
        #-----------
        
        iref = dom.xpath( "./mif:interactorRef/text()",
                          namespaces = self.ns )

        if len( iref ) == 1:
            iref = iref[0]
            p11t['iref'] = iref
            p11t['i10r'] = mif["i10r"][iref]
        else:
            p11t['iref'] = "N/A"
            
        if p11t['iref'] == "N/A":
            
            i10rl = entry.xpath( ".//mif:interactor",
                                 namespaces = self.ns )
            if i10rl:
                p11t['i10r'] =  self.__parseInteractor( i10rl[0] )

        # bio role
        #---------
        
        bioRoleDom = dom.xpath( "./mif:biologicalRole",
                                namespaces = self.ns )
        p11t['brole'] = self.__parseCV( bioRoleDom[0] ) 

    
        # exp role
        #---------
        
        erole = dom.xpath( ".//mif:experimentalRole",
                        namespaces = self.ns )
        
        if len(erole) > 0:
            p11t['erole'] = self.__parseCV( erole[0] )
            
         
        # exp prep role
        #--------------
        
        eprep = dom.xpath( ".//mif:experimentalPreparation",
                           namespaces = self.ns )
        if len( eprep ) > 0:
            p11t['eprep'] = self.__parseCV( eprep[0] )

        
        # features
        #---------
        
        ftrDom = dom.xpath( ".//mif:feature",
                            namespaces = self.ns )
        
        p11t['feature'] = self.__parseFeature( ftrDom )
        
        return p11t



    def __parseFeature( self, dom ):

        ftrl = []
    
        for f in dom:

            ftr = {}
            
            ftrSL = f.xpath( "./mif:names/mif:shortLabel/text()",
                             namespaces = self.ns )
            if ftrSL:
                ftr['label'] = ftrSL[0]
            else:
                ftr['label'] = ""
                
            ftrFN =  f.xpath( "./mif:names/mif:fullName/text()",
                               namespaces = self.ns )
            if ftrFN:
                ftr['name'] = ftrFN[0] 
            else:
                ftr['name'] = ""
                
            # feature type
            
            ftrTPDom = f.xpath( "./mif:featureType",
                                  namespaces = self.ns )
            
            ftr['type'] = self.__parseCV( ftrTPDom[0] )

            # feature det method
            
            ftrDMDom = f.xpath( "./mif:featureDetectionMethod",
                                  namespaces = self.ns )
            if len( ftrDMDom ) > 0:
                ftr['detmeth'] = self.__parseCV( ftrDMDom[0] )
            else:
                ftr['detmeth'] = {}
            # feature range(s)

            ftrRGDom = f.xpath( ".//mif:featureRange",
                                  namespaces = self.ns )
                        
            ftr['range'] = self.__parseRange( ftrRGDom )
            ftrl.append(ftr)

        return ftrl
                  
    def __parseRange( self, dom ):

        rng = []

        for r in dom:
            
            crng = {}

            sSTDom = r.xpath( "./mif:startStatus",
                                namespaces = self.ns )
            
            sST = self.__parseCV(  sSTDom[0])
            sPS = r.xpath( "./mif:begin/@position",
                              namespaces = self.ns )

            crng["begStatus"] = sST
            crng["beg"] = sPS

            eSTDom = r.xpath( "./mif:endStatus",
                                namespaces = self.ns )
            
            eST = self.__parseCV( eSTDom[0])
            ePS = r.xpath( "./mif:end/@position",
                             namespaces = self.ns )

            crng["endStat"] = sST
            crng["end"] = sPS
            
            
            lnk =  r.xpath( "./mif:isLink/text()",
                            namespaces = self.ns )
            
            if lnk:
                crng["link"] = True

            rng.append(crng)

        return rng
                

    def __parseInteractor( self, dom ):

        idix = {}
        i10r = {}
        
        i10r['rid'] =  dom.xpath( "./@id",
                                  namespaces = self.ns )[0]

    
        SL =  dom.xpath( "./mif:names/mif:shortLabel/text()",
                          namespaces = self.ns )
        if len(SL) > 0:
            i10r[ 'label' ] = SL[0]
        else:
            i10r[ 'label' ] = ''

            
        FN =  dom.xpath( "./mif:names/mif:fullName/text()",
                         namespaces = self.ns )
        if len( FN ) > 0:
            i10r[ 'name' ] = FN[0]
        else:
            i10r[ 'name' ] = ''
            
        i10r[ 'attr' ] =  dom.xpath( "./mif:names/mif:alias/text()",
                                     namespaces = self.ns )
                
        prXRefDom = dom.xpath( "./mif:xref/mif:primaryRef",
                               namespaces = self.ns )

        i10r[ 'pxref' ] = self.__parseRefList( prXRefDom )
        
        scXRefDom = dom.xpath( "./mif:xref/mif:secondaryRef",
                               namespaces = self.ns )
        i10r[ 'sxref' ] = self.__parseRefList( scXRefDom )
        
        itpDom= dom.xpath( "./mif:interactorType",
                                     namespaces = self.ns )
    
        i10r[ 'type' ] =  self.__parseCV( itpDom[0] )
        
        return i10r
        
    def __parseRefList( self, ref ):
        ret = []
        for r in ref:
            acID = r.xpath( "./@id",
                              namespaces = self.ns )
            ns = r.xpath( ".//@db",
                              namespaces = self.ns )
            nsID = r.xpath( ".//@dbAc",
                              namespaces = self.ns )
            verID  = r.xpath( "@version",
                              namespaces = self.ns )
            ref = {}
            if acID:
                ref['ac'] = acID[0]
            if ns:
                ref['ns'] = ns[0]
            if nsID:
                ref['nsAc'] = nsID[0]
            if verID:
                ref['ver'] = verID[0]

            ret.append( ref )
            
        return ret

    
    def __parseCV( self, cvdom ):

        cvID = cvdom.xpath( "./mif:xref/mif:primaryRef/@id",
                            namespaces = self.ns )
        cvDB = cvdom.xpath( "./mif:xref/mif:primaryRef/@db",
                            namespaces = self.ns )
        cvDBID = cvdom.xpath( "./mif:xref/mif:primaryRef/@dbAc",
                              namespaces = self.ns )
        cvSL = cvdom.xpath( "./mif:names/mif:shortLabel/text()",
                            namespaces = self.ns )
        cvFN = cvdom.xpath( "./mif:names/mif:fullName/text()",
                            namespaces = self.ns )
    
        return {"ns":cvDB[0],"ac":cvID[0], "label":cvSL[0],"name":cvSL[0],"nsAc":cvDBID[0] }


    def __parseSource( self, entry ):

        source = {}
        
        # source
        #-------

        msrc = entry.xpath( "./mif:source",
                            namespaces = self.ns )[0]
        
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
                                   namespaces = self.ns )[0]

        name =  msrc.xpath( ".//mif:fullName/text()",
                                   namespaces = self.ns )
        if name:
            source['name'] = name[0]
        else:
            source['name'] = source['label']

        source['ns'] = msrc.xpath( ".//mif:primaryRef/@db",
                                   namespaces = self.ns )[0]
        source['ac'] = msrc.xpath( ".//mif:primaryRef/@id",
                                   namespaces = self.ns )[0]
                
        return source



    def __parseExperiment( self, dom ):


        e10t = {}
            
        # <experimentDescription id="11175">
        
        e10t['rid'] =  dom.xpath( "./@id",
                                  namespaces = self.ns )[0]
        
        self.crid = e10t['rid']
            
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
                                  namespaces = self.ns )[0]
            
        #     <xref>
        #         <primaryRef db="mint" dbAc="MI:0471" id="MINT-722978" refType="identity" refTypeAc="MI:0356"/>
        #         <secondaryRef db="pubmed" dbAc="MI:0446" id="12437929" refType="secondary-ac" refTypeAc="MI:0360"/>
        #         <secondaryRef db="doi" dbAc="MI:0574" id="10.1016/S0092-8674(02)01076-0" refType="secondary-ac" refTypeAc="MI:0360"/>
        #         <secondaryRef db="intact" dbAc="MI:0469" id="EBI-8558582" refType="identity" refTypeAc="MI:0356"/>
        #         <secondaryRef db="mint" dbAc="MI:0471" id="MINT-5215918" refType="primary-reference" refTypeAc="MI:0358"/>
        #     </xref>
        
            
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
                                        namespaces = self.ns )[0]
            
        e10t['host']['label'] = dom.xpath( ".//mif:hostOrganism//mif:shortLabel/text()",
                                           namespaces = self.ns )[0]

        name = dom.xpath( ".//mif:hostOrganism//mif:fullName/text()",
                          namespaces = self.ns ) 
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
                                 namespaces = self.ns )
            
        e10t['intMth'] = self.__parseCV( intIdMthDom[0] )
                        
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
                                 namespaces = self.ns )
        
        e10t['prtMth'] = self.__parseCV( prtIdMthDom[0] )            
         
        return e10t



