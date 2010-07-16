from coding import *

def lookup_refseq(symbols,org,geneius_db):
    '''
    Looks up a list of refseqs
    @param symbols a list of refseq id's
    @param org organism to filter by organism
    @param geneius_db mysql wrapper for genenius
    '''
    sids = []
    for symbol in symbols:
        sids.append("main.refseq_id=\"%s\"" % symbol)

    main_fields = ["main.id","main.refseq_id","main.chr","main.strand","main.start","main.end","main.cds_start","main.cds_end","main.num_exons"]
    exon_fields = ["exon.exon_num","exon.exon_start","exon.exon_end"]
    species_fields = ["species.name","species.build"]
    desc_fields = ["descrip.gid","descrip.description"]
    query = " select "+", ".join(main_fields+exon_fields+species_fields+desc_fields)+" from tbl_refMain as main "
    query += " left join tbl_refExon as exon on exon.ref_id=main.id "
    query += " left join tbl_species as species on species.tax_id=main.map_org ";
    query += " left join tbl_refDesc as descrip on main.refseq_id=descrip.refseq_id ";
    query += " where (species.name like \"%"+org+"%\" or species.build like \"%"+org+"%\") and "
    query += " ( %(refs)s ) " % {'refs':" or ".join(sids)}
    query += " order by main.refseq_id,main.id,exon_start ;"

    #store frequently used indices
    f_uid=0
    f_rfsq=1
    f_chr=2
    f_str=3
    f_start=4
    f_end=5
    f_cdsSt=6
    f_cdsEnd=7
    f_nExons=8
    f_exNum=9
    f_exSt=10
    f_exEnd=11
    f_maporg=12
    f_mapbuild=13
    f_gid=14
    f_desc=15

    results = []
    map_id = None
    mref = None

    lutr=[]
    rutr=[]
    cds=[]
    left=True
    for entry in geneius_db.query(query):
        new_exon = {
            "number":entry[f_exNum],
            "start":entry[f_exSt],
            "end":entry[f_exEnd]
            }
        new_mapping = {
            "uid":entry[f_uid],
            "map_org":entry[f_maporg],
            "map_build":entry[f_mapbuild],
            "chr":entry[f_chr],
            "strand":entry[f_str],
            "start":entry[f_start],
            "end":entry[f_end],
            "cds_start":entry[f_cdsSt],
            "cds_end":entry[f_cdsEnd],
            "num_exons":entry[f_nExons],
            "exons":[new_exon]
            }
        if entry[f_uid] == map_id:
            results[-1]["mappings"][-1]["exons"].append(new_exon)
        else:
            if entry[f_rfsq] == mref:
                map_id = entry[0]
                results[-1]["mappings"].append(new_mapping)
            else:
                map_id = entry[f_uid]
                mref = entry[f_rfsq]
                results.append({
                        "refseq_id":entry[f_rfsq],
                        "gid":entry[f_gid],
                        "desc":entry[f_desc],
                        "mappings":[new_mapping]
                        })
    return results

def get_gene_protein_lookup_table( org, geneius_db ):
    '''
    Get a dict of refseq gene to refseq protein mappings for the given organism 
    @param symbols a list of refseq id's
    @param geneius_db mysql wrapper for genenius
    '''

    query =  " select gref.refseq_rna, gref.refseq_protein from tbl_gene_refseq as gref "
    query += " left join tbl_entrez_xref as entrez on gref.entrez_id = entrez.entrez_id "
    query += " left join tbl_species as species on species.tax_id = entrez.species "
    query += " where (species.name like \"%"+org+"%\" or species.build like \"%"+org+"%\") "
    query += " order by gref.refseq_rna, gref.refseq_protein; "

    ref_map = {}

    for entry in geneius_db.query(query):
        ref_map[ entry[0] ] = entry[1]

    return ref_map

def get_symbols_for_refseqs( org, geneius_db ):
    '''
    Get a dict of refseq gene to refseq protein mappings for the given organism 
    @param symbols a list of refseq id's
    @param geneius_db mysql wrapper for genenius
    '''

    query =  " select gref.refseq_rna, entrez.official_symbol from tbl_gene_refseq as gref "
    query += " left join tbl_entrez_xref as entrez on gref.entrez_id = entrez.entrez_id "
    query += " left join tbl_species as species on species.tax_id = entrez.species "
    query += " where (species.name like \"%"+org+"%\" or species.build like \"%"+org+"%\") "
    query += " order by gref.refseq_rna, gref.refseq_protein; "

    ref_map = {}

    for entry in geneius_db.query(query):
        ref_map[ entry[0] ] = entry[1]

    return ref_map

def get_symbol_for_refseq( refseq, org, geneius_db ):
    '''
    Get a dict of refseq gene to refseq protein mappings for the given organism 
    @param symbols a list of refseq id's
    @param geneius_db mysql wrapper for genenius
    '''

    query =  " select gref.refseq_rna, entrez.official_symbol from tbl_gene_refseq as gref "
    query += " left join tbl_entrez_xref as entrez on gref.entrez_id = entrez.entrez_id "
    query += " left join tbl_species as species on species.tax_id = entrez.species "
    query += " where gref.refseq_rna = \""+refseq+ "\" and " 
    query += " (species.name like \"%"+org+"%\" or species.build like \"%"+org+"%\"); "

    for entry in geneius_db.query(query):
        return entry[1]

def get_symbols_for_refseqs_genomic( org, geneius_db ):
    '''
    Get a dict of refseq gene to refseq protein mappings for the given organism 
    @param symbols a list of refseq id's
    @param geneius_db mysql wrapper for genenius
    '''

    query =  " select distinct gref.refseq_rna, entrez.official_symbol from tbl_gene_refseq as gref "
    query += " left join tbl_entrez_xref as entrez on gref.entrez_id = entrez.entrez_id "
    query += " left join tbl_species as species on species.tax_id = entrez.species "
    query += " left join tbl_refMain as main on main.refseq_id = gref.refseq_rna "
    query += " where (species.name like \"%"+org+"%\" or species.build like \"%"+org+"%\") "
    query += " order by main.chr, main.start, gref.refseq_rna; "

    ref_map = {}

    for entry in geneius_db.query(query):
        ref_map[ entry[0] ] = entry[1]

    return ref_map

def lookup_refseq_with_utrs(symbols,org,geneius_db):
    res = lookup_refseq(symbols, org, geneius_db)
    for gene in res:
       gene = addcoding(gene)
    return res 

def get_refseq_by_uid(uid,geneius_db):
    '''
    Looks up a refseq mapping by uid 
    @param uid unique id from refseq table
    @param geneius_db mysql wrapper for genenius
    '''

    main_fields = ["main.id","main.refseq_id","main.chr","main.strand","main.start","main.end","main.cds_start","main.cds_end","main.num_exons"]
    exon_fields = ["exon.exon_num","exon.exon_start","exon.exon_end"]
    species_fields = ["species.name","species.build"]
    desc_fields = ["descrip.gid","descrip.description"]
    query = " select "+", ".join(main_fields+exon_fields+species_fields+desc_fields)+" from tbl_refMain as main "
    query += " left join tbl_refExon as exon on exon.ref_id=main.id "
    query += " left join tbl_species as species on species.tax_id=main.map_org ";
    query += " left join tbl_refDesc as descrip on main.refseq_id=descrip.refseq_id ";
    query += " where main.id="+uid
    query += " order by main.refseq_id,main.id,exon_start ;"

    #store frequently used indices
    f_uid=0
    f_rfsq=1
    f_chr=2
    f_str=3
    f_start=4
    f_end=5
    f_cdsSt=6
    f_cdsEnd=7
    f_nExons=8
    f_exNum=9
    f_exSt=10
    f_exEnd=11
    f_maporg=12
    f_mapbuild=13
    f_gid=14
    f_desc=15

    results = []
    map_id = None
    mref = None

    lutr=[]
    rutr=[]
    cds=[]
    left=True

    for entry in geneius_db.query(query):
        new_exon = {
            "number":entry[f_exNum],
            "start":entry[f_exSt],
            "end":entry[f_exEnd]
            }
        new_mapping = {
            "uid":entry[f_uid],
            "map_org":entry[f_maporg],
            "map_build":entry[f_mapbuild],
            "chr":entry[f_chr],
            "strand":entry[f_str],
            "start":entry[f_start],
            "end":entry[f_end],
            "cds_start":entry[f_cdsSt],
            "cds_end":entry[f_cdsEnd],
            "num_exons":entry[f_nExons],
            "exons":[new_exon]
            }
        if entry[f_uid] == map_id:
            results[-1]["mappings"][-1]["exons"].append(new_exon)
        else:
            if entry[f_rfsq] == mref:
                map_id = entry[0]
                results[-1]["mappings"].append(new_mapping)
            else:
                map_id = entry[f_uid]
                mref = entry[f_rfsq]
                results.append({
                        "refseq_id":entry[f_rfsq],
                        "gid":entry[f_gid],
                        "desc":entry[f_desc],
                        "mappings":[new_mapping]
                        })
        for gene in results:
           gene = addcoding(gene)
    return results[-1]["mappings"][-1]
