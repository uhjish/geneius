
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
    query = " select "+", ".join(main_fields+exon_fields+species_fields)+" from tbl_refMain as main "
    query += " inner join tbl_refExon as exon on exon.ref_id=main.id "
    query += " inner join tbl_species as species on species.tax_id=main.map_org ";
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
                        "mappings":[new_mapping]
                        })

    return results

