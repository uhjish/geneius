def whereami_flank(build, chr, pos, direction, geneius_db):

    flank_limit = 100

    main_fields = ["main.refseq_id","main.chr","main.strand","main.start","main.end"]
    species_fields = ["species.name","species.build"]

    if direction == "left":
        cmp = " < "
        ord = " DESC "
        cmp_to = " main.end "
        entry_cmp = 4
    else:
        cmp = " > "
        ord = " ASC "
        cmp_to = " main.start "
        entry_cmp = 3

    query_pos = "main.chr=\""+chr+"\" and "+cmp_to+cmp+pos+" "

    query = "select distinct "+ ", ".join(main_fields+species_fields) + " from tbl_refMain as main "
    query += "inner join tbl_species as species on main.map_org=species.tax_id "
    query += "where (species.build = \""+build+"\") and "
    query += query_pos+" order by "+cmp_to+ord+" LIMIT "+str(flank_limit)

    flank=[]
    closest = -1
    for entry in geneius_db.query(query):
        if closest != -1 and entry[entry_cmp] != closest:
            break
        closest=entry[entry_cmp]
        feature = {
                "map_org":entry[5],
                "map_build":entry[6],
                "chr":entry[1],
                "start":entry[3],
                "end":entry[4],
                "strand":entry[2],
                "refseq_id":entry[0],
                "type":"transcript",
                "distance":int(entry[entry_cmp])-int(pos)
                }
        flank.append(feature)
    return flank


def whereami_inside(build, chr, pos, geneius_db):

    main_fields = ["main.refseq_id","main.chr","main.strand","main.start","main.end"]
    exon_fields = ["exon.exon_num","exon.exon_start","exon.exon_end"]
    species_fields = ["species.name","species.build"]

    query_pos = "main.chr=\""+chr+"\" and main.start <= "+pos+" and main.end >="+pos

    query = "select "+", ".join( main_fields + exon_fields + species_fields) + " from tbl_refMain as main "
    query += " inner join tbl_refExon as exon on exon.ref_id=main.id "
    query += " inner join tbl_species as species on species.tax_id=main.map_org ";
    query += "where (species.build = \""+build+"\") and "
    query += " %s order by main.refseq_id,main.id,exon_start ;" % query_pos

    last_exon_end=-1
    last_exon_num=-1

    inside=[]
    pos = int(pos)
    for entry in geneius_db.query(query):
        if entry[6] <= pos and entry[7] >=pos:
            inside_exon = {
                "chr":entry[1],
                "start":entry[6],
                "end":entry[7],
                "strand":entry[2],
                "refseq_id":entry[0],
                "type":"exon",
                "exon_num":entry[5],
                "map_org":entry[8],
                "map_build":entry[9]
                }
            inside.append(inside_exon)
        if last_exon_end != -1 and last_exon_end <= pos and entry[6] >=pos:
            inside_intron = {
                "chr":entry[1],
                "start":last_exon_end,
                "end":entry[6],
                "strand":entry[2],
                "refseq_id":entry[0],
                "type":"intron",
                "flank_exons":[last_exon_num,entry[5]],
                "map_org":entry[8],
                "map_build":entry[9]
                }
            inside.append(inside_intron)
        last_exon_end = entry[7]
        last_exon_num = entry[5]
    return inside


def whereami(build,chr,pos,geneius_db):
    results= {
            "inside":whereami_inside(build,chr,pos,geneius_db),
            "left":whereami_flank(build,chr,pos,"left",geneius_db),
            "right":whereami_flank(build,chr,pos,"right",geneius_db),
            }
    return results #lookup_refseq("TP53",geneius_db)

def whereami_gene(build,chr,pos,geneius_db):
    
    main_fields = ["main.id","main.refseq_id","main.chr","main.strand","main.start","main.end","main.cds_start","main.cds_end","main.num_exons"]
    exon_fields = ["exon.exon_num","exon.exon_start","exon.exon_end"]
    species_fields = ["species.name","species.build"]
    
    query_pos = "main.chr=\""+chr+"\" and main.start <= "+pos+" and main.end >="+pos
    
    query = "select "+", ".join( main_fields + exon_fields + species_fields) + " from tbl_refMain as main "
    query += " inner join tbl_refExon as exon on exon.ref_id=main.id "
    query += " inner join tbl_species as species on species.tax_id=main.map_org ";
    query += "where (species.build = \""+build+"\") and "
    query += " %s order by main.refseq_id,main.id,exon_start ;" % query_pos

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
