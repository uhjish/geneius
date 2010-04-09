def whereami_flank(org, chr, pos, direction, geneius_db):

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
    query += "where (species.name like \"%"+org+"%\" or species.build like \"%"+org+"%\") and "
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


def whereami_inside(org, chr, pos, geneius_db):

    main_fields = ["main.refseq_id","main.chr","main.strand","main.start","main.end"]
    exon_fields = ["exon.exon_num","exon.exon_start","exon.exon_end"]
    species_fields = ["species.name","species.build"]

    query_pos = "main.chr=\""+chr+"\" and main.start <= "+pos+" and main.end >="+pos

    query = "select "+", ".join( main_fields + exon_fields + species_fields) + " from tbl_refMain as main "
    query += " inner join tbl_refExon as exon on exon.ref_id=main.id "
    query += " inner join tbl_species as species on species.tax_id=main.map_org ";
    query += "where (species.name like \"%"+org+"%\" or species.build like \"%"+org+"%\") and "
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


def whereami(organism,chr,pos,geneius_db):
    results= {
            "inside":whereami_inside(organism,chr,pos,geneius_db),
            "left":whereami_flank(organism,chr,pos,"left",geneius_db),
            "right":whereami_flank(organism,chr,pos,"right",geneius_db),
            }
    return results #lookup_refseq("TP53",geneius_db)
