
def search_for_refseq(qsymbol,organism,geneius_db):

    query = "select entrez.entrez_id,entrez.type,entrez.official_symbol,"
    query += "entrez.official_gene_name,entrez.other_id,entrez.other_symbols,"
    query += "gref.refseq_rna, species.name from tbl_entrez_xref as entrez "
    query += "left join tbl_species as species on "
    query += "entrez.species=species.tax_id "
    query += "left join tbl_gene_refseq as gref on "
    query += "gref.entrez_id=entrez.entrez_id "
    query += "where (entrez.entrez_id like \"%%%s%%\" " % qsymbol
    query += "or entrez.official_symbol like \"%%%s%%\" or entrez.official_gene_name "  % qsymbol
    query += "like \"%%%s%%\" or entrez.other_id like \"%%%s%%\" or entrez.other_symbols " %(qsymbol,qsymbol)
    query += "like \"%%%s%%\" or entrez.other_gene_names like \"%%%s%%\" " %(qsymbol,qsymbol)
    query += "or gref.refseq_rna like \"%%%s%%\" or gref.refseq_protein like \"%%%s%%\") "%(qsymbol,qsymbol)
    query += "and species.name like \"%%%s%%\" " % organism
    query += "order by entrez.entrez_id;"

    results = []
    last_eid = None
    for entry in geneius_db.query(query):
        if entry[0] == last_eid:
            results[-1]['refseq_ids'].append(entry[6])
        else:
            results.append({
                    "entrez_id":entry[0],
                    "type":entry[1],
                    "official_symbol":entry[2],
                    "official_gene_name":entry[3],
                    "other_id":entry[4],
                    "other_symbols":entry[5],
                    "refseq_ids":[entry[6]],
                    "species":entry[7],
                    })
            last_eid = entry[0]


    #filter results                                                                                                                                                               
    best = []
    good = []
    ok = []
    for result in results:
        if result['official_symbol'].lower() == qsymbol.lower():
            best.append(result)
        elif qsymbol.lower() in [i.lower() for i in result['other_symbols'].split("|")]:
            good.append(result)
        else:
            ok.append(result)
    
    if len(best)>0:
        retval = best
    elif len(good)>0:
        retval = good
    else:
        retval = ok
    return retval

def fetch_annotations(results, geneius_db):
    res_idx = {}
    for result in results:
        res_idx[str(result["entrez_id"])]=result
        result["annotation"] = {}
    query = "select distinct anno.entrez_id,src.source,term.term "
    query += " from tbl_anno as anno " 
    query += " inner join tbl_anno_src as src on anno.src_id = src.src_id "
    query += " inner join tbl_anno_term as term on anno.term_id = term.term_id "
    query += " where anno.entrez_id in(\'%s\') " % "','".join( list(res_idx.keys()) )
    query += " order by anno.entrez_id,src.source;"

    for entry in geneius_db.query(query):
        entrez = str(entry[0])
        if not res_idx[ entrez ]["annotation"].has_key( entry[1] ):
            res_idx[ entrez ]["annotation"][ entry[1] ] = []
        res_idx[ entrez ]["annotation"][ entry[1] ].append( entry[2] )

    return results

def search_by_annotation(qsymbol,organism,geneius_db):

    query = "select distinct entrez.entrez_id,entrez.type,entrez.official_symbol, "
    query += "entrez.official_gene_name,entrez.other_id,entrez.other_symbols,"
    query += "gref.refseq_rna, species.name, "
    query += "source.source, term.term "
    query += "from tbl_entrez_xref as entrez inner "
    query += "join tbl_gene_refseq_new as gref on "
    query += "gref.entrez_id=entrez.entrez_id inner join tbl_species as species on "
    query += "entrez.species=species.tax_id "
    query += " inner join tbl_anno as anno on entrez.entrez_id = anno.entrez_id "
    query += " inner join tbl_anno_src as source on anno.src_id = source.src_id "
    query += " inner join tbl_anno_term as term on anno.term_id = term.term_id "
    query += "where term.term like \"%%%s%%\" " % qsymbol
    query += "and species.name like \"%%%s%%\" and gref.refseq_rna IS NOT NULL " % organism
    query += "order by source.source, term.term, entrez.entrez_id;"

    res_set = {}
    results = []
    last_eid = None
    last_source = None
    last_term = None
    for entry in geneius_db.query(query):
        if not res_set.has_key(entry[8]):
            #add source to hashes
            res_set[entry[8]] = {}
        if not res_set[entry[8]].has_key(entry[9]):
            #add term to source hash
            res_set[entry[8]][entry[9]] = []
            
        if last_source == entry[8] and last_term == entry[9] and entry[0] == last_eid:
            #same source term and gene so just add the refseq
            res_set[entry[8]][entry[9]][-1]['refseq_ids'].append(entry[6])
        else:
            #add new gene
            res_set[entry[8]][entry[9]].append({
                    "entrez_id":entry[0],
                    "type":entry[1],
                    "official_symbol":entry[2],
                    "official_gene_name":entry[3],
                    "other_id":entry[4],
                    "other_symbols":entry[5],
                    "refseq_ids":[entry[6]],
                    "species":entry[7],
                    })
            last_eid = entry[0]
            last_source = entry[8]
            last_term = entry[9]


    #filter results                                                                                                                                                               
    #best = []
    #good = []
    #ok = []
    #for result in results:
    #    if result['official_symbol'].lower() == qsymbol.lower():
    #        best.append(result)
    #    elif qsymbol.lower() in [i.lower() for i in result['other_symbols'].split("|")]:
    #        good.append(result)
    #    else:
    #        ok.append(result)
    #
    #if len(best)>0:
    #    retval = best
    #elif len(good)>0:
    #    retval = good
    #else:
    #    retval = ok

    #retval = fetch_annotations(retval, geneius_db)
    return res_set
