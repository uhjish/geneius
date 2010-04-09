
def search_for_refseq(qsymbol,organism,geneius_db):

    query = "select entrez.entrez_id,entrez.type,entrez.official_symbol,"
    query += "entrez.official_gene_name,entrez.other_id,entrez.other_symbols,"
    query += "gref.refseq_rna, species.name from tbl_entrez_xref as entrez inner "
    query += "join tbl_gene_refseq_new as gref on "
    query += "gref.entrez_id=entrez.entrez_id inner join tbl_species as species on "
    query += "entrez.species=species.tax_id "
    query += "where (entrez.entrez_id like \"%%%s%%\" " % qsymbol
    query += "or entrez.official_symbol like \"%%%s%%\" or entrez.official_gene_name "  % qsymbol
    query += "like \"%%%s%%\" or entrez.other_id like \"%%%s%%\" or entrez.other_symbols " %(qsymbol,qsymbol)
    query += "like \"%%%s%%\" or entrez.other_gene_names like \"%%%s%%\" " %(qsymbol,qsymbol)
    query += "or gref.refseq_rna like \"%%%s%%\" or gref.refseq_protein like \"%%%s%%\") "%(qsymbol,qsymbol)
    query += "and species.name like \"%%%s%%\" and gref.refseq_rna IS NOT NULL " % organism
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
        return best
    elif len(good)>0:
        return good

    return ok
