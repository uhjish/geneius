import math
from warnings import warn
from stats import *
from mysql import GeneiusDb
from search import fetch_annotations_by_symbols

def count_annotations(gene2annos):
    anno2genes = {}
    for (gene, annos) in gene2annos.items():
        for anno in annos:
            if not anno2genes.has_key(anno):
                anno2genes[anno] = set()
            anno2genes[anno].add(gene)
    return anno2genes

def get_significance( qct, bct, qsize, bsize):
    pval = prob(qct, float(bct)/float(bsize), qsize)
    return pval

def mine_annotations( qs, bs, organism, id_type, geneius_db ):
    qset = set()
    bset = set()
    for item in qs:
        qset.add(item)
        bset.add(item)
    for item in bs:
        bset.add(item)
    qsize = len(qset)
    bsize = len(bset)
    qannos = fetch_annotations_by_symbols( qset, organism, id_type, True, geneius_db )
    bannos = fetch_annotations_by_symbols( bset, organism, id_type, True, geneius_db )

    qannos = count_annotations(qannos)
    bannos = count_annotations(bannos)

    res = []

    for (anno, genes) in qannos.items():
        qct = len(genes)
        bct = len( bannos[anno] )
        pval = get_significance( qct, bct, qsize, bsize )
        res.append( {   "anno": anno,
                        "qct": qct,
                        "bct": bct,
                        "pval": pval,
                        "items": genes } )

    result = {  "organism": organism,
                "id_type": id_type,
                "query_size": qsize,
                "background_size": bsize,
                "scores": res }

    return result
 
        
   
