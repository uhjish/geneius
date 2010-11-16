#!/usr/bin/env python
#########################
##Imports
#########################

import sys
import traceback
import warnings

def Fatal_Error(msg):
    print "{ \"error\":\"%s\"}" % msg
    sys.exit(0)

def missing_package(pkg):
    Fatal_Error("Missing Package: %s" % pkg)

#import third party packages
#and check if they're missing
try:
    import MySQLdb
except:
    missing_package("MySQLdb")

#import geneius libraries
#wrap our imports incase we have a problem
#importing
from libgeneius.error import GeneiusError
try:
    import settings
    from mysql import GeneiusDb
    from search import *
    from whereami import *
    from lookup import *
    from sequence import *
    from translate import *
    from rightwords import *
except GeneiusError,ge:
    Fatal_Error(str(ge))

##########################################
###Done imports
#########################################
warnings.simplefilter('ignore',DeprecationWarning)

#parse options

class SimpleGeneius:
    def __init__(self):
        try:
            self.geneius_db = GeneiusDb(settings.DB_DATABASE,settings.DB_SERVER,
                                   settings.DB_USER,settings.DB_PASSWORD)
            self.genomes_rule = settings.GENOME_PATH
        except:
            raise
    def search(self,organism, qsymbol, annotate=False, fullsearch=True):
        dbresults = search_for_refseq(qsymbol,organism,self.geneius_db,fullsearch)
        if annotate:
            dbresults = fetch_annotations(dbresults, self.geneius_db)
        return dbresults
    def get_all_refseqs_for_symbol(self,organism, qsymbol):
        dbresults = search_for_refseq(qsymbol,organism,self.geneius_db,False)
        refseqs = []
        for res in dbresults:
            refseqs.extend(res["refseq_ids"])
        return refseqs
    def get_symbol_for_refseq(self,organism, refseq):
        dbresults = search_for_refseq(refseq,organism,self.geneius_db,False)
        if len(dbresults) > 1:
            raise Exception("Why does this refseq have more than one symbol?!?!?!");
        return dbresults[0]["official_symbol"];
    def pathway(self,organism, qsymbol):
        dbresults = search_by_annotation(qsymbol,organism,self.geneius_db)
        return dbresults
    def annotate_genes(self,qsymbols, organism, id_type, simplify=False):
        dbresults = fetch_annotations_by_symbols(qsymbols,organism,id_type,simplify,self.geneius_db)
        return dbresults
    def fetch_annotations_by_symbols(self, qsymbols, organism, id_type, simplify=True):
        qannos = fetch_annotations_by_symbols( qsymbols, organism, id_type, simplify, self.geneius_db )
        return qannos
    def mine_annotations(self, qset, bset, organism, id_type ):
        return mine_annotations(qset, bset, organism, id_type, self.geneius_db)
    def lookup(self, organism, refids, sequence=[]):
        dbresults = lookup_refseq_with_utrs(refids,organism,self.geneius_db)
        if "dna" in sequence:
            dbresults = fetch_dna_for_genes(self.genomes_rule, dbresults)
        if "rna" in sequence:
            dbresults = fetch_rna_for_genes(self.genomes_rule, dbresults)
        if "protein" in sequence:
            dbresults = fetch_protein_for_genes(self.genomes_rule, dbresults)
        return dbresults
    def whereami(self,organism, chr, pos, as_gene=False):
        if as_gene:
            dbresults = whereami_gene(organism, chr, pos, self.geneius_db)
        else:
            dbresults = whereami(organism, chr, pos, self.geneius_db)
        return dbresults
    def get_genes_in_region(self,organism, chr, start, end):
        dbresults = get_genes_in_region(organism, chr, start, end, self.geneius_db)
        return dbresults
    def get_refseq_by_uid(self,uid,sequence=[]):
        dbresults = get_refseq_by_uid(uid, self.geneius_db)
        sequence = map( lambda x: x.lower(), sequence )
        if "dna" in sequence:
            dbresults["dna"] = fetch_mapping_dna(self.genomes_rule, dbresults)
        if "rna" in sequence:
            dbresults["rna"] = fetch_mapping_rna(self.genomes_rule, dbresults)
        if "protein" in sequence:
            dbresults["protein"] = fetch_mapping_protein(self.genomes_rule, dbresults)
        return dbresults
    def get_all_mappings_for_organism(self,org,sequence=[]):
        dbresults = get_all_mappings_for_organism(org, self.geneius_db)
        sequence = map( lambda x: x.lower(), sequence )
        if "dna" in sequence:
            dbresults = fetch_dna_for_genes(self.genomes_rule, dbresults)
        if "rna" in sequence:
            dbresults = fetch_rna_for_genes(self.genomes_rule, dbresults)
        if "protein" in sequence:
            dbresults = fetch_protein_for_genes(self.genomes_rule, dbresults)
        res_mat = {}
        for gene in dbresults:
            for mapping in gene["mappings"]:
                res_mat[mapping["uid"]] = mapping
        return res_mat
    def get_gene_protein_lookup_table( self, org ):
        return get_gene_protein_lookup_table( org, self.geneius_db )
    def get_refseq_uniprot_lookup_table( self, org ):
        return get_refseq_uniprot_lookup_table( org, self.geneius_db )
    def get_symbols_for_refseqs( self, org ):
        return get_symbols_for_refseqs( org, self.geneius_db )
    def get_symbols_for_entrez(self, org = ""):
        return get_symbols_for_entrez(self.geneius_db, org)
    def get_descriptions_for_official_symbols(self, org=""):
        return get_descriptions_for_official_symbols(self.geneius_db, org)
    def get_synonyms_for_official_symbols(self, org = ""):
        return get_synonyms_for_official_symbols(self.geneius_db, org)
    def get_symbol_for_refseq( self, refseq, org ):
        return get_symbol_for_refseq( refseq, org, self.geneius_db )
    def get_symbols_for_refseqs_genomic( self, org ):
        return get_symbols_for_refseqs_genomic( org, self.geneius_db )
