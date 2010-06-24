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
    from libgeneius import settings
    from libgeneius.mysql import GeneiusDb
    from libgeneius.search import *
    from libgeneius.whereami import *
    from libgeneius.lookup import *
    from libgeneius.sequence import *
    from libgeneius.translate import *
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
    def pathway(self,organism, qsymbol):
        dbresults = search_by_annotation(qsymbol,organism,self.geneius_db)
        return dbresults
    def lookup(self,organism, refids, sequence=[]):
        dbresults = lookup_refseq_with_utrs(refids,organism,self.geneius_db)
        if "dna" in sequence:
            dbresults = fetch_dna_for_genes(self.genomes_rule, dbresults)
        if "rna" in sequence:
            dbresults = fetch_rna_for_genes(self.genomes_rule, dbresults)
        if "protein" in sequence:
            dbresults = fetch_protein_for_genes(self.genomes_rule, dbresults)
        return dbresults
    def whereami(self,build, chr, pos, as_gene=False):
        if as_gene:
            dbresults = whereami_gene(organism, chr, pos, self.geneius_db)
        else:
            dbresults = whereami(organism, chr, pos, self.geneius_db)
        return dbresults
    def get_gene_protein_lookup_table( self, org ):
        return get_gene_protein_lookup_table( org, self.geneius_db )
    def get_symbols_for_refseqs( self, org ):
        return get_symbols_for_refseqs( org, self.geneius_db )
