#!/usr/bin/env python
#########################
##Imports
#########################

import sys
import cgi
import traceback
import warnings

def Fatal_Error(msg):
    print "Content-Type: application/json\n"
    print "{ \"error\":\"%s\"}" % msg
    sys.exit(0)

def missing_package(pkg):
    Fatal_Error("Missing Package: %s" % pkg)

#import third party packages
#and check if they're missing
try:
    import simplejson as json
except:
    missing_package("simplejson")
try:
    import MySQLdb
except:
    missing_package("MySQLdb")

#import geneius libraries
#wrap our imports incase we have a problem
#importing
from libgeneius.error import GeneiusError
try:
    from gsettings import settings
    from libgeneius.mysql import GeneiusDb
    from libgeneius.search import *
    from libgeneius.whereami import *
    from libgeneius.lookup import *
    from libgeneius.sequence import *
    from libgeneius.translate import *
except GeneiusError,ge:
    Fatal_Error(str(ge))

#if debug
if settings.DEBUG:
    import cgitb
    cgitb.enable()
##########################################
###Done imports
#########################################

####
###Return obj that will be converted to json
####
class ReturnObj:
    def __init__(self):
        self.usage = None
        self.results = {}
        self.error = None
    
    def to_json(self):
        obj = {}
        obj['results'] = self.results
        obj['error'] = self.error
        obj['usage'] = self.usage
        return json.dumps(obj)

def print_return(robj):
    print "Content-Type: application/json\n"
    print robj.to_json()
    sys.exit()

def returnobj_error(obj,error):
    obj.error = error + "\n" + traceback.format_exc()
    print_return(obj)

def get_required_var(var,form,return_obj):
    if not form.has_key(var):
        returnobj_error(return_obj,"missing variable %s" % var)
    return form.getvalue(var)

def get_optional_var(var,form,return_obj):
    if not form.has_key(var):
        return False
    return form.getvalue(var)


warnings.simplefilter('ignore',DeprecationWarning)
return_obj = ReturnObj()

###if statements for the different Geneius web actions
try:
    geneius_db = GeneiusDb(settings.DB_DATABASE,settings.DB_SERVER,
                           settings.DB_USER,settings.DB_PASSWORD)
    genomes_rule = settings.GENOME_PATH
    form = cgi.FieldStorage()
    action = get_required_var("action",form,return_obj)
    allowable_actions = ["search","lookup","whereami","sequence","codon","mutate","getmapping","pathway"]
    if not action in allowable_actions:
        returnobj_error(return_obj,"action must be of %s" % allowable_actions)

    if action == "search":
        qsymbol = get_required_var("qsymbol",form,return_obj)
        organism = get_required_var("organism",form,return_obj)
        try:
            dbresults = search_for_refseq(qsymbol,organism,geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))

        if len(dbresults) < 1:
            returnobj_error(return_obj,"No Matches Found for %s" % qsymbol)
        return_obj.results = dbresults
    elif action == "pathway":
        qsymbol = get_required_var("qsymbol",form,return_obj)
        organism = get_required_var("organism",form,return_obj)
        try:
            dbresults = search_by_annotation(qsymbol,organism,geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))

        if len(dbresults) < 1:
            returnobj_error(return_obj,"No Matches Found for %s" % qsymbol)
        return_obj.results = dbresults
    elif action=="lookup":
        organism=get_optional_var("organism",form,return_obj)
        if not organism:
            organism="%"
        sequence=get_optional_var("sequence",form,return_obj)
        jsonrefids = get_required_var("refseq_ids",form,return_obj)
        try:
            refids = json.loads(jsonrefids)
        except:
            returnobj_error(return_obj,"problem decoding refseq_ids should be json array")
        try:
            dbresults = lookup_refseq_with_utrs(refids,organism,geneius_db)
            if sequence:
                sequence = json.loads(sequence)
                sequence = map( lambda x: x.lower(), sequence )
                if "dna" in sequence:
                    dbresults = fetch_dna_for_genes(genomes_rule, dbresults)
                if "rna" in sequence:
                    dbresults = fetch_rna_for_genes(genomes_rule, dbresults)
                if "protein" in sequence:
                    dbresults = fetch_protein_for_genes(genomes_rule, dbresults)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results = dbresults
    elif action=="getmapping":
        uid = get_required_var("uid", form, return_obj)
        sequence=get_optional_var("sequence",form,return_obj)
        try:
            dbresults = get_refseq_by_uid(uid,geneius_db)
            if sequence:
                sequence = json.loads(sequence)
                sequence = map( lambda x: x.lower(), sequence )
                if "dna" in sequence:
                    dbresults["dna"] = fetch_mapping_dna(genomes_rule, dbresults)
                if "rna" in sequence:
                    dbresults["rna"] = fetch_mapping_rna(genomes_rule, dbresults)
                if "protein" in sequence:
                    dbresults["protein"] = fetch_mapping_protein(genomes_rule, dbresults)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results = dbresults
    elif action=="whereami":
        organism=get_required_var("organism",form,return_obj)
        chr=get_required_var("chr",form,return_obj)
        pos=get_required_var("pos",form,return_obj)
        as_gene=get_optional_var("as_gene",form,return_obj)
        try:
        #dbresults = search_for_refseq("TP53","sapiens",geneius_db)
        #dbresults = lookup_refseq(refids,geneius_db)
            if as_gene and (as_gene == "1" or as_gene.upper().startswith("T")):
                dbresults = whereami_gene(organism, chr, pos, geneius_db)
            else:
                dbresults = whereami(organism, chr, pos, geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results = dbresults
    elif action == "sequence":
        build=get_required_var("build", form, return_obj)
        jranges=get_required_var("ranges", form, return_obj)
        spliced=get_optional_var("spliced", form, return_obj)
        try:
            ranges = json.loads(jranges)
        except:
            returnobj_error(return_obj,"problem decoding ranges should be json array of dicts")
        genome = genomes_rule.replace("%",build)
        try:
            if spliced and (spliced == "1" or spliced.upper().startswith("T") ):
                    seq = fetch_spliced_sequence(genome, ranges)
            else:
                seq = fetch_sequences(genome, ranges)
        except Exception, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results = seq
    elif action == "codon":
        uid = get_required_var("uid", form, return_obj)
        pos = int(get_required_var("pos", form, return_obj))
        try:
            mapping = get_refseq_by_uid(uid,geneius_db)
            res = getCodonFromSequence(genomes_rule, mapping, pos)
        except Exception, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results=res
    elif action == "mutate":
        uid = get_required_var("uid", form, return_obj)
        mutalleles = json.loads( get_required_var("alleles", form, return_obj) )
        genomic = get_required_var("genomic", form, return_obj)
        try:
            mapping = get_refseq_by_uid(uid,geneius_db)
            if genomic.upper().startswith("T") or genomic == "1":
                chr = get_required_var("chr", form, return_obj)
                start = int(get_required_var("start", form, return_obj))
                end = int(get_required_var("end",form,return_obj))
                res = getGenomicMutationEffects(genomes_rule, mapping, chr, start, end, mutalleles)
            else:
                pos = int(get_required_var("pos", form, return_obj))
                res = getMutationEffects(genomes_rule, mapping, pos, mutalleles)
        except Exception, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results=res
except GeneiusError, ge:
    returnobj_error(return_obj,str(ge))

print_return(return_obj)

