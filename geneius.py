#!/usr/bin/env python
#########################
##Imports
#########################

import sys
import cgi

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
    from libgeneius.search import search_for_refseq
    from libgeneius.whereami import *
    from libgeneius.lookup import lookup_refseq
    from libgeneius.sequence import *
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
    obj.error = error
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
    allowable_actions = ["search","lookup","whereami","sequence"]
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

    elif action=="lookup":
        organism=get_optional_var("organism",form,return_obj)
        if not organism:
            organism="%"
        sequence=get_optional_var("sequence",form,return_obj)
        #with_definition=get_optional_var("with_definition",form,return_obj)
        #if str(with_sequence).upper().startswith("T"):
        #    with_sequence=True
        #else:
        #    with_sequence = False
        #if str(with_definition).upper().startswith("T"):
        #    with_definition=True
        #else:
        #    with_definition=False
        jsonrefids = get_required_var("refseq_ids",form,return_obj)
        try:
            refids = json.loads(jsonrefids)
        except:
            returnobj_error(return_obj,"problem decoding refseq_ids should be json array")
        try:
            dbresults = lookup_refseq(refids,organism,geneius_db)
            if sequence and (sequence == "1" or sequence.upper().startswith("T") ):
                    dbresults = fetch_gene_sequences(genomes_rule, dbresults)
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
except GeneiusError, ge:
    returnobj_error(return_obj,str(ge))

print_return(return_obj)

