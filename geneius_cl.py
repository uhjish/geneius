#!/usr/bin/env python
#########################
##Imports
#########################

import sys
import traceback
import warnings
from optparse import OptionParser

def Fatal_Error(msg):
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

##########################################
###Done imports
#########################################

####
###Return obj that will be converted to json
####
class ReturnObj:
    def __init__(self):
        self.usage = None
        self.results = None
        self.error = None
    
    def to_json(self):
        obj = {}
        obj['results'] = self.results
        obj['error'] = self.error
        obj['usage'] = self.usage
        return json.dumps(obj)

def print_return(robj):
    ret_str = robj.to_json()
    print ret_str
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

#parse options

usage = "usage: %prog action [options]"
try:
    action = sys.argv[1]
except:
    raise Exception("No action defined! Usage python geneius_cl.py action [options] \n where action in [search, lookup, whereami, sequence, codon, mutate, getmapping, pathway]")

parser = OptionParser(usage=usage)

###if statements for the different Geneius web actions
try:
    geneius_db = GeneiusDb(settings.DB_DATABASE,settings.DB_SERVER,
                           settings.DB_USER,settings.DB_PASSWORD)
    genomes_rule = settings.GENOME_PATH
    #form = cgi.FieldStorage()
    #action = get_required_var("action",form,return_obj)
    allowable_actions = ["search","lookup","whereami","sequence","codon","mutate","getmapping","pathway"]
    if not action in allowable_actions:
        returnobj_error(return_obj,"action must be of %s" % allowable_actions)

    if action == "search":
        parser.add_option("-o", "--organism",
                          action="store", dest="organism",
                          help="the species to search in (ex: sapiens or hg19)")
        parser.add_option("-q", "--query",
                          action="store", dest="qsymbol",
                          help="what to search for -- a string or list of strings depending on the action")
        parser.add_option("-a", "--annotate",
                          action="store_true", dest="annotate", default=False,
                          help="add the pathway annotations for each Entrez ID [False]")
        parser.add_option("-f", "--fullsearch",
                          action="store_true", dest="fullsearch", default=True,
                          help="do a full regex search for the query -- disables the tiered search [True]")
        (opts,args)=parser.parse_args(sys.argv[2:])
        try:
            (opts.qsymbol,opts.organism)
        except:
            parser.print_help()
        try:
            dbresults = search_for_refseq(opts.qsymbol,opts.organism,geneius_db,opts.fullsearch)
            if opts.annotate:
                dbresults = fetch_annotations(dbresults, geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))

        if len(dbresults) < 1:
            returnobj_error(return_obj,"No Matches Found for %s" % qsymbol)
        return_obj.results = dbresults
    elif action == "pathway":
        parser.add_option("-o", "--organism",
                          action="store", dest="organism",
                          help="the species to search in (ex: sapiens or hg19)")
        parser.add_option("-q", "--query",
                          action="store", dest="qsymbol",
                          help="which (partial) pathway names to search for")
        (opts,args)=parser.parse_args(sys.argv[2:])
        try:
            if not opts.qsymbol:
                raise Exception("Argument -q pathway_name is required and missing!")
            if not opts.organism:
                raise Exception("Argument -o organism is required and missing!")
        except:
            parser.print_help()
        try:
            dbresults = search_by_annotation(opts.qsymbol,opts.organism,geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))

        if len(dbresults) < 1:
            returnobj_error(return_obj,"No Matches Found for %s" % opts.qsymbol)
        return_obj.results = dbresults
    elif action=="lookup":
        parser.add_option("-o", "--organism",
                          action="store", dest="organism", default="%",
                          help="the species to search in (ex: sapiens or hg19) [default: all]")
        parser.add_option("-q", "--query",
                          action="store", dest="qsymbol",
                          help="comma-separated list of refseq ids to search for")
        parser.add_option("-s", "--sequence",
                          action="store", dest="sequence", default="",
                          help="comma separated list of sequences to fetch dna,rna,protein")
        (opts,args)=parser.parse_args(sys.argv[2:])
        try:
            if not opts.qsymbol:
                raise Exception("Argument -q refseq_ids is required and missing!")
        except:
            parser.print_help()
            raise Exception("Argument -q refseq_ids is required and missing!")
        try:
            refids = map(lambda x: x.strip(), opts.qsymbol.split(","))
        except:
            parser.print_help()
            returnobj_error(return_obj,"problem decoding refseq_ids should be comma separated list")
        try:
            dbresults = lookup_refseq_with_utrs(refids,opts.organism,geneius_db)
            sequence = map( lambda x: x.lower(), opts.sequence.split(",") )
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
        raise Exception("action getmapping not yet implemented!")
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
        parser.add_option("-b", "--build",
                          action="store", dest="build",
                          help="the build to search in (ex:hg19)")
        parser.add_option("-c", "--chr",
                          action="store", dest="chr",
                          help="chromosome to search (ex: chr1)")
        parser.add_option("-p", "--position",
                          action="store", dest="pos", type="int",
                          help="position on the chromosome")
        parser.add_option("-g", "--gene",
                          action="store_true", dest="as_gene", default=False,
                          help="get each matched item as a full-gene not just a fragment")
        (opts,args)=parser.parse_args(sys.argv[2:])
        try:
            if not opts.build:
                raise Exception("Argument build -b required and missing!")
            if not opts.chr:
                raise Exception("Argument chromosome -c required and missing!")
            if not opts.pos:
                raise Exception("Argument position -p required and missing!")
        except:
            parser.print_help()
            raise
        try:
            if as_gene:
                dbresults = whereami_gene(organism, chr, pos, geneius_db)
            else:
                dbresults = whereami(organism, chr, pos, geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results = dbresults
    elif action == "sequence":
        raise Exception("action sequence not yet implemented!")
        parser.add_option("-b", "--build",
                          action="store", dest="build",
                          help="the build to search in (ex:hg19)")
        parser.add_option("-c", "--chr",
                          action="store", dest="chr",
                          help="chromosome to search (ex: chr1)")
        parser.add_option("-p", "--position",
                          action="store", dest="pos", type="int",
                          help="position on the chromosome")
        parser.add_option("-g", "--gene",
                          action="store_true", dest="as_gene", default=False,
                          help="get each matched item as a full-gene not just a fragment")
        (opts,args)=parser.parse_args(sys.argv[2:])
        try:
            if not opts.build:
                raise Exception("Argument build -b required and missing!")
            if not opts.chr:
                raise Exception("Argument chromosome -c required and missing!")
            if not opts.pos:
                raise Exception("Argument position -p required and missing!")
        except:
            parser.print_help()
            raise
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
        raise Exception("action codon not yet implemented!")
        uid = get_required_var("uid", form, return_obj)
        pos = int(get_required_var("pos", form, return_obj))
        try:
            mapping = get_refseq_by_uid(uid,geneius_db)
            res = getCodonFromSequence(genomes_rule, mapping, pos)
        except Exception, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results=res
    elif action == "mutate":
        raise Exception("action mutate not yet implemented!")
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

