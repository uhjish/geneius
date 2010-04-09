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
    from suds.client import Client
except:
    missing_package("suds")

#import geneius libraries
from libgeneius.error import GeneiusError
#wrap our imports incase we have a problem
#importing
try:
    from gsettings import settings
    from libgeneius.mysql import GeneiusDb
    from libgeneius.search import search_for_refseq
    from libgeneius.whereami import whereami
    from libgeneius.lookup import lookup_refseq
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
        return("%");    
    return form.getvalue(var)


warnings.simplefilter('ignore',DeprecationWarning)
return_obj = ReturnObj()

###if statements for the different Geneius web actions
try:
    geneius_db = GeneiusDb(settings.DB_DATABASE,settings.DB_SERVER,
                           settings.DB_USER,settings.DB_PASSWORD)

    form = cgi.FieldStorage()
    action = get_required_var("action",form,return_obj)
    allowable_actions = ["search","lookup","whereami"]
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
        jsonrefids = get_required_var("refseq_ids",form,return_obj)
        try:
            refids = json.loads(jsonrefids)
        except:
            returnobj_error(return_obj,"problem decoding refseq_ids should be json array")
        try:
            dbresults = lookup_refseq(refids,organism,geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))
        
        return_obj.results = dbresults

    elif action=="whereami":
        organism=get_required_var("organism",form,return_obj)
        chr=get_required_var("chr",form,return_obj)
        pos=get_required_var("pos",form,return_obj)
        try:
        #dbresults = search_for_refseq("TP53","sapiens",geneius_db)
        #dbresults = lookup_refseq(refids,geneius_db)
            dbresults = whereami(organism, chr, pos, geneius_db)
        except MySQLdb.ProgrammingError, pe:
            returnobj_error(return_obj,str(pe))
        return_obj.results = dbresults
except GeneiusError, ge:
    returnobj_error(return_obj,str(ge))

print_return(return_obj)

