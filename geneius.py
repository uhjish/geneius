#!/usr/bin/env python -W ignore::DeprecationWarning

from gsettings import settings

import sys
import cgi
import MySQLdb
import json

if settings.DEBUG:
    import cgitb
    cgitb.enable()

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
    
class GeneiusDb:
    
    def __init__(self,return_obj):
        self.return_obj = return_obj
        try:
            self.db = MySQLdb.connect(host=settings.DB_SERVER,
                                 user=settings.DB_USER,
                                 passwd=settings.DB_PASSWORD,
                                 db = settings.DB_DATABASE)
        except MySQLdb.OperationalError as oe:
            self.return_obj.error = str(oe)
            print_return(self.return_obj)

        self.cursor = self.db.cursor()

    def query(self,qstring):
        self.cursor.execute(qstring)
        return self.cursor.fetchall()


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


def lookup_refseq(symbols,geneius_db):

    sids = []
    for symbol in symbols:
        sids.append("main.refseq_id=\"%s\"" % symbol)
    
    query = "select * from tbl_refMain as main inner join "
    query += "tbl_refExon as exon on exon.ref_id=main.id "
    query += "where %s order by main.refseq_id,main.id,exon_start ;" % " or ".join(sids)

    results = []
    map_id = None
    mref = None
    for entry in geneius_db.query(query):
        new_exon = {
            "number":entry[10],
            "start":entry[11],
            "end":entry[12]
            }
        new_mapping = {
            "chr":entry[2],
            "strand":entry[3],
            "start":entry[4],
            "end":entry[5],
            "cds_start":entry[6],
            "cds_end":entry[7],
            "num_exons":entry[8],
            "exons":[new_exon],
            }
        #if entry[3] == "+":
        #    new_mapping["utr5_start"]= entry[4]
        #    new_mapping["utr5_end"]=entry[6]
        #    new_mapping["utr3_start"]=entry[7]
        #    new_mapping["utr3_end"]= entry[5]
        #else:
        #    new_mapping["utr5_start"]=entry[7]
        #    new_mapping["utr5_end"]=entry[5]
        #    new_mapping["utr3_start"]=entry[4]
        #    new_mapping["utr3_end"]= entry[6]


        if entry[0] == map_id:
            #if entry[8] == entry[10] and entry[3] == "+":
            #    new_exon["end"] = entry[7]
            #elif int(entry[10]) == 1 and entry[3] == "-":
            #    new_exon["end"] = entry[7]
            results[-1]["mappings"][-1]["exons"].append(new_exon)

        elif entry[1] == mref:
            map_id = entry[0]
            results[-1]["mappings"].append(new_mapping)
        else:
            map_id = entry[0]
            mref = entry[1]
            #move exon start to cds start for first exon
            #new_mapping["exons"][0]["start"] = entry[6]
            results.append({
                    "refseq_id":entry[1],
                    "mappings":[new_mapping]
                    })

    return results

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
    
    
    return results


return_obj = ReturnObj()
geneius_db = GeneiusDb(return_obj)

form = cgi.FieldStorage()
action = get_required_var("action",form,return_obj)
allowable_actions = ["search","lookup"]
if not action in allowable_actions:
    returnobj_error(return_obj,"action must be of %s" % allowable_actions)

if action == "search":
    qsymbol = get_required_var("qsymbol",form,return_obj)
    organism = get_required_var("organism",form,return_obj)
    try:
        dbresults = search_for_refseq(qsymbol,organism,geneius_db)
    except MySQLdb.ProgrammingError as pe:
        returnobj_error(return_obj,str(pe))
   
    return_obj.results = dbresults
elif action=="lookup":
    jsonrefids = get_required_var("refseq_ids",form,return_obj)
    try:
        refids = json.loads(jsonrefids)
    except:
        returnobj_error(return_obj,"problem decoding refseq_ids should be json array")
    try:
        dbresults = lookup_refseq(refids,geneius_db)
    except MySQLdb.ProgrammingError as pe:
        returnobj_error(return_obj,str(pe))
        
    return_obj.results = dbresults

print_return(return_obj)

