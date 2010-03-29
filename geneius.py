#!/usr/bin/env python2.6

from gsettings import settings

import sys
import cgi
import MySQLdb
import simplejson as json
import warnings
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
        except MySQLdb.OperationalError, oe:
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


def whereami_flank(org, chr, pos, direction, geneius_db):
    
    flank_limit = 100
    
    if direction == "left":
        cmp = " < "
        ord = " DESC "
        cmp_to = " main.end "
        entry_cmp = 5
    else:
        cmp = " > "
        ord = " ASC "
        cmp_to = " main.start "
        entry_cmp = 4

    query_pos = "main.chr=\""+chr+"\" and "+cmp_to+cmp+pos+" "
    
    query = "select distinct main.* from tbl_refMain as main "
    query += "inner join tbl_gene_refseq as gref on gref.refseq_rna=main.refseq_id "
    query += "inner join tbl_entrez_xref as entrez on gref.entrez_id=entrez.entrez_id "
    query += "inner join tbl_species as species on entrez.species=species.tax_id "
    query += "where species.name like \"%"+org+"%\" and "
    query += query_pos+" order by "+cmp_to+ord+" LIMIT "+str(flank_limit)

    flank=[]
    closest = -1 
    for entry in geneius_db.query(query):
        if closest != -1 and entry[entry_cmp] != closest:
            break
        closest=entry[entry_cmp]
        feature = {
                "chr":entry[2],
                "start":entry[4],
                "end":entry[5],
                "strand":entry[3],
                "refseq_id":entry[1],
                "type":"transcript",
                "distance":int(entry[entry_cmp])-int(pos)
                }
        flank.append(feature)
    return flank

def whereami_inside(org, chr, pos, geneius_db):

    query_pos = "main.chr=\""+chr+"\" and main.start <= "+pos+" and main.end >="+pos
    query = "select main.*,exon.* from tbl_refMain as main inner join "
    query += "tbl_refExon as exon on exon.ref_id=main.id "
    query += "inner join tbl_gene_refseq_new as gref on gref.refseq_rna=main.refseq_id "
    query += "inner join tbl_entrez_xref as entrez on gref.entrez_id=entrez.entrez_id "
    query += "inner join tbl_species as species on entrez.species=species.tax_id "
    query += "where species.name like \"%"+org+"%\" and "
    query += " %s order by main.refseq_id,main.id,exon_start ;" % query_pos

    last_exon_end=-1
    last_exon_num=-1

    inside=[]
    pos = int(pos)
    for entry in geneius_db.query(query):
        if entry[11] <= pos and entry[12] >=pos:
            inside_exon = {
                "chr":entry[2],
                "start":entry[11],
                "end":entry[12],
                "strand":entry[3],
                "refseq_id":entry[1],
                "type":"exon",
                "exon_num":entry[10]
                }
            inside.append(inside_exon)
        if last_exon_end != -1 and last_exon_end <= pos and entry[12] >=pos:
            inside_intron = {
                "chr":entry[2],
                "start":entry[11],
                "end":entry[12],
                "strand":entry[3],
                "refseq_id":entry[1],
                "type":"intron",
                "flank_exons":[last_exon_num,entry[10]]
                }
            inside.append(inside_intron)
        last_exon_end = entry[12]
        last_exon_num = entry[10]
    return inside

def whereami(organism,chr,pos,geneius_db):
    results= {
            "inside":whereami_inside(organism,chr,pos,geneius_db),
            "left":whereami_flank(organism,chr,pos,"left",geneius_db),
            "right":whereami_flank(organism,chr,pos,"right",geneius_db),
            }
    return results #lookup_refseq("TP53",geneius_db)

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

        if entry[0] == map_id:
            results[-1]["mappings"][-1]["exons"].append(new_exon)

        elif entry[1] == mref:
            map_id = entry[0]
            results[-1]["mappings"].append(new_mapping)
        else:
            map_id = entry[0]
            mref = entry[1]
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

warnings.simplefilter('ignore',DeprecationWarning)

return_obj = ReturnObj()
geneius_db = GeneiusDb(return_obj)

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
   
    return_obj.results = dbresults
elif action=="lookup":
    jsonrefids = get_required_var("refseq_ids",form,return_obj)
    try:
        refids = json.loads(jsonrefids)
    except:
        returnobj_error(return_obj,"problem decoding refseq_ids should be json array")
    try:
        dbresults = lookup_refseq(refids,geneius_db)
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

print_return(return_obj)

