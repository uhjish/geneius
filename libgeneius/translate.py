import codons
from sequence import fetch_mapping_rna 
from sequence import fetch_sequence 

leftSplice = { "+": "ss5", "-": "ss3" }
rightSplice = { "+": "ss3", "-": "ss5" }
leftGenic = {"+": "upstream", "-": "downstream"}
rightGenic = {"+": "downstream", "-": "upstream"}



def getGenomicMutationEffects(genomes_rule, map, chr, start, end, newbases):
    genome = genomes_rule.replace("%", map["map_build"])
    allele = fetch_sequence(genome, chr, start, end, map["strand"])
    result = {  "mapping_uid": map["uid"],
                "pos":None,
                "allele":allele,
                "codon":None,
                "aa":None, 
		        "aa_pos":None,
                "offset":None,
             }
    mutations = []
    effect = []
    if chr != map["chr"]:
        effect=["off-chromosome"]
    elif end <= map["start"]:
        effect=[leftGenic[map["strand"]]]
    elif start >= map["end"]:
        effect=[rightGenic[map["strand"]]]
    elif end - start > 1:
        effect=["complex"]
    else:
        #check exons
        mrna_pos = 0
        for exon in map["exons"]:
            if start >= exon["start"] and end <= exon["end"]:
                mrna_pos += end - exon["start"]
                if map["strand"] == "-":
                    mrna_pos = map["mrna_length"]-mrna_pos+1
                return getMutationEffects(genomes_rule, map, mrna_pos, newbases)
            elif end == exon["start"]:
                effect=leftSplice[map["strand"]]
                break
            elif start == exon["end"]:
                effect=rightSplice[map["strand"]]
                break
            else:
                mrna_pos += exon["end"] - exon["start"]
    if not effect:
        #intronic
        effect = ["intronic"]
    for nb in newbases:
        mutations.append( { "allele" : nb,
                            "codon" : None,
                            "aa" : None,
                            "effect" : effect } )
    result["mutations"]=mutations
    return result

def getMutationEffects(genomes_rule, map, pos, newbases):
    original = getCodonFromSequence(genomes_rule, map, pos)
    mutations = []
    for newbase in newbases:
        ncod = None
        naa = None
        if newbase == "-":
            if original["codon"].startswith("utr"):
                effect = ["deletion","untranslated"]
            else:
                effect = ["deletion", "frameshift"]
        elif len(newbase) > 1 and (len(newbase)-1) % 3 != 0:
            if original["codon"].startswith("utr"):
                effect = ["insertion","untranslated"]
            else:
                effect = ["insertion", "frameshift"]
        elif len(newbase) > 1 and (len(newbase)-1) % 3 == 0:
            if original["codon"].startswith("utr"):
                effect = ["insertion","untranslated"]
            else:
                effect = ["insertion", "inframe"]
        elif newbase == original["allele"]:
            effect = ["reference"]
            ncod = original["codon"]
            naa = original["aa"]
        elif len(newbase) == 1:
            effect = []
            ofs = original["offset"]
            ocod = original["codon"]
            if ocod.startswith("utr"):
                effect.append("untranslated")
                ncod=ocod
                naa = original["aa"]
            else:
                ncod = ocod[0:ofs] + newbase + ocod[ofs+1:3]
                naa = codons.translate(ncod)
                if naa == original["aa"]:
                    effect.append("synonymous")
                elif original["aa"] == "*":
                    effect.append("runon")
                elif naa == "*":
                    effect.append("nonsense")
                else:
                    effect.append("missense")
        mutations.append( { "allele" : newbase,
                            "codon" : ncod,
                            "aa" : naa,
                            "effect" : effect } )
    original["mutations"] = mutations
    return original


def getCodonFromSequence(genomes_rule, map, pos):
    seq = fetch_mapping_rna(genomes_rule, map)
    if pos < 1 or pos > len(seq):
        raise Exception( "translate.py:getCodonFromSequence - pos out of range [1,len(sequence)]")
    ofst=None
    cod_num = None
    len5=0
    for lexon in map["utr5"]:
        len5 += lexon[1]-lexon[0]
    lencds = 0
    for cexon in map["cds"]:
        lencds += cexon[1]-cexon[0]
    if pos <= len5:
        codon = "utr5"
        aa= None
    elif pos > len5+lencds:
        codon= "utr3"
        aa=None
    else:
        ofst = (pos-len5-1) % 3
	cod_num = (pos-len5-1) / 3 + 1 
        codon =  seq[pos-ofst-1:pos-ofst+2]
        aa = codons.translate(codon)
    base = seq[pos-1]
    result = {  "mapping_uid": map["uid"],
                "pos":pos,
                "allele":base,
                "codon":codon,
                "aa":aa, 
		        "aa_pos":cod_num,
                "offset":ofst }
    return result
