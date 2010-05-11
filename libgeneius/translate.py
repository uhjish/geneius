import codons
from sequence import fetch_mapping_rna 

def getMutationEffects(genomes_rule, map, pos, newbases):
    original = getCodonFromSequence(genomes_rule, map, pos)
    mutations = []
    for newbase in newbases:
        ncod = None
        naa = None
        if newbase == "-":
            effect = ["deletion", "frameshift"]
        elif len(newbase) > 1 and len(newbase) % 3 != 0:
            effect = ["insertion", "frameshift"]
        elif len(newbase) % 3 == 0:
            effect = ["insertion", "in frame"]
        elif newbase == original["allele"]:
            effect = ["reference"]
            ncod = original["codon"]
            naa = original["aa"]
        elif len(newbase) == 1:
            effect = ["SNP"]
            ofs = original["offset"]
            ocod = original["codon"]
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
	cod_num = (pos-len5-1) / 3 
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
