from sequence import fetch_spliced_sequence 

class Codons:
    codon_tbl = {
        "aaa" : "K", "aac" : "N", "aag" : "K", "aat" : "N",
        "aca" : "T", "acc" : "T", "acg" : "T", "act" : "T",
        "aga" : "R", "agc" : "S", "agg" : "R", "agt" : "S",
        "ata" : "I", "atc" : "I", "atg" : "M", "att" : "I",
        "caa" : "Q", "cac" : "H", "cag" : "Q", "cat" : "H",
        "cca" : "P", "ccc" : "P", "ccg" : "P", "cct" : "P",
        "cga" : "R", "cgc" : "R", "cgg" : "R", "cgt" : "R",
        "cta" : "L", "ctc" : "L", "ctg" : "L", "ctt" : "L",
        "gaa" : "E", "gac" : "D", "gag" : "E", "gat" : "D",
        "gca" : "A", "gcc" : "A", "gcg" : "A", "gct" : "A",
        "gga" : "G", "ggc" : "G", "ggg" : "G", "ggt" : "G",
        "gta" : "V", "gtc" : "V", "gtg" : "V", "gtt" : "V",
        "taa" : "*", "tac" : "Y", "tag" : "*", "tat" : "Y",
        "tca" : "S", "tcc" : "S", "tcg" : "S", "tct" : "S",
        "tga" : "*", "tgc" : "C", "tgg" : "W", "tgt" : "C",
        "tta" : "L", "ttc" : "F", "ttg" : "L", "ttt" : "F"  }
    def translate(cdn):
        try:
            return codon_tbl[ cdn.lower() ]
        except:
            return None

def getCodonFromSequence(genomes_rule, map, pos):
    seq = fetch_spliced_sequence(genomes_rule, map)
    for lexon in map["utr5"]:
        len5 += lexon[1]-lexon[0]
    lencds = 0
    for cexon in map["cds"]:
        lencds += cexon[1]-cexon[0]
    if pos < len5:
        return "utr5"
    if pos >= len5+lencds:
        return "utr3"
    ofst = (pos-len5-1)/3
    return seq[ofst,ofst+2]
