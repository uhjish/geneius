
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
