import codons
from sequence import fetch_mapping_sequence 

def getCodonFromSequence(genomes_rule, map, pos):
    #raise Exception(str(map))
    seq = fetch_mapping_sequence(genomes_rule, map)
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
        codon =  seq[pos-ofst-1:pos-ofst+2]
        aa = codons.translate(codon)
    base = seq[pos-1]
    result = {  "base":base,
                "codon":codon,
                "aa":aa, 
                "offset":ofst }
    return result 
