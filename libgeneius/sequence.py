import codons
from libgeneius.error import GeneiusError
try:
    import bx.seq.twobit
except:
    pass
    #raise GeneiusError("Missing Package: bx / bx-python / bx-omelogic http://bitbucket.org/james_taylor/bx-python")

def reverse_complement( s ):
    complement_dna = {"A":"T", "T":"A", "C":"G", "G":"C", "a":"t", "t":"a", "c":"g", "g":"c", "N":"N", "n":"n" }
    reversed_s = []
    for i in s:
        reversed_s.append( complement_dna[i] )
    reversed_s.reverse()
    return "".join( reversed_s )

def fetch_sequence(genome, chrom, start, end, strand=1):
    try:
        twobitfile = bx.seq.twobit.TwoBitFile( file( genome ) )
    except:
        return "error: sequence file "+ genome +" unavailable; "
    s = twobitfile[chrom][start:end]
    if strand < 0 or strand == "-":
        return reverse_complement(s)
    else:
        return s

def fetch_sequences(genome, ranges):
    for rng in ranges:
        if not(rng.has_key("strand")):
            rng["strand"]="."
        rng["sequence"] = fetch_sequence(genome, rng["chr"], rng["start"], rng["end"], rng["strand"])
    return ranges

def fetch_spliced_sequence(genome, ranges):
    seq = ""
    for rng in ranges:
        if not(rng.has_key("strand")):
            rng["strand"]="."
        strand = rng["strand"]
        block = fetch_sequence(genome, rng["chr"], rng["start"], rng["end"], rng["strand"])
        if strand == "-" or strand < 0:
            seq = block+seq
        else:
            seq = seq+block
    return seq

def get_rna_ranges_for_mapping(mapping):
    ranges = []
    for exon in mapping["exons"]:
        ranges.append( {    "chr": mapping["chr"],
                            "start": exon["start"],
                            "end": exon["end"],
                            "strand": mapping["strand"] } )
    return ranges

def get_coding_ranges_for_mapping(mapping):
    ranges = []
    for cexon in mapping["cds"]:
        ranges.append( {    "chr": mapping["chr"],
                            "start": cexon[0],
                            "end": cexon[1],
                            "strand": mapping["strand"] } )
    return ranges

def fetch_mapping_dna(genomes_rule, mapping):
    genome = genomes_rule.replace("%",mapping["map_build"])
    seq = fetch_sequence(genome, mapping["chr"], mapping["start"], mapping["end"], mapping["strand"])
    return seq

def fetch_mapping_rna(genomes_rule, mapping):
    genome = genomes_rule.replace("%",mapping["map_build"])
    ranges = get_rna_ranges_for_mapping(mapping)
    rna = fetch_spliced_sequence( genome, ranges )
    if rna.startswith("error"):
        return rna
    rna = rna.upper().replace("T","U")
    return rna

def fetch_mapping_protein(genomes_rule, mapping):
    genome = genomes_rule.replace("%",mapping["map_build"])
    pranges = get_coding_ranges_for_mapping(mapping)
    cds = fetch_spliced_sequence( genome, pranges )
    if cds.startswith("error"):
        return cds
    return codons.translateAll( cds )
    
def fetch_protein_for_genes(genomes_rule, results):
    for rfsq in results:
        for mapping in rfsq["mappings"]:
            mapping["protein"] = fetch_mapping_protein( genomes_rule, mapping )
    return results

def fetch_rna_for_genes(genomes_rule, results):
    for rfsq in results:
        for mapping in rfsq["mappings"]:
            mapping["rna"] = fetch_mapping_rna( genomes_rule, mapping )
    return results

def fetch_dna_for_genes(genomes_rule, results):
    for rfsq in results:
        for mapping in rfsq["mappings"]:
            mapping["dna"] = fetch_mapping_dna( genomes_rule, mapping )
    return results

