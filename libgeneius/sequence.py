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
    twobitfile = bx.seq.twobit.TwoBitFile( file( genome ) )
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

def get_ranges_for_mapping(mapping):
    ranges = []
    for exon in mapping["exons"]:
        ranges.append( {    "chr": mapping["chr"],
                            "start": exon["start"],
                            "end": exon["end"],
                            "strand": mapping["strand"] } )
    return ranges

def fetch_mapping_sequence(genomes_rule, mapping):
    genome = genomes_rule.replace("%",mapping["map_build"])
    ranges = get_ranges_for_mapping(mapping)
    return fetch_spliced_sequence( genome, ranges )


def fetch_gene_sequences(genomes_rule, results):
    for rfsq in results:
        for mapping in rfsq["mappings"]:
            mapping["sequence"] = fetch_mapping_sequence( genome, ranges )
    return results
