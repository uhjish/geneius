from libgeneius.error import GeneiusError
try:
    import bx.seq.twobit
except:
    raise GeneiusError("Missing Package: bx / bx-python / bx-omelogic")

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

