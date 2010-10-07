from libgeneius.simple_geneius import SimpleGeneius

sg = SimpleGeneius()

res = sg.get_symbols_for_entrez()

print res[34]
