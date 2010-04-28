def parsecodingexons(mapping):
    cds_st = mapping["cds_start"]
    cds_end = mapping["cds_end"]
    strand = mapping["strand"]
    lutr = []
    cds = []
    rutr = []
    mLength = 0
    for exon in mapping["exons"]: 
        num=exon["number"]
        st=exon["start"]
        en=exon["end"]
        mLength += en - st
        if(en <= cds_st): 
            ## completely in 5 utr
            lutr.append([st,en])
        elif st < cds_st and en < cds_end:
            ## partially in 5'utr
            lutr.append([st,cds_st])
            cds.append([cds_st,en])
        elif st < cds_st and en > cds_end:
            ## partially in 5'utr and partially in 3' utr
            lutr.append([st,cds_st])
            cds.append([cds_st,cds_end])
            rutr.append([cds_end,en])
        elif st >= cds_st and en <= cds_end:
            ## compeletely in cds
            cds.append([st,en])
        elif st < cds_end and en > cds_end:
            ## partially in utr3
            cds.append([st,cds_end])
            rutr.append([cds_end,en])
        elif st > cds_end:
            rutr.append([st,en])
    if strand == "+":
        utr5=lutr
        utr3=rutr
    else:
        utr5=rutr
        utr3=lutr
    pLength=mLength
    for u in lutr:
        pLength -= u[1]-u[0]
    for u in rutr:
        pLength -= u[1]-u[0]
    pLength = pLength/3 -1
    gstruct = { 
            "utr5" : utr5,
            "cds" : cds,
            "utr3" : utr3,
            "mrna_length": mLength,
            "protein_length": pLength
             }
    return gstruct

def addcoding(gene_result):
    for mping in gene_result["mappings"]:
        gst = parsecodingexons(mping)
        mping.update(gst)
    return gene_result
