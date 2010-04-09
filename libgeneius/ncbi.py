from libgeneius.error import GeneiusError

def add_sequence_to_refseqs(refseqs,ncbi_eutils,ncbi_sequence):
    """Looks up each refseq and adds the corresponding sequence
    using ncbi's webservices"""

    for refseq in refseqs:
        refseq_id = refseq['refseq_id']
        try:
            response = ncbi_eutils.service.run_eSearch("Nucleotide",refseq_id)
        except:
            raise GeneiusError("Trouble Contacting NCBI eSearch")
        try:
            IdList = response.IdList.Id
        except:
            raise GeneiusError("Trouble parsing NCBI eSearch for %s" % refseq_id)
        if len(IdList) < 1:
            raise GeneiusError("Not results for accession %s" % refseq_id)
        if len(IdList) > 1:
            raise GeneiusError("Found multiple ID's for %s:%s" % (refseq_id,IdList))
        
        GID = IdList[0]

        try:
            sresponse = ncbi_sequence.service.run_eFetch("Nucleotide",GID)
        except:
            raise GeneiusError("Trouble Contacting NCBI eFetch")
        try:
            Sequence = sresponse.GBSet.GBSeq.GBSeq_sequence
        except:
            raise GeneiusError("Trouble Parsing NCBI eFetch for %s" % (refseq_id,sresponse))
        refseq['sequence'] = Sequence
        
