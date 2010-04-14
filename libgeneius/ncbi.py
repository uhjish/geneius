from libgeneius.error import GeneiusError


def get_gid_for_refseq(refseq_id,ncbi_eutils):
    """looks up gid for refseq"""
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

    return IdList[0]

def add_sequence_to_refseqs(refseqs,ncbi_eutils,ncbi_sequence):
    """Looks up each refseq and adds the corresponding sequence
    using ncbi's webservices"""

    for refseq in refseqs:
        refseq_id = refseq['refseq_id']
                
        GID = get_gid_for_refseq(refseq_id,ncbi_eutils)

        try:
            sresponse = ncbi_sequence.service.run_eFetch("Nucleotide",GID)
        except:
            raise GeneiusError("Trouble Contacting NCBI eFetch")
        try:
            Sequence = sresponse.GBSet.GBSeq.GBSeq_sequence
        except:
            raise GeneiusError("Trouble Parsing NCBI eFetch for %s" % (refseq_id,sresponse))
        refseq['sequence'] = Sequence
        

def get_ncbi_entry_for_gid(gid,ncbi_sequence):
    try:
        sresponse = ncbi_sequence.service.run_eFetch("Nucleotide",gid)
    except:
        raise GeneiusError("Trouble Contacting NCBI eFetch")
    return sresponse
        

def seq_from_ncbi_data(ncbi_data):
    try:
        Sequence = ncbi_data.GBSet.GBSeq.GBSeq_sequence
    except:
        raise GeneiusError("Trouble getting sequence from  NCBI eFetch data")
    return Sequence


def definition_from_ncbi_data(ncbi_data):
    try:
        definition = ncbi_data.GBSet.GBSeq.GBSeq_definition
    except:
        raise GeneiusError("Trouble getting definition from NCBI eFetch data")
    return definition
