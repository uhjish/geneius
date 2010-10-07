
use geneius;

-- drop table if exists tbl_entrez_xref;
drop table if exists tbl_uniprot2refseq;


create table tbl_uniprot2refseq(
	uniprot varchar(100) not null,
    index uniprot_idx(uniprot),
	refseq_protein varchar(50),
	index refprot_idx (refseq_protein)	);

load data local infile "uniprot2refseq" into table tbl_uniprot2refseq;
