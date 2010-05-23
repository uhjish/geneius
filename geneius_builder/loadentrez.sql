
use geneius;

-- drop table if exists tbl_entrez_xref;
drop table if exists tbl_gene_refseq;

/*
create table tbl_entrez_xref(
	entrez_id INT NOT NULL,
	PRIMARY KEY (entrez_id),
	species INT NOT NULL,
	index sp_idx (species),
	chr_locus VARCHAR(100),
	type VARCHAR(100),
	index typ_idx (type),
	official_symbol VARCHAR(100),
	official_gene_name VARCHAR(500),
	other_id VARCHAR(100),
	other_symbols VARCHAR(200),
	index names_idx (official_symbol, official_gene_name, other_id, other_symbols),
	other_gene_names TEXT
	fulltext index (other_gene_names)	);

*/

create table tbl_gene_refseq(
	entrez_id int not null,
    index entrez_idx(entrez_id),
	status varchar(50),
	refseq_rna varchar(50),
	refseq_protein varchar(50),
	index names_idx (refseq_rna, refseq_protein)	);

-- load data local infile "tbl_entrez_xref" into table tbl_entrez_xref;
load data local infile "tbl_gene_refseq" into table tbl_gene_refseq;
