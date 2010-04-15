use geneius;
warnings;
drop table if exists tbl_refDesc;

create table tbl_refDesc(
	gid int not null,
	refseq_id varchar(100),
	primary key(refseq_id),
	index ref_idx (refseq_id),
	description text);

load data local infile "refDesc" into table tbl_refDesc;

