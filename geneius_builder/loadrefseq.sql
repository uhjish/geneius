use geneius;
warnings;
drop table if exists tbl_refMain;
drop table if exists tbl_refExon;

create table tbl_refMain(
	id int not null,
	primary key(id),
	refseq_id varchar(100),
	index ref_idx (refseq_id),
	chr varchar(100),
	strand char(1),
	start int,
	end int,
	cds_start int,
	cds_end int,
	num_exons int);

create table tbl_refExon(
	ref_id int not null,
	foreign key (ref_id) REFERENCES tbl_refMain(id),
	exon_num int,
	exon_start int,
	exon_end int);

load data local infile "refMain2" into table tbl_refMain;
load data local infile "refExon2" into table tbl_refExon;


