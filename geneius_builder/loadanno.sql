
use geneius;

warnings;

drop table if exists tbl_anno;
drop table if exists tbl_anno_src;
drop table if exists tbl_anno_term;


create table tbl_anno_src(
    src_id INT NOT NULL,
    PRIMARY KEY (src_id),
    source VARCHAR(100),
    index source_idx(source)
);

create table tbl_anno_term(
    term_id INT NOT NULL,
    PRIMARY KEY (term_id),
    term VARCHAR(512),
    index term_idx(term)
);

create table tbl_anno(
	entrez_id INT NOT NULL,
    index entrez (entrez_id),
	src_id INT NOT NULL,
    term_id INT NOT NULL
);

LOAD DATA LOCAL INFILE "tbl_anno_src" into table tbl_anno_src;
LOAD DATA LOCAL INFILE "tbl_anno_term" into table tbl_anno_term;
LOAD DATA LOCAL INFILE "tbl_anno" into table tbl_anno;

