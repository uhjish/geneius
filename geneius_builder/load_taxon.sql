
use geneius;

drop table if exists tbl_species;

create table tbl_species(
	tax_id INT NOT NULL,
	PRIMARY KEY (tax_id),
	name VARCHAR(500) ,
    build VARCHAR(100));

load data local infile "tbl_species" into table tbl_species;
