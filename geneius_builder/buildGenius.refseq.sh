
#rm -f refGene.all.ucsc
#rm -f refMain.txt
#rm -f refExons.txt

echo "Downloading Entrez Gene table"
#./aspera/getNCBI.sh ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/All_Data.gene_info.gz   

echo "Parsing Entrez Gene table"
#zcat All_Data.gene_info.gz | gawk 'BEGIN{FS="\t";OFS="\t"}{print $2,$1,$8,$10,$3,$9,$6,$5,$14}' - | sed -e's/\t-/\tNULL/g' >tbl_entrez_xref


echo "Downloading Entrez to Refseq mappings"
#./aspera/getNCBI.sh ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/gene2refseq.gz

echo "Parsing Entrez to Refseq mappings"
#zcat gene2refseq.gz | gawk 'BEGIN{FS="\t";OFS="\t";}(NR>1){print $2,$3,$4,$6}' |  sed -e's/-/NULL/g' | sed -e's/\.[0-9]*//g' >tbl_gene_refseq


#get species master list
echo "Getting UCSC species list"
#curl ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/ | awk '(NR>2){gsub(/\.\.\//,"",$11);print $9 "\t" $11}' >spList

echo "Getting Refseq annos for all species"
#cut -f 1 spList | xargs -I % echo 'echo % && curl ftp://hgdownload.cse.ucsc.edu/goldenPath/currentGenomes/%/database/refFlat.txt.gz | zcat | mawk "{print \"%\t\" \$0;}" >>refGene.all.ucsc'| sh 
#$found = `wc -l refGene.all.ucsc`
echo "Found $found Refseq records"

echo "Parsing Refseq annos"
#sh split_ref.sh refGene.all.ucsc refMain.txt refExons.txt

echo "Updating taxon table"
#curl ftp://hgdownload.cse.ucsc.edu/goldenPath/uniProt/database/taxon.txt.gz | zcat | cut -f 1,2 >taxon

echo "Building species table"
#cat spList | tr "_" " " | python matchids.py taxon - 1 >tbl_species

echo "Building getting sequences & headers"
#mkdir tmp_fa
#cd tmp_fa
#../aspera/getNCBIbatch.sh ftp://ftp.ncbi.nlm.nih.gov/refseq/release/complete/ rna.fna.gz
#cd ..
zcat tmp_fa/*.rna.fna.gz | grep -e'^>' | sed -e 's/^>gi|//' -e 's/ref|//' | tr "|" "\t" | mawk 'BEGIN{FS="\t"; OFS="\t"}{gsub(/\.[0-9]+/,"",$2); print $0;}' >refDesc

