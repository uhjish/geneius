echo "Fetching GO Terms from NCBI"
curl ftp://ftp.ncbi.nlm.nih.gov/gene/DATA/gene2go.gz | zcat | grep -v -e'^#' | cut -f 2,3,6,8 | mawk 'BEGIN{FS="\t";OFS="\t";}{print $1, "GO-" $4, $2 "|" $3;}' >> entrez_anno

echo "Fetching KEGG Pathways from genome.jp"
curl ftp://ftp.genome.jp/pub/kegg/pathway/map_title.tab >kegg_path
curl ftp://ftp.genome.jp/pub/kegg/linkdb/genes/ncbi-geneid.ixl | grep "path:" | sed -e's/path://' -e's/ncbi-geneid://' >kegg_gene

echo "Collating KEGG Pathways to Entrez"
mawk '
BEGIN{
    FS="\t";
    OFS="\t";
    while(getline<"kegg_path" > 0){ 
        names[$1]=$2;
    } 
}
{ 
    kp=substr($2,4); 
    print $1, "KEGG",  $2 "|" names[kp]; 
}' kegg_gene >> entrez_anno

rm kegg_gene kegg_path

echo "Normalizing annotation tables"
rm tbl_anno_src tbl_anno_term tbl_anno
mawk '
BEGIN{ 
    FS="\t";
    OFS="\t";
    src_idx=0; 
    term_idx=0;
}
{ 
    if(!src[$2]){
        src_idx++;
        src[$2]=src_idx;
        print src_idx, $2 >>"tbl_anno_src";
    } 
    if(!term[$3]){
        term_idx++; 
        term[$3]=term_idx; 
        print term_idx, $3 >>"tbl_anno_term";
    } 
    print $1, src[$2], term[$3];
}' entrez_anno > tbl_anno
