#usage: split_ref.sh inputUCSCFile usefulColFile exonFile
awk '
BEGIN{
	FS="\t";
	OFS="\t";
	while (getline < "taxon"){
		gsub(" ","_",$2);
		species[$2]=$1;
	}
}
{
	print NR,species[$1],$3,$4,$5,$6,$7,$8,$9,$10 >>refMain; 
	gsub("/,$//",$11); 
	gsub("/,$//",$12); 
	split($11,starts,","); 
	split($12,ends,","); 
	i=1; 
	while (i <=$10){
		exnum=i;
		if ($5 == "-"){
			exnum=$10-(i-1);
		}	
		 print NR, exnum, starts[i], ends[i] >>refExon;
		 i++;
	}
}' refMain="$2" refExon="$3" $1
