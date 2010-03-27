#usage: split_ref.sh inputUCSCFile usefulColFile exonFile
awk '
BEGIN{
	FS="\t";
	OFS="\t";
}
{
	print NR,$2,$3,$4,$5,$6,$7,$8,$9 >>refMain; 
	gsub("/,$//",$10); 
	gsub("/,$//",$11); 
	split($10,starts,","); 
	split($11,ends,","); 
	i=1; 
	while (i <=$9){
		exnum=i;
		if ($4 == "-"){
			exnum=$9-(i-1);
		}	
		 print NR, exnum, starts[i], ends[i] >>refExon;
		 i++;
	}
}' refMain="$2" refExon="$3" $1
