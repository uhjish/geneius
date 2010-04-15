#usage: split_ref.sh inputUCSCFile usefulColFile exonFile
awk '
BEGIN{
	FS="\t";
	OFS="\t";
	while (getline < "refMain.txt"){
		refids[$3]=$1;
        print $3
	}
}
{
    if (refids[$2] != 0){
        print refids[$2], $1, $3;
    }
}'  $1
