if [ -f "/bio/ucsc/$1/$1.2bit" ]
then
    echo "$1 2bit exists, skipping\n"
else
    echo "Fetching $1.2bit\n"
    `curl ftp://hgdownload.cse.ucsc.edu/gbdb/$1/$1.2bit >/bio/ucsc/$1/$1.2bit`
fi
