#!/usr/bin/env bash
genome=$1
echo $genome
#if [ -f "/bio/ucsc/${genome}/${genome}.2bit" ]
#then
#    echo ${genome} 2bit exists, skipping
#else
    echo "Fetching ${genome}.2bit to /bio/ucsc/${genome}/${genome}.2bit"
    `curl ftp://hgdownload.cse.ucsc.edu/gbdb/${genome}/${genome}.2bit >/bio/ucsc/${genome}/${genome}.2bit`
fi
