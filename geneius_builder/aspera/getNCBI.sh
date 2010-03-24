
originalName=$1
relName=${originalName##ftp*gov/}

retVal=-1
tries=0
while [ $retVal -ne 0 ]
do
	echo $originalName
	/bio/aspera/bin/ascp  -QT -m 20M -l 1000M -i /bio/aspera/etc/asperaweb_id_dsa.putty anonftp@ftp-private.ncbi.nlm.nih.gov:$relName .
	retVal=$?
	if [ $retVal -ne 0 ]
	then
		echo "failed... trying again in 10 seconds."
		sleep 10
	fi
	
done
