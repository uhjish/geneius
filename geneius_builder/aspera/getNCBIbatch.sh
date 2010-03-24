originalName=$1
curl $1 | awk '{print $9}' | grep "$2" |  xargs -I % /bio/aspera/getNCBI.sh $1/%
