import sys

if (sys.argv[2] == "-"):
	listFile = sys.stdin
else:
	listFile = open(sys.argv[2],'r')

idList = set()
for line in listFile:
	idList.add(line.rstrip())

idCol = int(sys.argv[3])

idFile =  open(sys.argv[1],'r')
for line in idFile:
	cols = line.rstrip().split("\t")
	#print cols[idCol]
	try:
		if (cols[idCol] in idList):
			print line.rstrip()
	except:
		print "#Err: "+line

