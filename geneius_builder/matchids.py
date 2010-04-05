import sys

if (sys.argv[2] == "-"):
    listFile = sys.stdin
else:
    listFile = open(sys.argv[2],'r')

otherFlds={}

idList = set()
for line in listFile:
    flds = line.rstrip().split("\t")
    idList.add(flds[0])
    otherFlds[flds[0]]=flds[1:]

idCol = int(sys.argv[3])

idFile =  open(sys.argv[1],'r')
for line in idFile:
    cols = line.rstrip().split("\t")
    #print cols[idCol]
    try:
        if (cols[idCol] in idList):
            outln = line.rstrip()
            if (len(otherFlds[cols[idCol]])>0):
                outln = outln + "\t" + "\t".join(otherFlds[cols[idCol]])
            print outln
    except:
        print "#Err: "+line

