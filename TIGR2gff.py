import re
import sys

input=sys.argv[1]
output=sys.argv[2]

inFile = open(input,"r")
outFile = open(output,"w")

for lines in inFile:
        lineList = lines.split()
#	print lineList
        if lineList[0] != "ORF_ID":
        #if lineList[11] == "protein homolog model":
                posCapt = lineList[0].split('|')
#               print(posCapt[6])
                _seqid = posCapt[6]
                _source = "TIGR"
                _type = "func_gene"
                _start = posCapt[4]
                _stop = posCapt[5]
                _score = lineList[5]
                _strand = posCapt[3]
                _phase = "0"
                _attributes = "gene_id=" + posCapt[0] + "; " +  "TIGRfam=" + lineList[2]

                outFile.write( _seqid + "\t" + _source + "\t" + _type + "\t" + _start + "\t" + _stop + "\t" + _score + "\t" + _strand + "\t" + _phase + "\t" + _attributes + "\n")

inFile.close()
outFile.close()
