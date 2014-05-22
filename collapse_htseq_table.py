import sys
import os
#Because the counts table from htseq-count comes transcript ids originally, there are multiple instances 
#of the same sequence name. The counts table is thus incompatible with DESeq2 and edgeR. Therefore,
#the counts from transcripts of the same ref sequence  much be combined or collapsed together as belonging 
#to one sequence. This script does that. Run "sudo python collapse_htseq_table.py <htseq.counts.txt> 
new_table = open("htseq.counts.collapsed.txt", "wr")
genes = {}
table = open(sys.argv[1])
lines = []


for line in table:
	line = line.split("\t")
	lines.append(line)

newList = []
for x in lines[0]:
	newList.append(x.split("\r"))

atoms = []
for x in newList:
	for y in x:
		atoms.append(y)

ROWS = []
for i in range(0, len(atoms)-37, 37):
	ROWS.append(atoms[i:i+37])

i=0
for x in ROWS:
	if i>0:
		if x[0] not in genes.keys():
			genes[x[0]] = x[1:]
		else:
			for i in range(1, len(genes[x[0]]), 1):
				genes[x[0]][i] = genes[x[0]][i] + x[i]
	else:
		i += 1

firstRow = ""
for x in ROWS[0]:
	if x != ROWS[0][len(ROWS[0])-1]:
		firstRow += x + "\t"
	else:
		firstRow += x

last_rows = []
last_rows.append(firstRow + "\n")

for key in sorted(genes):
	row = key + "\t"
	i = 1
	for x in genes[key]:
		if i == 36:
			row += x
		else:
			row += x + "\t"
            		i += 1
	last_rows.append(row + "\n") 

new_table.writelines(last_rows)




