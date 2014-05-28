hypo = open("N24vHH24_consensus_genes.txt", "r")
hyper = open("N24vH24_consensus_genes.txt", "r")
newHypo = open("N24vHH24_relation_with_H24", "w")

hypoList = hypo.readline().split("\t")
for x in hypoList:
	hypoList[hypoList.index(x)] = x.split("\r")

hyperList = hyper.readline().split("\t")
for x in hyperList:
	hyperList[hyperList.index(x)] = x.split("\r")

newHypoList = []
for x in hypoList:
	for y in x:
		newHypoList.append(y)

newHyperList = []
for x in hyperList:
	for y in x:
		newHyperList.append(y)

hypoTable = []
hyperTable = []
for i in range(0, len(newHypoList),7):
	hypoTable.append(newHypoList[i:i+7])

for i in range(0, len(newHyperList),7):
	hyperTable.append(newHyperList[i:i+7])

firstRow = "contigs\tGO\tH24_cuff_FC\tH24_deseq_FC\tRelation\tHH24_cuff_FC\tHH24_deseq_FC\n"

newHypo.write(firstRow)

hyperContigs = {}
for i in range(1,len(hyperTable)):
	hyperContigs[hyperTable[i][0]] = hyperTable[i][1:]


for i in range(1, len(hypoTable)):
	row = hypoTable[i][0] + "\t" + hypoTable[i][5] + "\t" + hypoTable[i][2] + "\t" + hypoTable[i][4] + "\t"
	result = ""
	cuffFC = ""
	deseq = ""
	contig = hypoTable[i][0]
	if contig not in hyperContigs.keys():
		result = "Not shared"
		cuffFC = "NA"
		deseq = "NA"
	elif float(hypoTable[i][2]) > 0 and float(hyperContigs[contig][1]) > 0:
		result = "+:+"
		cuffFC = hyperContigs[contig][1]
		deseq = hyperContigs[contig][3]
	elif float(hypoTable[i][2]) < 0 and float(hyperContigs[contig][1]) < 0:
		result = "-:-"
		cuffFC = hyperContigs[contig][1]
		deseq = hyperContigs[contig][3]
	elif float(hypoTable[i][2]) > 0 and float(hyperContigs[contig][1]) < 0:
		result = "+:-"
		cuffFC = hyperContigs[contig][1]
		deseq = hyperContigs[contig][3]
	else:
		result = "-:+"
		cuffFC = hyperContigs[contig][1]
		deseq = hyperContigs[contig][3]

	row += result + "\t" + cuffFC + "\t" + deseq + "\n"
	newHypo.write(row)










