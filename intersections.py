import sys
import os
import math
import shutil
import csv

folder = sys.argv[1]

if os.path.isdir(folder + "_con_genes") == True:
	shutil.rmtree(folder + "_con_genes")

cuff = open("gene_exp.diff", "r")

PW_anals = {}
PW_anals["N4vH4"] = {}
PW_anals["N4vHH4"] = {}
PW_anals["N24vH24"] = {}
PW_anals["N24vHH24"] = {}
PW_anals["H4vHH4"] = {}
PW_anals["H24vHH24"] = {}
PW_anals["N4vN24"] = {}
PW_anals["H4vH24"] = {}
PW_anals["HH4vHH24"] = {}

cuffList = {}

i = 0
for x in cuff:
	row = x.split("\t")
	if row[4] == "N4" and row[5] == "H4":
		PW_anals["N4vH4"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "N4" and row[5] == "HH4":
		PW_anals["N4vHH4"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "N24" and row[5] == "H24":
		PW_anals["N24vH24"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "N24" and row[5] == "HH24":
		PW_anals["N24vHH24"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "H4" and row[5] == "HH4":
		PW_anals["H4vHH4"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "H24" and row[5] == "HH24":
		PW_anals["H24vHH24"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "N4" and row[5] == "N24":
		PW_anals["N4vN24"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "H4" and row[5] == "H24":
		PW_anals["H4vH24"][row[3][0:row[3].find(":")]] = row[9]
	elif row[4] == "HH4" and row[5] == "HH24":
		PW_anals["HH4vHH24"][row[3][0:row[3].find(":")]] = row[9]
	else:
		l = 0


sigGenes = {}
path=folder + "/"
for dir in os.listdir(path):	
	for sub in os.listdir(path + dir):
		for list in os.listdir(path + dir + "/" +  sub):
			newpath = path + dir + "/" +  sub + "/" + list
			if list[0:5] == "genes" or list[len(list)-8:len(list)] == "cuff.txt":
				if list.find('edgeR') > -1:
					name= list[6:list.find('txt')-1]
				elif list[len(list)-8:len(list)] == "cuff.txt" :
					name = list[:list.find('txt')-1]
				else:
					end = list.find('results')
					name= list[6:end-1]
	
				genes = []
				for line in open(newpath):
					if list[len(list)-8:len(list)] == "cuff.txt":
						line = line[0:line.find(":")] + "\n"
					genes.append(line)
				sigGenes[name] = genes

os.makedirs(folder + "_con_genes")

output = open(folder + "_con_genes/dge-report.txt", "w")
for key in sorted(sigGenes):
	list = sigGenes[key]
	if len(list) == 1 and list[0]  == "\n":
		list.pop()
	output.write(key + " has " + str(len(list)) + " significant genes. \n")
output.write("\n")


pairs = []
anals = {}
anals["N4vH4"] = set()
anals["N4vHH4"] = set()
anals["N24vH24"] = set()
anals["N24vHH24"] = set()
anals["H4vHH4"] = set()
anals["H24vHH24"] = set()
anals["N4vN24"] = set()
anals["H4vH24"] = set()
anals["HH4vHH24"] = set()

for key in sigGenes:
	for lock in sigGenes:
		if key != lock and (key + lock) not in pairs:
			keyName = key[0:key.find(".")]
			lockName = lock[0:lock.find(".")]
			pairs.append(key+lock)
			pairs.append(lock+key)
			list1 = sigGenes[key]
			list2 = sigGenes[lock]
			inter = set(list1).intersection(set(list2))
			inter2 = []
			for x in inter:
				inter2.append(x)
			if len(inter) > 0:			
				distrib = str(len(list1)-len(inter)) + ":" + str(len(inter)) + ":" + str(len(list2) - len(inter))
				output.write("The distribution of " + key + " and " + lock + " is " + distrib + "\n") 
				genes = ""
				for atom in inter:
					genes = genes + atom			
				output.write("The shared contigs are:\n" +  genes + "\n")
			if keyName == lockName:
				if anals[keyName] == None:
					anals[keyName] = set(inter)
				else:
					anals[keyName] = anals[keyName] | set(inter)

annot = open("macro.csv", "rU")
annotReader = csv.reader(annot)

annots = {}
for row in annotReader:
	annots[row[1]] = [row[2], row[11]]



for key in anals:
	table = open(folder + "_con_genes/" + key + "_consensus_genes.txt", "w")
	rsemTable = open(folder + "/deseq/rsem/" + key + ".deseq.rsem.results.txt", "r")
	expressTable=open(folder + "/deseq/express/" + key + ".deseq.express.results.txt", "r")
	table.write("contigs\tcuffdiff_FC\tdeseq.express_FC\tdeseq.RSEM_FC\tGO\tDescription\n")
	rsemlogLookUp = {}
	explogLookUp = {}
	for y in rsemTable:
		row = y.split("\t")
		rsemlogLookUp[row[0].strip("\"")] = row[2].strip("\"")
	for y in expressTable:
		row = y.split("\t")
		explogLookUp[row[0].strip("\"")] = row[2].strip("\"")
	
	
	for x in anals[key]:
		cufflog2 = PW_anals[key][x.strip("\n")]
		rsemlog2 = rsemlogLookUp[x.strip("\n")]
		explog2 = explogLookUp[x.strip("\n")]
		
		if cufflog2 != "NA":
			if float(cufflog2) < 0:
				cuffFoldChange = -1 * (math.pow(2, abs(float(cufflog2))))
			else:
				cuffFoldChange = math.pow(2, float(cufflog2))
		else:
			cuffFoldChange = "NA"

		if rsemlog2 != "NA":
			if float(rsemlog2) < 0:
				rsemFoldChange = -1 * (math.pow(2, abs(float(rsemlog2))))
			else:
				rsemFoldChange = math.pow(2, float(rsemlog2))
		else:
			rsemFoldChange = "NA"

		if explog2 != "NA":
			if float(explog2) < 0:
				expFoldChange = -1 * (math.pow(2, abs(float(explog2))))
			else:
				expFoldChange = math.pow(2, float(explog2))
		else:
			expFoldChange = "NA"

		lock = x.strip("\n")
		if lock in annots.keys():
			row = "" + x.strip("\n") + "\t" + str(cuffFoldChange) +  "\t" + str(expFoldChange) + "\t" + str(rsemFoldChange) + "\t" + annots[x.strip("\n").strip("\"")][0] + "\t" + annots[x.strip("\n").strip("\"")][1] + "\n"
		else:
			row = "" + x.strip("\n") + "\t"  + str(cuffFoldChange) + "\t" + str(expFoldChange) + "\t" + str(rsemFoldChange) + "\t" + "ERROR: look up annot." + "\t" + "ERROR: look up annot." + "\n"

		
		table.write(row)






