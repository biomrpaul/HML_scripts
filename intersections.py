import sys
import os

sigGenes = {}
path="shrimp_dge/"
for dir in os.listdir(path):	
	for sub in os.listdir(path + dir):
		for list in os.listdir(path + dir + "/" +  sub):
			newpath = path + dir + "/" +  sub + "/" + list
			if list[0:5] == "genes" or list[len(list)-4:len(list)] == "cuff":
				if list.find('edgeR') > -1:
					name= list[6:list.find('t')-1]
				elif list[len(list)-4:len(list)] != "cuff":
					end = list.find('results')
					name= list[6:end-1]
				else:
					name = list
				genes = []
				for line in open(newpath):
					if list[len(list)-4:len(list)] == "cuff":
						line = line[0:line.find(":")] + "\n"
					genes.append(line)
				sigGenes[name] = genes
							
output = open("dge-report.txt", "w")				
for key in sorted(sigGenes):
	list = sigGenes[key]
	if len(list) == 1 and list[0]  == "\n":
		list.pop()
	output.write(key + " has " + str(len(list)) + " significant genes. \n")
output.write("\n")
pairs = []
for key in sigGenes:
	for lock in sigGenes:
		if key != lock and (key + lock) not in pairs:
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
