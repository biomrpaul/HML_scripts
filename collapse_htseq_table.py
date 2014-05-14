import sys

new_table = open("htseq.counts.collapsed.txt", "w")
genes = {}
with open(sys.argv[1]) as table:
	i = 0
	for line in table:
		if i>0:
			if line[0] not in genes.keys():
				genes[line[0]] = line[1:]
			else:
				for i in range(len(genes[line[0]])):
					genes[line[0]][i] = genes[line[0]][i] + line[i]
		else:
			new_table.write(line + "\n")
			i += 1

for key in genes:
	row = key
	for x in genes[key]:
		row += x + "\t"
	new_table.write(row + "\n") 

