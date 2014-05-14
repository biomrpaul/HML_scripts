import sys

output = open("locus.stripped.txt", 'w')
with open(sys.argv[1]) as file:
	i = 0
	for line in file:
		line = line.split()
		if i>0:
			output.write(line[1][1:line[1].find(":")]+"\n")
		i += 1
print "loci stripped"
