import sys
import csv
import numpy
import array
import pandas as pd
#Produced by Matt Paul, 5.5.2014
#This tool converts the expected counts from gene.results file outputs from RSEM and creates a counts table to
# be used for DESeq2 and edgeR. Create a folder that contains all of the RSEM outputs (sampleName.gene.results)
# , copy this script into it, and run this command "python counts_create.py *.gene.results"

#Finds the first word seperated by a '.', assumes that is the sample name
period_pos = sys.argv[1].find('.')
name = sys.argv[1][:period_pos]

#opens first file
firstFile = open(sys.argv[1], 'rb')
#reads the file as a csv
firstData = csv.reader(firstFile, delimiter='\t')

#Create 2 arrays that serve as the first collumn in the table
gene_ids = []
first_counts = []

#This block populates the first two collumns (arrays), had a i=0 case which is equivalent to the header
i = 0
for line in firstData:
    if i != 0:
        gene_ids.append(line[0])
        first_counts.append(int(round(float(line[4]), 0)))
    else:
        gene_ids.append("gene_id")
        first_counts.append(name)
    i += 1

#Created a two-dimensional array that contains the first collumns (arrays)
counts_table= [gene_ids, first_counts]

#Now that there is a 2D array with the gene_ids and the first sample's expected counts, this block creates and populates and array of exp. counts from each sample then adds it to the 2D array, effectively adding a new column
for i in range(2, len(sys.argv)):
    
    period_pos = sys.argv[i].find('.')
    name = sys.argv[i][:period_pos]
    file = open(sys.argv[i], 'rb')
    data = csv.reader(file, delimiter='\t')
    counts=[]

    m = 0
    for line in data:
        if m != 0:
            counts.append(int(round(float(line[4]), 0)))
        else:
            counts.append(name)
            m += 1

    counts_table.append(counts)

#Make the 2D array into a numpy array (maybe not necessary?)
ar = numpy.array(counts_table)
#Create a new csv file and initiate a csv writer for that file
output = open('rsem.counts.txt', 'w')


#To a print the count table as a tab delimited file, ehe for loop iterates though every element in the 2D arrays via rows (i). The collumns are accessed in the array (x). Therefore, x[i] is the expected count value for sample x.
for i in range(0, len(gene_ids)-1):
    row = ""
    for x in ar:
        row += x[i] + "\t"
    output.write(row + "\n")
output.close()

print "well it ran"
