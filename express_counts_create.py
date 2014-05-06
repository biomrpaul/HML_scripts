import sys
import csv
import numpy
import array
import re
import pandas as pd
#Produced by Matt Paul, 5.5.2014
#This tool converts the expected counts from results.xprs file outputs from eXpress and creates a counts table to be used for DESeq2 and edgeR. Create a folder that contains all of the express result folders (sampleName.bam_xprs) and run this command "python counts_create.py /filepath_to_outputs/*.bam_xprs"


#Finds the first word seperated by a '.', assumes that is the sample name
period_pos = sys.argv[1].find('.')
name = sys.argv[1][:period_pos]

#Creates a dictionary to store gene_id/count pairs
geneToCount = {}

#Goes into the first result directory and opens up the results file
firstFile = open("" + sys.argv[1] + "/results.xprs")

#Populates the dictionary with gene/count pairs, by splitting up the row by white spaces and calling the 2 and 5 columns (gene and count, resp.)
i = 0
for line in firstFile:
    if i>0:
        columns = line.split()
        geneToCount[columns[1]] = columns[4]
    else:
        i +=1

#Creates an initializes the first two arrays which serve as the first two columns
gene_ids = []
gene_ids.append("gene_id")
first_counts = []
first_counts.append(name)

#Populates the first two arrays with the sorted gene/count pairs
for key in sorted(geneToCount):
    gene_ids.append(key)
    first_counts.append(geneToCount[key])

#Created a two-dimensional array that contains the first collumns (arrays)
counts_table= [gene_ids, first_counts]

#Now that there is a 2D array with the gene_ids and the first sample's expected counts, this block creates and populates and array of exp. counts from each sample then adds it to the 2D array, effectively adding a new column
for i in range(2, len(sys.argv)):
    period_pos = sys.argv[i].index('.')
    name = sys.argv[i][:period_pos]
    file = open("" + sys.argv[i] + "/results.xprs")
    gene_dict = {}
    i = 0
    for line in file:
        if i>0:
            columns = line.split()
            gene_dict[columns[1]] = columns[4]
        else:
            i +=1

    counts=[]
    counts.append(name)

    for key in sorted(gene_dict):
        counts.append(gene_dict[key])

    counts_table.append(counts)

#Make the 2D array into a numpy array (maybe not necessary?)
ar = numpy.array(counts_table)
#Create a new csv file and initiate a csv writer for that file
output = open('express.counts.txt', 'w')


#The for loop iterates though every element in the 2D arrays via rows (i). The collumns are accessed in the array (x). Therefore, x[i] is the expected count value for sample x.
for i in range(0, len(gene_ids)-1):
    row = ""
    for x in ar:
        row += x[i] + "\t"
    output.write(row + "\n")

output.close()

print "well it ran"
