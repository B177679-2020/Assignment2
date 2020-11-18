#!/usr/bin/python3

##FIRST STEP: GET ALL THE SEQUENCES

import os, sys
import subprocess
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

'''
#Example query to get the birds G6P protein sequences in fasta format
#ESearchBirds = 'esearch -db protein -query "txid8782[Organism] AND glucose-6-phosphatase[Protein]" | efetch -db protein -format fasta > birds_G6P.fasta'

#ESearchBText = subprocess.check_output(ESearchBirds, shell=True)

TaxonID = input("Please type in your Taxon ID (4-number digit):")
Protein = input("Please type in the protein of interest (ex.: phosphatase, kinase...):")

GeneralES = 'esearch -db protein -query "{}[Organism] AND {}*[Protein] NOT partial NOT predicted" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)

GeneralESout = subprocess.check_output(GeneralES, shell=True)

#Check the user file can be open in Python and it is correct
UserFile = open('{}.fasta'.format(TaxonID))
UserFile_Contents = UserFile.read() 
print(UserFile_Contents)

num = len([1 for line in UserFile_Contents if line.startswith(">")])
print("The number of protein sequences for your query is {}".format(num))

##MULTIPLE SEQUENCE ALIGNMENT WITH CLUSTALO
Clustalo = 'clustalo -i {}.fasta -o {}MA.fasta -v'.format(TaxonID, TaxonID)
ClustaloOut = subprocess.check_output(Clustalo, shell=True)

#Check the ClustalO output
UserClustal = open('{}MA.fasta'.format(TaxonID))
UserClustal_Contents = UserClustal.read()
print(UserClustal_Contents)
numClustal = len([1 for line in UserClustal_Contents if line.startswith(">")])
print("The number of aligned sequences with clustal is {}".format(numClustal))

#Get consensus sequence
if numClustal > 250:
	Cons = 'cons -sequence {}MA.fasta -outseq {}Cons.fasta'.format(TaxonID, TaxonID)
	ConsOut = subprocess.check_output(Cons, shell=True)
	BlastDB = 'makeblastdb -in {}.fasta -dbtype prot -out {}'.format(TaxonID, TaxonID)

#Try with blast
#blastp -query txid8782Cons.fasta -num_alignments 250 -db txid8782.fasta -outfmt 7 -out blastout.txt

##Sort the blast output by the value of identity (3 column)
#sort -k3,3nr blastout.txt


"""
##Sort using pandas
df = pd.read_csv('/localdisk/home/s2100021/Assignment2/try2.txt',sep='\t', skiprows=(0,1,2,3,4))
df.columns = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
dfsorted = df.sort_values('3', ascending = False)
AccNumbs = str(dfsorted.iloc[1:250, 1])

UserFile = open('txid8782.fasta')
UserFile_Contents = UserFile.read()
print(UserFile_Contents)


seqs250=[]
for AccNumb in AccNumbs:
	for line in UserFile_Contents:
		if AccNumb in line:
			seqs250.append(line)



		

