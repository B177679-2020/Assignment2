#!/usr/bin/python3

##FIRST STEP: GET ALL THE SEQUENCES

#Import all modules that the program needs to run
import os, sys
import subprocess
import numpy as np
import pandas as pd


#Make sure we are running in the user workspace of the server
chdir = ("cd /localdisk/home/$USER")
subprocess.check_output(chdir, shell=True)
#Create a new directory to store the results
Directory = input("Type in a name for the directory that will be created in your workspace with the results:")
os.mkdir(Directory)
os.chdir(Directory)
#Create source directory to keep results tidy
os.mkdir("src")
os.chdir("src")

#Example query to get the birds G6P protein sequences in fasta format
#ESearchBirds = 'esearch -db protein -query "txid8782[Organism] AND glucose-6-phosphatase[Protein]" | efetch -db protein -format fasta > birds_G6P.fasta'

#ESearchBText = subprocess.check_output(ESearchBirds, shell=True)

##GET USER INPUT w. some error traps for typing errors
TaxonID = input("Please type in your Taxon ID (txidxxxx):")
if not TaxonID.startswith("txid"):
	print("Please, enter the Taxon ID in the correct format (txid followed by 4 or 5 numbers)")
	TaxonID = input("Please type in your Taxon ID again (remember, txidxxxx):")
Protein = input("Please type in the protein of interest (ex.: phosphatase, kinase...):")

##Search in the NCBI database
GeneralES = 'esearch -db protein -query "{}[Organism] AND {}*[Protein] NOT partial NOT predicted" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)

GeneralESout = subprocess.check_output(GeneralES, shell=True)

#Check the user file can be open in Python and it is correct
UserFile = open('{}.fasta'.format(TaxonID))
UserFile_Contents = UserFile.read() 
print(UserFile_Contents)
num = len([1 for line in UserFile_Contents if line.startswith(">")])
print("The number of protein sequences for your query is {}".format(num))
if num == 0:
	print("Your query did not yield any results! Are you sure you typed in the gene name and the taxon ID correctly...? Please, start again")
	exit()

UserFile.close()

##MULTIPLE SEQUENCE ALIGNMENT WITH CLUSTALO
print("ClustalO multiple alignment is being performed. Please be a little bit patient...")
Clustalo = 'clustalo -i {}.fasta -o {}MA.fasta -v --output-order tree-order'.format(TaxonID, TaxonID)
ClustaloOut = subprocess.check_output(Clustalo, shell=True)

#Check the ClustalO output
UserClustal = open('{}MA.fasta'.format(TaxonID))
UserClustal_Contents = UserClustal.read()
print(UserClustal_Contents)
numClustal = len([1 for line in UserClustal_Contents if line.startswith(">")])
print("The number of aligned sequences with clustal is {}".format(numClustal))

#Get consensus sequence
if numClustal > 250:
	##Get a consensus sequence to BLAST against
	print("There were more than 250 sequences. BLAST analysis will be performed to keep only the 250 most similar ones")
	Cons = 'cons -sequence {}MA.fasta -outseq {}Cons.fasta'.format(TaxonID, TaxonID)
	ConsOut = subprocess.check_output(Cons, shell=True)
	print("A consensus sequence was generated")
	#BLAST the generated consensus seq against all the prot sequences found in the query
	#First generate a BLAST database to search against
	BlastDB = 'makeblastdb -in {}.fasta -dbtype prot'.format(TaxonID)
	#Perform the BLAST search
	BlastP = 'blastp -query {}Cons.fasta -num_alignments 250 -db {}.fasta -outfmt 7 -out blastout.txt'.format(TaxonID, TaxonID)
	BlastDBOut = subprocess.check_output(BlastDB, shell=True)
	BlastPOut = subprocess.check_output(BlastP, shell=True)
	print("BLAST analysis performed")


	#Sort BLAST results by similarity using Pandas
	df = pd.read_csv('blastout.txt',sep='\t', skiprows=(0,1,2,3,4))
	print(df)
	df.columns = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
	pd.concat([pd.DataFrame(df.columns), df], ignore_index=True)
	dfsorted = df.sort_values('3', ascending = False)
	print("Sorted by similarity")
	print(dfsorted)

	#This gets the accession number of the 250 more highly related proteins
	#Get accession numbers as a list
	listID = df['2'].to_list()
	ListID250 = listID[0:250]
	print("Acc numbers of the most similar sequences are shown")
	print(len(ListID250))

	FinalFASTA = open("SimilarSeqs.txt", mode="w")
	AI_DICT = {}
	for line in ListID250:
    		AI_DICT[line[:-1]] = 1

	UserFile = open('{}.fasta'.format(TaxonID), 'r')
	

	skip = 0
	for line in UserFile:
        	if line[0] == '>':
                	_splitline = line.split(' ')
                	AccNumbArrow = _splitline[0]
                	accessorID = AccNumbArrow[1:-1]
                	print(accessorID)
                	if accessorID in AI_DICT:
                        	FinalFASTA.write(line)
                        	skip = 0
                	else:
                        	skip = 1
        	else:
                	if not skip:
                        	FinalFASTA.write(line)
	UserFile.close()
	FinalFASTA.close()
	Plot = 'plotcon SimilarSeqs.txt -graph svg -goutfile plotcon'
	subprocess.check_output(Plot, shell=True)
else:
##Generate the plotcon with original dataset if the number of sequences was <250
	Plot2 = 'plotcon {}.fasta -graph svg -goutfile plotcon'.format(TaxonID)
	subprocess.check_output(Plot2, shell=True)

##LEVEL OF SIMILARITY BETWEEN SEQUENCES


UserFile.close()

