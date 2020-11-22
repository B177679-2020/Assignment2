#!/usr/bin/python3

##FIRST STEP: GET ALL THE SEQUENCES

#Import all modules that the program needs to run
import os, sys
import subprocess
import numpy as np
import pandas as pd
import shutil
import re

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

##Some options for the user
yes = {'yes','y', 'ye', 'YES', 'Yes'}
no = {'no','n', 'No', 'NO'}

PartialsPred = input("Do you want to include partial and predicted sequences of the protein in your search, if there are any? yes or no")

##Search in the NCBI database

if PartialsPred in yes:
	GeneralES = 'esearch -db protein -query "{}[Organism] AND {}*[Protein]" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)
elif PartialsPred in no:
	GeneralES = 'esearch -db protein -query "{}[Organism] AND {}*[Protein] NOT partial NOT predicted" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)
else:
	print("No valid answer!")
	exit()

GeneralESout = subprocess.check_output(GeneralES, shell=True)


#Check the user file can be open in Python and it has a reasonable number of sequences
UserFile = open('{}.fasta'.format(TaxonID))
UserFile_Contents = UserFile.read() 
print(UserFile_Contents)
num = len([1 for line in UserFile_Contents if line.startswith(">")])
print("The number of protein sequences for your query is {}".format(num))
if num == 0:
	print("Your query did not yield any results! Are you sure you typed in the gene name and the taxon ID correctly...? Please, start again")
	exit()
elif num >= 1000:
	print("Your query yielded more than 1000 results! That's a lot to handle... Maybe try another time?")
	exit()
elif num != 0 and num < 1000:
	UserNumbSeqs = input("Do you want to keep going now that you know the number of sequences that will be analysed?")
	if UserNumbSeqs in yes:
		print("Great! Next step is performing a multiple alignment")
	elif UserNumbSeqs in no:
		exit()
	else:
		sys.stdout.write("please respond yes or no!")
		exit()

##Check for the species in the file
for line in UserFile:
	if line[0] == '>':
		splitline = line.split('[')
		NameSpecies = _splitline[1]
		#accessorID = AccNumbArrow[1:-1]
		print(NameSpecies)



UserFile.close()

##MULTIPLE SEQUENCE ALIGNMENT WITH CLUSTALO
print("ClustalO multiple alignment is being performed. Please be a little bit patient...")
Clustalo = 'clustalo -i {}.fasta -o {}MA.msf --outfmt msf --output-order tree-order'.format(TaxonID, TaxonID)
ClustaloOut = subprocess.check_output(Clustalo, shell=True)

#Check the ClustalO output
UserClustal = open('{}MA.msf'.format(TaxonID))
UserClustal_Contents = UserClustal.read()
print(UserClustal_Contents)
print("The number of aligned sequences with clustal is {}".format(num))

#IF THERE ARE MORE THAN 250 SEQUENCES, THE CODE REDUCES THE NUMBER TO 250 TO KEEP GOING WITH THE ANALYSIS
if num >= 250:
	##Get a consensus sequence to BLAST against
	print("There were more than 250 sequences. BLAST analysis will be performed to keep only the 250 most similar ones")
	Cons = 'cons -sequence {}MA.msf -outseq {}Cons.fasta'.format(TaxonID, TaxonID)
	ConsOut = subprocess.check_output(Cons, shell=True)
	print("A consensus sequence was generated")
	#BLAST the generated consensus seq against all the prot sequences found in the query
	#First generate a BLAST database to search against
	BlastDB = 'makeblastdb -in {}.fasta -dbtype prot'.format(TaxonID)
	#Perform the BLAST search
	BlastP = 'blastp -query {}Cons.fasta -db {}.fasta -outfmt 7 -out blastout.txt'.format(TaxonID, TaxonID)
	BlastDBOut = subprocess.check_output(BlastDB, shell=True)
	BlastPOut = subprocess.check_output(BlastP, shell=True)
	print("BLAST analysis performed")


	#Sort BLAST results by similarity using Pandas
	df = pd.read_csv('blastout.txt',sep='\t', skiprows=(0,1,2,3,4))
	df.columns = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "10", "11", "12"]
	pd.concat([pd.DataFrame(df.columns), df], ignore_index=True)
	##Column 3 contains the similarity values
	dfsorted = df.sort_values('3', ascending = False)
	print("Sorted by similarity")

	#This gets the accession number of the 250 more highly related proteins
	#Get accession numbers as a list
	listID = df['2'].to_list()
	ListID250 = listID[0:250]
	print("Acc numbers of the most similar sequences are shown")
	print(len(ListID250))

	#Create dictionary with the top 250 IDs
	FinalFASTA = open("SimilarSeqs.txt", mode="w")
	AI_DICT = {}
	for line in ListID250:
    		AI_DICT[line[:-1]] = 1

	UserFile = open('{}.fasta'.format(TaxonID), 'r')
	
	#Look for the sequences with an ID that is present in the top 250 IDs with higher similarity
	#Write them in a new file
	skip = 0
	for line in UserFile:
        	if line[0] == '>':
                	_splitline = line.split(' ')
                	AccNumbArrow = _splitline[0]
                	accessorID = AccNumbArrow[1:-1]
                	if accessorID in AI_DICT:
                        	FinalFASTA.write(line)
                        	skip = 0
                	else:
                        	skip = 1
        	else:
                	if not skip:
                        	FinalFASTA.write(line)
	
	FinalFASTA.close()
	
        ##Perform alignment again, this time with just the 250 most similar sequences
	print("ClustalO multiple alignment is being performed. Please be patient...")
	Clustalo = 'clustalo -i SimilarSeqs.txt -o SimilarSeqsMA.msf --outfmt msf -v --output-order tree-order'.format(TaxonID, TaxonID)
	ClustaloOut = subprocess.check_output(Clustalo, shell=True)

	##Plot similarity of the aligned sequences
	Plot = 'plotcon SimilarSeqsMA.msf -graph svg -goutfile plotcon'
	subprocess.check_output(Plot, shell=True)
	

else:
##Generate the plotcon with original dataset if the number of sequences was <250
	Plot = 'plotcon {}MA.msf -graph svg -goutfile plotcon'.format(TaxonID)
	subprocess.check_output(Plot, shell=True)

##The plot has been generated. Let the user know.
print("A similarity plot has been created and stored in your workspace! It is called plotcon.svg")
UserFile.close()


##PART 3: PROSITE MOTIF SEARCH

#Ask the user how many sequences to study for motifs
REQUESTED_LINES = input("How many sequences would you like to send for PROSITE motif analysis?:")
if num > 250:
	SelectedSeqs = 'awk "/^>/ {n++} n> %s {exit} {print}" SimilarSeqs.txt > ToPROSITE.fasta' % REQUESTED_LINES
	SelSeqs = subprocess.check_output(SelectedSeqs, shell=True)
else:
	SelectedSeqs = 'awk "/^>/ {n++} n> %s {exit} {print}" %s.fasta > ToPROSITE.fasta' % (REQUESTED_LINES, TaxonID) 
	SelSeqs = subprocess.check_output(SelectedSeqs, shell=True)

UserFile = open('ToPROSITE.fasta', 'r')

#Generate a single FASTA file for each sequence that will be analysed
#Necessary for individual output to patmatmotifs
outfile = []
for line in UserFile:
    if line.startswith(">"):
        if (outfile != []): outfile.close()
        AccIDarrow = line.strip().split(' ')[0]
        AccID = AccIDarrow[1:-1]
        filename = AccID+".fasta"
        outfile = open(filename,'w')
        outfile.write(line)
    else:
        outfile.write(line)
outfile.close()

#Run a search against the PROSITE motif database for the sequence number indicated by the user
#search for all files that end in ..fasta as indicated with the previous loop
for filename in os.listdir('.'):
	if filename.endswith("..fasta"):
		PROSITE = 'patmatmotifs -sequence {} -sprotein1 -sformat1 fasta -rname2 {} | cat'.format(filename, filename)
		PROSITE_OUT = subprocess.check_output(PROSITE, shell=True)
		continue
	else:
		continue

##Check if there is any hit (any motif found):
##The ID of any Interesting outputs (those showing hits) is indicated
for filename in os.listdir('.'):
	if filename.endswith(".patmatmotifs"):
		if re.search(r'HitCount: 0', filename):
			continue
		else:
			print("The sequence" + filename + "shows at least one motif from the PROSITE database")



