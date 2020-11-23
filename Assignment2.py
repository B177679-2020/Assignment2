#!/usr/bin/python3

##GETTING THE PROGRAM READY

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

##GETTING THE SEQUENCES REQUESTED BY THE USER

#Example query to get the birds G6P protein sequences in fasta format
#ESearchBirds = 'esearch -db protein -query "txid8782[Organism] AND glucose-6-phosphatase[Protein]" | efetch -db protein -format fasta > birds_G6P.fasta'
#ESearchBText = subprocess.check_output(ESearchBirds, shell=True)

##GET USER INPUT w. some error traps for typing errors
TaxonID = input("Please type in your Taxon ID (txidxxxx):")
if not TaxonID.startswith("txid"):
	print("Please, enter the Taxon ID in the correct format (txid followed by 4 or 5 numbers)")
	TaxonID = input("Please type in your Taxon ID again (remember, txidxxxx):")
#Get the name of the species to which the taxon refers to
TaxonInput = 'efetch -db taxonomy -id {} -format xml |   xtract -pattern Taxon     -element TaxId ScientificName GenbankCommonName Division | cat'.format(TaxonID)
txinputout = subprocess.check_output(TaxonInput, shell = True)
print("This is the name of your chosen taxon")
print(txinputout)


Protein = input("Please type in the protein of interest (ex.: phosphatase, kinase...):")

##Some options for the user
yes = {'yes','y', 'ye', 'YES', 'Yes'}
no = {'no','n', 'No', 'NO'}

##Search in the NCBI database

#Error traps to only allow for answers in yes or no, if not ask again
#Options for including predicted/partial sequences
#Options for include protein name variations in the query
i=0
while (i==0):
    PartialsPred = input("Do you want to include partial and predicted sequences of the protein in your search, if there are any? yes or no:")
    Variations = input("Do you want to include variations of your protein name? For example, if your input was Glucose-6-phosphatase, would you also like to see the isoforms X1, X2...? yes or no:")
    if PartialsPred in yes:
            if Variations in yes:
               i=1
               GeneralES = 'esearch -db protein -query "{}[Organism] AND {}*[Protein]" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)
            elif Variations in no:
               i=1
               GeneralES = 'esearch -db protein -query "{}[Organism] AND {}[Protein]" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)
            else:
               print("Please answer yes or no!") 
    elif PartialsPred in no:
            if Variations in yes:
               i=1
               GeneralES = 'esearch -db protein -query "{}[Organism] AND {}*[Protein] NOT partial NOT predicted" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)
            elif Variations in no:
               i=1
               GeneralES = 'esearch -db protein -query "{}[Organism] AND {}[Protein] NOT partial NOT predicted" | efetch -db protein -format fasta > {}.fasta'.format(TaxonID, Protein, TaxonID)
            else:
               print("Please answer yes or no!")
    else:
            print("Sorry, I did not understand that!")

GeneralESout = subprocess.check_output(GeneralES, shell=True)


#Check the user file can be open in Python and it has a reasonable number of sequences
UserFile = open('{}.fasta'.format(TaxonID))
UserFile_Contents = UserFile.read() 
print(UserFile_Contents)
#Count the number of sequences
num = len([1 for line in UserFile_Contents if line.startswith(">")])
print("The number of protein sequences for your query is {}".format(num))

#User options and input control for very large or empty files
#If sequences between 1 and 1000, give user the option to choose if run or not.
c = 0
while (c == 0):
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
			c = 1
		elif UserNumbSeqs in no:
			exit()
		else:
			sys.stdout.write("please respond yes or no!")
			

UserFile.close()

##MULTIPLE SEQUENCE ALIGNMENT WITH CLUSTALO
print("ClustalO multiple alignment is being performed. Please be patient...")
Clustalo = 'clustalo -i {}.fasta -o {}MA.msf --outfmt msf --output-order tree-order'.format(TaxonID, TaxonID)
ClustaloOut = subprocess.check_output(Clustalo, shell=True)

#Check the ClustalO output
#UserClustal = open('{}MA.msf'.format(TaxonID))
#UserClustal_Contents = UserClustal.read()
#print(UserClustal_Contents)
#print("The number of aligned sequences with clustal is {}".format(num))

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
	AI_DICT = {}
	for line in ListID250:
    		AI_DICT[line[:-1]] = 1

	UserFile = open('{}.fasta'.format(TaxonID), 'r')
	
	#Look for the sequences with an ID that is present in the top 250 IDs with higher similarity
	#Write them in a new file called SimilarSeqs.txt
	FinalFASTA = open("SimilarSeqs.txt", mode="w")
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
	##File obtained with 250 most similar sequences. Close it after finishing writing 
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

#Open the FASTA files with the sequences that will be searched for motifs
if num > 250:
	UserFile = open('SimilarSeqs.txt', 'r')
else:
	UserFile = open('{}.fasta'.format(TaxonID), 'r')

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

#Run a search against the PROSITE motif database for the sequences
#search for all files that end in ..fasta as indicated with the previous loop
for filename in os.listdir('.'):
	if filename.endswith("..fasta"):
		PROSITE = 'patmatmotifs -sequence {} -sprotein1 -sformat1 fasta -rname2 {} -auto | cat'.format(filename, filename)
		PROSITE_OUT = subprocess.check_output(PROSITE, shell=True)
		continue
	else:
		continue


##Check if there is any hit (any motif found):
##The ID of any Interesting outputs (those showing hits) is indicated
##All reports showing a hit are written in a new file PROSITEreport.txt
PAT = open('PROSITEreport.txt', 'w')
for filename in os.listdir('.'):
	if filename.endswith(".patmatmotifs"):
		report = open(filename, 'r').read()
		if re.search(r'HitCount: 0', report):
			continue
		else:
			print("The sequence" + filename + "shows at least one motif from the PROSITE database")
			PAT.write(report)

#Ask the user if they wish to display the reports on screen
print("A new report file has been created with information about the motifs found on the proteins. The file is called PROSITEreport.txt")
Display = input("Do you want to see the reports of the sequences that showed motifs from the PROSITE database? (yes/no):")

k = 0
while (k == 0):
	if Display in yes:
		print(PAT)
		k = 1
	if Display in no:
		print("That's OK! The file is saved in your directory, you can go back to it later.")
		k = 1
	else:
		print("Please respond yes or no")
PAT.close()

#Create a tidy summary output with the ID of sequences with motifs, and the name of the motifs
mfind = open('PROSITEreport.txt', 'r')
summary = open('PROSITESummary.txt', 'w')

#Loop that finds the ID of sequences and the motifs:
for line in mfind:
	if re.search(r'Sequence', line):
		print(line)
		summary.write(line)
	if re.search(r'Motif', line):
		print(line)
		summary.write(line)

mfind.close()
summary.close()

##SECTION 4: EXTRA ANALYSIS
#HYDROPHOBICITY PLOT
        ##Plot hydrophobicity of the aligned sequences
if num > 250:
	HPlot = 'pepwindowall SimilarSeqsMA.msf -graph svg -goutfile HydrophobicityPlot'
        subprocess.check_output(HPlot, shell=True)

else:
##Generate the hydrophobicity plot with original dataset if the number of sequences was <250
        HPlot = 'pepwindowall {}MA.msf -graph svg -goutfile HydrophobicityPlot'.format(TaxonID)
        subprocess.check_output(HPlot, shell=True)




##Get the output files of interest for the user out of the src directory so that they can be easily found
print("The analysis of your query is done. The directory that you created at the beginning contains the plotcon image and the full and summarized PROSITE reports. All the intermediate files are in the src directory, in case you need then for further analysis. Thanks!")

shutil.move('PROSITESummary.txt', '../PROSITESummary.txt')
shutil.move('PROSITEreport.txt', '../PROSITEreport.txt')
shutil.move('plotcon.svg', '../plotcon.svg')
shutil.move('HidrophobicityPlot.svg', '../HidrophobicityPlot.svg')
